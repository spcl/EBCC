/* The allocator reads OpenJPEG's tile encoder state, which is private API. */
#define OPJ_SKIP_POISON
#include "opj_includes.h"
#include "ebcc_base.h"
#include "ebcc_codec.h"
#include "ebcc_mct.h"
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include "log.h"

#define LEVELS 65535                /* samples are scaled to 16 bits */
#define PRECISION 16
#define CODEBLOCK_SIDE 64
#define PLACEHOLDER_PSNR 60.0f      /* required by OpenJPEG's parameter check; the allocator ignores it */

/* Decomposition depth */
#define MAX_RESOLUTIONS 8
#define MIN_RESOLUTIONS 4
#define MIN_COARSE_SIDE 9           /* the coarsest subband keeps at least this many samples per side */
#define PRECISION_MARGIN 2

/* Paper allocation weight; the residual enforces the final error bound. */
#define DEFAULT_KAPPA 2.25

static double read_kappa(void) {
    const char *value = getenv("EBCC_KAPPA");
    if (!value) return DEFAULT_KAPPA;
    char *end;
    errno = 0;
    double kappa = strtod(value, &end);
    if (end != value && !*end && errno != ERANGE && isfinite(kappa) && kappa > 0)
        return kappa;
    log_warn("Invalid EBCC_KAPPA; expected a finite positive number, using %.2f", DEFAULT_KAPPA);
    return DEFAULT_KAPPA;
}

/* ---- allocation ---------------------------------------------------------- */

typedef struct {
    double price;   /* bytes per unit of Tier-1 (synthesis-weighted) squared error */
    unsigned calls;
    int use_fixed_quantization;
} allocation_t;

/* Historical fixed-step experiment: scale QCD steps by 1.5 times the bound.
 * Tier-1 rounds the integer index before coding; every coding pass is retained. */
static OPJ_BOOL scale_quantization(opj_tccp_t *tccp, double scale) {
    if (!isfinite(scale) || scale <= 0 || tccp->qntsty != J2K_CCP_QNTSTY_SEQNT) return OPJ_FALSE;
    for (OPJ_UINT32 b = 0; b < 3 * tccp->numresolutions - 2; ++b) {
        opj_stepsize_t *step = &tccp->stepsizes[b];
        int shift;
        double fraction = frexp((1.0 + (double)step->mant / 2048) * scale, &shift);
        int exponent = step->expn - shift + 1;
        int mantissa = (int)floor((2 * fraction - 1.0) * 2048 + 0.5);
        if (mantissa == 2048) { mantissa = 0; --exponent; }
        if (exponent < 0 || exponent > 31) return OPJ_FALSE;
        step->expn = exponent;
        step->mant = (OPJ_INT32)mantissa;
    }
    return OPJ_TRUE;
}

/* Replaces OpenJPEG's rate allocation. Each codeblock keeps the prefix of its
 * coding passes that minimizes bytes plus the priced distortion the prefix
 * leaves; blocks are independent, so one scan of the Tier-1 pass statistics
 * suffices. The residual layer corrects whatever still exceeds the bound. */
static OPJ_BOOL allocate(opj_tcd_t *tcd, void *context) {
    allocation_t *allocation = context;
    opj_tcd_tile_t *tile = tcd->tcd_image->tiles;
    ++allocation->calls;
    if (!tile->numcomps || tcd->tcp->numlayers != 1) return OPJ_FALSE;
    tile->distolayer[0] = 0;
    for (OPJ_UINT32 c = 0; c < tile->numcomps; ++c) {
        opj_tcd_tilecomp_t *component = &tile->comps[c];
        for (OPJ_UINT32 r = 0; r < component->numresolutions; ++r) {
            opj_tcd_resolution_t *resolution = &component->resolutions[r];
            for (OPJ_UINT32 b = 0; b < resolution->numbands; ++b) {
                opj_tcd_band_t *band = &resolution->bands[b];
                for (OPJ_UINT32 p = 0; p < resolution->pw * resolution->ph; ++p) {
                    opj_tcd_precinct_t *precinct = &band->precincts[p];
                    for (OPJ_UINT32 k = 0; k < precinct->cw * precinct->ch; ++k) {
                        opj_tcd_cblk_enc_t *block = &precinct->cblks.enc[k];
                        /* distortiondec is the squared error removed up to a pass;
                         * the empty prefix costs nothing and removes nothing. */
                        double best = 0;
                        OPJ_UINT32 selected = allocation->use_fixed_quantization ? block->totalpasses : 0;
                        for (OPJ_UINT32 j = 0; !allocation->use_fixed_quantization && j < block->totalpasses; ++j) {
                            double cost = block->passes[j].rate - allocation->price * block->passes[j].distortiondec;
                            /* Ties go to the longer, more accurate prefix. */
                            if (cost <= best) { best = cost; selected = j + 1; }
                        }
                        opj_tcd_layer_t *layer = &block->layers[0];
                        block->numpassesinlayers = selected;
                        layer->numpasses = selected;
                        layer->len = selected ? block->passes[selected - 1].rate : 0;
                        layer->data = block->data;
                        layer->disto = selected ? block->passes[selected - 1].distortiondec : 0;
                        tile->distolayer[0] += layer->disto;
                        log_trace("Block: r=%u b=%u x0=%d y0=%d x1=%d y1=%d removed=%.17g rate=%u passes=%u/%u",
                            r, band->bandno, block->x0, block->y0, block->x1, block->y1, layer->disto, layer->len, selected, block->totalpasses);
                    }
                }
            }
        }
    }
    return OPJ_TRUE;
}

static void base_error(const char *message, void *context) {
    (void)context; log_error("OpenJPEG: %s", message);
}

/* ---- encoder setup ------------------------------------------------------- */

/* As deep as the geometry allows, but no deeper than the bound can use:
 * tight bounds leave the deepest levels nothing to discard. */
static int resolution_depth(size_t height, size_t width, double bound) {
    size_t side = height < width ? height : width;
    int geometry = 1;
    for (; side / 2 >= MIN_COARSE_SIDE && geometry < MAX_RESOLUTIONS; side /= 2) ++geometry;
    int precision = (int)floor(log2(fmax(1.0, bound * LEVELS))) + PRECISION_MARGIN;
    if (precision < MIN_RESOLUTIONS) precision = MIN_RESOLUTIONS;
    return precision < geometry ? precision : geometry;
}

typedef struct {
    uint8_t *data;
    size_t size, capacity, offset;
} memory_stream_t;

static OPJ_SIZE_T stream_write(void *data, OPJ_SIZE_T size, void *user) {
    memory_stream_t *s = user;
    if (size > SIZE_MAX - s->offset) return (OPJ_SIZE_T)-1;
    size_t end = s->offset + size;
    if (end > s->capacity) {
        size_t capacity = s->capacity ? s->capacity : 4096;
        while (capacity < end) capacity = capacity > SIZE_MAX / 2 ? end : capacity * 2;
        uint8_t *grown = realloc(s->data, capacity);
        if (!grown) return (OPJ_SIZE_T)-1;
        s->data = grown; s->capacity = capacity;
    }
    memcpy(s->data + s->offset, data, size);
    s->offset = end;
    if (s->size < end) s->size = end;
    return size;
}

static OPJ_SIZE_T stream_read(void *data, OPJ_SIZE_T size, void *user) {
    memory_stream_t *s = user;
    if (s->offset == s->size) return (OPJ_SIZE_T)-1;
    if (size > s->size - s->offset) size = s->size - s->offset;
    memcpy(data, s->data + s->offset, size);
    s->offset += size;
    return size;
}

static OPJ_BOOL stream_seek(OPJ_OFF_T offset, void *user) {
    memory_stream_t *s = user;
    if (offset < 0 || (OPJ_UINT64)offset > s->size) return OPJ_FALSE;
    s->offset = (size_t)offset;
    return OPJ_TRUE;
}

static OPJ_OFF_T stream_skip(OPJ_OFF_T distance, void *user) {
    memory_stream_t *s = user;
    if (distance < 0 || (OPJ_UINT64)distance > s->size - s->offset) return -1;
    s->offset += (size_t)distance;
    return distance;
}

static size_t encode_components(const float *data, size_t components, size_t height, size_t width,
        float lo, float hi, double bound, uint8_t **codestream, int use_fixed_quantization, float base_cr) {
    *codestream = NULL;
    size_t count = height * width, result = 0;
    double bound2 = bound * LEVELS * (bound * LEVELS);
    if (!isfinite(bound2) || bound2 <= 0) return 0;
    allocation_t allocation = {use_fixed_quantization || base_cr > 0 ? 0 : read_kappa() / bound2, 0, use_fixed_quantization};
    if (!isfinite(allocation.price)) return 0;

    opj_cparameters_t parameters;
    opj_set_default_encoder_parameters(&parameters);
    parameters.numresolution = resolution_depth(height, width, bound);
    parameters.cblockw_init = parameters.cblockh_init = CODEBLOCK_SIDE;
    parameters.irreversible = 1;
    parameters.tcp_numlayers = 1;
    if (base_cr > 0) {
        /* OpenJPEG measures rates against 16-bit samples, EBCC against float32. */
        parameters.cp_disto_alloc = 1;
        parameters.tcp_rates[0] = fmaxf(1.0f, base_cr / 2);
    } else {
        /* One quality layer whose contents are chosen by allocate(). */
        parameters.cp_fixed_quality = 1;
        parameters.tcp_distoratio[0] = PLACEHOLDER_PSNR;
    }
    /* Draw spatial grid offsets in [0, 7] to vary the base-error pattern.
     * JPEG2000 SIZ metadata carries them to the decoder. */
    OPJ_UINT32 x0 = (OPJ_UINT32)(rand() % 8);
    OPJ_UINT32 y0 = (OPJ_UINT32)(rand() % 8);
    parameters.image_offset_x0 = (int)x0;
    parameters.image_offset_y0 = (int)y0;

    opj_image_cmptparm_t component = {0};
    component.dx = component.dy = 1;
    component.w = (OPJ_UINT32)width; component.h = (OPJ_UINT32)height;
    component.prec = PRECISION;
    memory_stream_t memory = {0};
    opj_image_cmptparm_t descriptors[EBCC_MAX_MCT_COMPONENTS];
    if (!components || components > EBCC_MAX_MCT_COMPONENTS) return 0;
    for (size_t c = 0; c < components; ++c) descriptors[c] = component;
    opj_image_t *image = opj_image_create((OPJ_UINT32)components, descriptors,
        components == 1 ? OPJ_CLRSPC_GRAY : OPJ_CLRSPC_UNSPECIFIED);
    opj_codec_t *codec = opj_create_compress(OPJ_CODEC_J2K);
    opj_stream_t *stream = opj_stream_default_create(OPJ_FALSE);
    if (!image || !codec || !stream) goto done;
    opj_set_error_handler(codec, base_error, NULL);
    image->x0 = x0; image->x1 = x0 + (OPJ_UINT32)width;
    image->y0 = y0; image->y1 = y0 + (OPJ_UINT32)height;
    for (size_t c = 0; c < components; ++c)
    for (size_t i = 0; i < count; ++i)
        image->comps[c].data[i] = (OPJ_INT32)floor(((double)data[c * count + i] - lo) / ((double)hi - lo) * LEVELS + 0.5);
    opj_stream_set_user_data(stream, &memory, NULL);
    opj_stream_set_write_function(stream, stream_write);
    opj_stream_set_seek_function(stream, stream_seek);

    const char *disable_mct = getenv("EBCC_DISABLE_MCT");
    if (components > 1 && !(disable_mct && !strcmp(disable_mct, "1"))) {
        float matrix[EBCC_MAX_MCT_COMPONENTS * EBCC_MAX_MCT_COMPONENTS];
        OPJ_INT32 shifts[EBCC_MAX_MCT_COMPONENTS];
        for (size_t c = 0; c < components; ++c) shifts[c] = 1 << (PRECISION - 1);
        if (!ebcc_mct_matrix(components, matrix) ||
            !opj_set_MCT(&parameters, matrix, shifts, (OPJ_UINT32)components)) goto done;
    }
    if (!opj_setup_encoder(codec, &parameters, image)) goto done;
    opj_j2k_t *j2k = ((opj_codec_private_t *)codec)->m_codec;
    /* The library-version comment carries nothing the decoder needs. */
    opj_free(j2k->m_cp.comment);
    j2k->m_cp.comment = NULL;
    if (use_fixed_quantization) {
        for (size_t c = 0; c < components; ++c) {
            opj_tccp_t *tccp = &j2k->m_cp.tcps[0].tccps[c];
            tccp->ebcc_nearest = OPJ_TRUE;
            if (!scale_quantization(tccp, fmax(1.0, 1.5 * bound * LEVELS))) goto done;
        }
    }
    if (!opj_start_compress(codec, image, stream)) goto done;
    if (base_cr == 0) {
        j2k->m_tcd->ebcc_allocate = allocate;
        j2k->m_tcd->ebcc_context = &allocation;
    }
    if (!opj_encode(codec, stream) || (base_cr == 0 && allocation.calls != 1) ||
            !opj_end_compress(codec, stream)) goto done;
    log_trace("Base layer: resolutions=%d, grid_offset=(%u, %u), bytes=%zu",
        parameters.numresolution, x0, y0, memory.size);
    *codestream = memory.data;
    memory.data = NULL;
    result = memory.size;
done:
    opj_stream_destroy(stream); opj_destroy_codec(codec); opj_image_destroy(image);
    free(memory.data);
    opj_free(parameters.mct_data);
    return result;
}

size_t ebcc_base_encode(const float *data, size_t components, size_t height, size_t width,
        float lo, float hi, double bound, uint8_t **codestream, int use_fixed_quantization) {
    return encode_components(data, components, height, width, lo, hi, bound, codestream, use_fixed_quantization, 0);
}

size_t ebcc_base_encode_fixed_cr(const float *data, size_t components, size_t height, size_t width,
        float lo, float hi, float base_cr, uint8_t **codestream) {
    *codestream = NULL;
    if (!isfinite(base_cr) || base_cr <= 0) return 0;
    /* Use the geometry's full resolution depth, independently of an error bound. */
    return encode_components(data, components, height, width, lo, hi, 1, codestream, 0, base_cr);
}

int ebcc_base_decode_components(const uint8_t *codestream, size_t size, size_t components, float lo, float hi,
        float **field, size_t *height, size_t *width, int nearest) {
    memory_stream_t memory = {(uint8_t *)codestream, size, size, 0};
    opj_dparameters_t parameters;
    opj_set_default_decoder_parameters(&parameters);
    opj_codec_t *codec = opj_create_decompress(OPJ_CODEC_J2K);
    opj_stream_t *stream = opj_stream_default_create(OPJ_TRUE);
    opj_image_t *image = NULL;
    int ok = 0;
    if (!codec || !stream) goto done;
    opj_set_error_handler(codec, base_error, NULL);
    opj_stream_set_user_data(stream, &memory, NULL);
    opj_stream_set_user_data_length(stream, size);
    opj_stream_set_read_function(stream, stream_read);
    opj_stream_set_seek_function(stream, stream_seek);
    opj_stream_set_skip_function(stream, stream_skip);
    if (!opj_setup_decoder(codec, &parameters) || !opj_read_header(stream, codec, &image)) goto done;
    if (components > 1) {
        opj_j2k_t *j2k = ((opj_codec_private_t *)codec)->m_codec;
        const opj_tcp_t *tcp = j2k->m_specific_param.m_decoder.m_default_tcp;
        /* The codestream selects MCT independently of the encoder environment. */
        if (!tcp || (tcp->mct != 0 &&
                (tcp->mct != 2 || !tcp->m_mct_decoding_matrix))) goto done;
    }
    const opj_image_comp_t *samples = image->comps;
    if (!components || components > EBCC_MAX_MCT_COMPONENTS || image->numcomps != components) goto done;
    for (size_t c = 0; c < components; ++c) {
    samples = &image->comps[c];
    if (samples->sgnd || samples->prec != PRECISION || samples->dx != 1 || samples->dy != 1 ||
        !samples->w || samples->w > EBCC_MAX_INTERNAL_IMAGE_DIM ||
        !samples->h || samples->h > EBCC_MAX_INTERNAL_IMAGE_DIM ||
        samples->w != image->comps[0].w || samples->h != image->comps[0].h) goto done;
    }
    if (nearest) {
        opj_j2k_t *j2k = ((opj_codec_private_t *)codec)->m_codec;
        opj_tcp_t *defaults = j2k->m_specific_param.m_decoder.m_default_tcp;
        for (size_t c = 0; c < components; ++c) defaults->tccps[c].ebcc_nearest = OPJ_TRUE;
        for (OPJ_UINT32 t = 0; t < j2k->m_cp.tw * j2k->m_cp.th; ++t)
            for (size_t c = 0; c < components; ++c)
                j2k->m_cp.tcps[t].tccps[c].ebcc_nearest = OPJ_TRUE;
    }
    samples = image->comps;
    if (!opj_decode(codec, stream, image) || !opj_end_decompress(codec, stream)) goto done;
    size_t count = (size_t)samples->w * samples->h;
    if (!*field) *field = malloc(count * components * sizeof(**field));
    if (!*field) goto done;
    for (size_t c = 0; c < components; ++c)
    for (size_t i = 0; i < count; ++i)
        (*field)[c * count + i] = ((float)image->comps[c].data[i] / LEVELS) * (hi - lo) + lo;
    *height = samples->h; *width = samples->w;
    ok = 1;
done:
    opj_stream_destroy(stream); opj_destroy_codec(codec); opj_image_destroy(image);
    return ok;
}

int ebcc_base_decode(const uint8_t *codestream, size_t size, float lo, float hi,
        float **field, size_t *height, size_t *width) {
    return ebcc_base_decode_components(codestream, size, 1, lo, hi, field, height, width, 0);
}
