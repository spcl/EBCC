/* Re-implementation of René Puchinger's SPIHT program: https://github.com/gilson27/imshrinker */
#include "spiht_re.h"
#include "dwt.h"
#include "ml.h"
#include "bitio.h"

typedef enum {
    SPIHT_ENCODE,
    SPIHT_DECODE,
} SPIHTMode;

struct SPIHTState {
    DWTData* dwt_data;
    MaskedList* lip; /* List of insignificant pixels */
    MaskedList* lsp; /* List of significant pixels */
    MaskedList* lis; /* List of insignificant sets */
    int_t step; /* Quantization step */
    SPIHTMode mode;
};
typedef struct SPIHTState SPIHTState;

void free_spiht_state(SPIHTState* state) {
    ml_free(state->lip);
    ml_free(state->lsp);
    ml_free(state->lis);
    free(state);
}

SPIHTState* spiht_encode_init(BitIOStreamState* bio, DWTData* dwt_data, size_t trunc_bits) {
    SPIHTState* state = (SPIHTState*) malloc(sizeof(SPIHTState));
    elem_t max = 2.0; /* step = log2(max) is at least 0 */
    size_t size_x = dwt_data->size_x;
    size_t size_y = dwt_data->size_y;
    size_t extra_x = dwt_data->extra_x;
    size_t extra_y = dwt_data->extra_y;
    size_t stride = dwt_data->stride;
    size_t num_stages = dwt_data->num_stages;
    size_t first_stage_size_x = (size_x + extra_x) / (1 << num_stages);
    size_t first_stage_size_y = (size_y + extra_y) / (1 << num_stages);
    size_t image_size = (size_x + extra_x) * (size_y + extra_y);
    size_t capacity = (first_stage_size_x * first_stage_size_y) + trunc_bits;
    if (capacity > image_size * 8) {
        /* trunc_bits can go unlimited (2^28), cap it with image data size */
        capacity = image_size * 8;
    }
    assert (capacity <= MAXINT); /* such that x, y, set flag can be encoded into 64 bit int*/
    assert (dwt_data->type == DWTCOEFF);
    assert (bio->mode == BITWRITE);
    state->lip = ml_init(capacity);
    state->lsp = ml_init(capacity);
    state->lis = ml_init(capacity);
    state->mode = SPIHT_ENCODE;
    for (size_t i = 0; i < image_size; i++) {
        elem_t absval = fabs(dwt_data->data[i]);
        if (absval > max) {
            max = absval;
        }
    }
    state->step = (int_t) floor(log(max) / log(2.0));
    assert(state->step <= MAXSTEPS);

    bitio_put_bits(bio, (uint64_t) state->step, 8);

    for (size_t y = 0; y < first_stage_size_y; y++) {
        for (size_t x = 0; x < first_stage_size_x; x++) {
            int_t pix_item = x + y * stride;
            ml_push(state->lip, pix_item);
            if ((x % 2 != 0) || (y % 2 != 0)) {
                /* set_item = (pix_item, A)  
                (pix_item, A) -> pix_item + 1
                (pix_item, B) -> - (pix_item + 1)  */
                int_t set_item = pix_item + 1; 
                ml_push(state->lis, set_item);
            }
        }
    }
    return state;
}

SPIHTState* spiht_decode_init(BitIOStreamState* bio, DWTData* dwt_data, size_t trunc_bits) {
    SPIHTState* state = (SPIHTState*) malloc(sizeof(SPIHTState));
    size_t size_x = dwt_data->size_x;
    size_t size_y = dwt_data->size_y;
    size_t extra_x = dwt_data->extra_x;
    size_t extra_y = dwt_data->extra_y;
    size_t stride = dwt_data->stride;
    size_t num_stages = dwt_data->num_stages;
    size_t first_stage_size_x = (size_x + extra_x) / (1 << num_stages);
    size_t first_stage_size_y = (size_y + extra_y) / (1 << num_stages);
    size_t image_size = (size_x + extra_x) * (size_y + extra_y);
    size_t capacity = (first_stage_size_x * first_stage_size_y) + trunc_bits;
    if (capacity > image_size * 8) capacity = image_size * 8;
    assert (capacity <= MAXINT); /* such that x, y, set flag can be encoded into 64 bit int*/
    assert (dwt_data->type == DWTCOEFF);
    assert (bio->mode == BITREAD);
    state->lip = ml_init(capacity);
    state->lsp = ml_init(capacity);
    state->lis = ml_init(capacity);
    state->mode = SPIHT_DECODE;
    memset(dwt_data->data, 0, image_size * sizeof(elem_t)); /* Clear DWT Coeff data */

    state->step = (int_t) bitio_get_bits(bio, 8);
    assert(state->step <= MAXSTEPS);

    for (size_t y = 0; y < first_stage_size_y; y++) {
        for (size_t x = 0; x < first_stage_size_x; x++) {
            int_t pix_item = x + y * stride;
            ml_push(state->lip, pix_item);
            if ((x % 2 != 0) || (y % 2 != 0)) {
                int_t set_item = pix_item + 1;
                ml_push(state->lis, set_item);
            }
        }
    }
    return state;
}

static inline bool_t is_significant_pixel(int_t step, elem_t val) {
    /* TODO: preprocess val into integer */
    /* Version 1: return (((int_t) fabs(val)) >= (1 << step)); */
    /* Version 2: return (( fabs(val)) >= (double) (1 << step)); */
    /* Version 3: (used in imshrinker)*/
    return abs((int_t) val) >= (1 << step);
}

static inline void get_successor(int_t x, int_t y, DWTData* dwt_data, int_t* sx, int_t* sy) {
    size_t size_x_pad = dwt_data->size_x + dwt_data->extra_x;
    size_t size_y_pad = dwt_data->size_y + dwt_data->extra_y;
    int_t lx = size_x_pad / (1 << dwt_data->num_stages);
    int_t ly = size_y_pad / (1 << dwt_data->num_stages);
    int_t sx_, sy_;
    if (x < lx && y < ly) {
        if (x % 2 == 1) {
            sx_ = x + lx - 1;
        } else {
            sx_ = x;
        }
        if (y % 2 == 1) {
            sy_ = y + ly - 1;
        } else {
            sy_ = y;
        }
        if (sx_ == x && sy_ == y) {
            sx_ = -1;
            sy_ = -1;
        }
    } else {
        sx_ = 2 * x;
        sy_ = 2 * y;
        if (sx_ >= size_x_pad || sy_ >= size_y_pad) {
            sx_ = -1;
            sy_ = -1;
        }
    }
    *sx = sx_;
    *sy = sy_;
}

static inline bool_t is_significant_set_A(int_t step, DWTData* dwt_data, int_t pix_item, size_t count) {
    size_t stride = dwt_data->stride;
    int_t sx, sy, x, y;
    if (count > 1 && is_significant_pixel(step, dwt_data->data[pix_item])) {
        return 1;
    }
    x = pix_item % stride;
    y = pix_item / stride;
    get_successor(x, y, dwt_data, &sx, &sy);
    if (sx == -1 || sy == -1) {
        return 0;
    }
    if (is_significant_set_A(step, dwt_data, sx + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_A(step, dwt_data, sx + 1 + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_A(step, dwt_data, sx + (sy + 1) * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_A(step, dwt_data, sx + 1 + (sy + 1) * stride, count + 1)) {
        return 1;
    }
    return 0;
}

static inline bool_t is_significant_set_B(int_t step, DWTData* dwt_data, int_t pix_item, size_t count) {
    size_t stride = dwt_data->stride;
    int_t sx, sy, x, y;
    if (count > 2 && is_significant_pixel(step, dwt_data->data[pix_item])) {
        return 1;
    }
    x = pix_item % stride;
    y = pix_item / stride;
    get_successor(x, y, dwt_data, &sx, &sy);
    if (sx == -1 || sy == -1) {
        return 0;
    }
    if (is_significant_set_B(step, dwt_data, sx + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_B(step, dwt_data, sx + 1 + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_B(step, dwt_data, sx + (sy + 1) * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_B(step, dwt_data, sx + 1 + (sy + 1) * stride, count + 1)) {
        return 1;
    }
    return 0;
}

void spiht_encode_process(size_t bits, BitIOStreamState* bio, DWTData* dwt_data, SPIHTState* state) {
    size_t bit_cnt = 0;
    MaskedList* lip = state->lip;
    MaskedList* lsp = state->lsp;
    MaskedList* lis = state->lis;
    size_t stride = dwt_data->stride;
    assert(dwt_data->type == DWTCOEFF);
    assert(state->mode == SPIHT_ENCODE);
    /* TODO: check num bits only in the end of each loop and in put_bit */
    for (int_t step = state->step; step >= 0; step--) {
        /* Sorting pass */
        /* First process LIP */
        for (size_t i = 0; i < lip->curr; i++) {
            int_t pix_item = ml_get(lip, i);
            elem_t val = dwt_data->data[pix_item];
            bool_t sig = is_significant_pixel(step, val);
            bitio_put_bit(bio, (uint8_t) sig);
            if (++bit_cnt > bits) return;
            if (sig) {
                ml_push(lsp, pix_item);
                /* Encode the sign bit */
                bitio_put_bit(bio, (uint8_t) ((val > 0) ? 0 : 1));
                if (++bit_cnt > bits) return;
                ml_remove(lip, i);
            }
        }
        ml_consolidate(lip);

        /* now process LIS */
        for (size_t i = 0; i < lis->curr; i++) {
            int_t set_item = ml_get(lis, i);
#ifdef DEBUG
            assert(set_item != 0);
#endif
            bool_t is_set_A = set_item > 0;
            int_t sx, sy;
            if (is_set_A) {
                int_t pix_item = set_item - 1;
                int_t x = pix_item % stride;
                int_t y = pix_item / stride;
                bool_t sig = is_significant_set_A(step, dwt_data, pix_item, 1);
                bitio_put_bit(bio, (uint8_t) sig);
                if (++bit_cnt > bits) return;
                if (sig) {
                    get_successor(x, y, dwt_data, &sx, &sy);
                    /* process the four offsprings */
                    
                    for (int_t dy = 0; dy < 2; dy++) {
                        for (int_t dx = 0; dx < 2; dx++) {
                            int_t pix_item = sx + dx + (sy + dy) * stride;
                            elem_t val = dwt_data->data[pix_item];
                            sig = is_significant_pixel(step, val);
                            bitio_put_bit(bio, (uint8_t) sig);
                            if (++bit_cnt > bits) return;
                            if (sig) {
                                ml_push(lsp, pix_item);
                                bitio_put_bit(bio, (uint8_t) ((val > 0) ? 0 : 1));
                                if (++bit_cnt > bits) return;
                            } else {
                                ml_push(lip, pix_item);
                            }
                        }
                    }

                    /* test if L(i, j) != 0 */
                    get_successor(sx, sy, dwt_data, &sx, &sy);
                    if (sx != -1) {
                        /* push (x, y, B) in lis */
                        int_t set_item_B = - (x + y * stride + 1);
                        ml_push(lis, set_item_B);
                    }

                    ml_remove(lis, i);
                }
            } else {
                /* Set B */
                int_t pix_item = - set_item - 1;
                int_t x = pix_item % stride;
                int_t y = pix_item / stride;
                bool_t sig = is_significant_set_B(step, dwt_data, pix_item, 1);
                bitio_put_bit(bio, (uint8_t) sig);
                if (++bit_cnt > bits) return;
                if (sig) {
                    int_t set_item_A;
                    get_successor(x, y, dwt_data, &sx, &sy);
                    set_item_A = sx + sy * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + 1 + sy * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + (sy + 1) * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + 1 + (sy + 1) * stride + 1;
                    ml_push(lis, set_item_A);
                    ml_remove(lis, i);
                }
            }
        }
        ml_consolidate(lis);

        /* Refinement pass */
        for (size_t i = 0; i < lsp->curr; i++) {
            elem_t val = dwt_data->data[ml_get(lsp, i)];
            if (is_significant_pixel(step+1, val)) {
                bitio_put_bit(bio, (uint8_t) (( ((int_t) abs((int_t) val)) >> step) & 1));
                if (++bit_cnt > bits) return;
            }
        }
        
    }
}

void spiht_decode_process(size_t bits, BitIOStreamState* bio, DWTData* dwt_data, SPIHTState* state) {
    size_t bit_cnt = 0;
    MaskedList* lip = state->lip;
    MaskedList* lsp = state->lsp;
    MaskedList* lis = state->lis;
    size_t stride = dwt_data->stride;
    assert(dwt_data->type == DWTCOEFF);
    assert(state->mode == SPIHT_DECODE);
    /* TODO: check num bits only in the end of each loop and in get_bit */
    for (int_t step = state->step; step >= 0; step--) {
        /* Sorting pass */
        /* First process LIP */
        for (size_t i = 0; i < lip->curr; i++) {
            int_t pix_item = ml_get(lip, i);
            bool_t sig = bitio_get_bit(bio);
            if (++bit_cnt > bits) return;
            if (sig) {
                ml_push(lsp, pix_item);
                /* Decode the sign bit */
                dwt_data->data[pix_item] = (elem_t) ((bitio_get_bit(bio)? -1 : 1) * (1 << step));
                if (++bit_cnt > bits) return;
                ml_remove(lip, i);
            }
        }
        ml_consolidate(lip);

        /* now process LIS */
        for (size_t i = 0; i < lis->curr; i++) {
            int_t set_item = ml_get(lis, i);
#ifdef DEBUG
            assert(set_item != 0);
#endif
            bool_t is_set_A = set_item > 0;
            int_t sx, sy;
            if (is_set_A) {
                int_t pix_item = set_item - 1;
                int_t x = pix_item % stride;
                int_t y = pix_item / stride;
                bool_t sig = bitio_get_bit(bio);
                if (++bit_cnt > bits) return;
                if (sig) {
                    get_successor(x, y, dwt_data, &sx, &sy);
                    /* process the four offsprings */
                    
                    for (int_t dy = 0; dy < 2; dy++) {
                        for (int_t dx = 0; dx < 2; dx++) {
                            int_t pix_item = sx + dx + (sy + dy) * stride;
                            sig = bitio_get_bit(bio);
                            if (++bit_cnt > bits) return;
                            if (sig) {
                                ml_push(lsp, pix_item);
                                dwt_data->data[pix_item] = (elem_t) ((bitio_get_bit(bio)? -1 : 1) * (1 << step));
                                if (++bit_cnt > bits) return;
                            } else {
                                ml_push(lip, pix_item);
                            }
                        }
                    }

                    /* test if L(i, j) != 0 */
                    get_successor(sx, sy, dwt_data, &sx, &sy);
                    if (sx != -1) {
                        /* push (x, y, B) in lis */
                        int_t set_item_B = - (x + y * stride + 1);
                        ml_push(lis, set_item_B);
                    }

                    ml_remove(lis, i);
                }
            } else {
                /* Set B */
                int_t pix_item = - set_item - 1;
                int_t x = pix_item % stride;
                int_t y = pix_item / stride;
                bool_t sig = bitio_get_bit(bio);
                if (++bit_cnt > bits) return;
                if (sig) {
                    int_t set_item_A;
                    get_successor(x, y, dwt_data, &sx, &sy);
                    set_item_A = sx + sy * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + 1 + sy * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + (sy + 1) * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + 1 + (sy + 1) * stride + 1;
                    ml_push(lis, set_item_A);
                    ml_remove(lis, i);
                }
            }
        }
        ml_consolidate(lis);

        /* Refinement pass */
        for (size_t i = 0; i < lsp->curr; i++) {
            int_t pix_item = ml_get(lsp, i);
            elem_t val = dwt_data->data[pix_item];
            int_t val_int = (int_t) val;
            if (is_significant_pixel(step+1, val)) {
                if (bitio_get_bit(bio)) {
                    if (val_int >= 0)
                        dwt_data->data[pix_item] = (elem_t) (val_int | (1 << step));
                    else
                        dwt_data->data[pix_item] = (elem_t) (-( (-val_int) | (1 << step)));
                } else {
                    dwt_data->data[pix_item] = (elem_t) (val_int & (~(1 << step)));
                }
                if (++bit_cnt > bits) return;
            }
        }
    }
}

void spiht_encode(elem_t *buffer, size_t height, size_t width, uint8_t** buffer_out, size_t* output_size, size_t trunc_bits, size_t num_stages) {
    size_t buffer_size = (trunc_bits == 0) ? height * width * sizeof(elem_t) : trunc_bits / sizeof(uint8_t) + 1 ;
    if (buffer_size < 16) buffer_size = 16;
    DWTData* dwt_data = load_image(buffer, height, width, num_stages);
    *buffer_out = (uint8_t*) calloc(buffer_size, sizeof(uint8_t));
    BitIOStreamState* bio = bitio_init(*buffer_out, buffer_size, BITWRITE);
    assert(num_stages <= 32 && num_stages >= 1 );
    size_t size_x = dwt_data->size_x;
    size_t size_y = dwt_data->size_y;
    size_t extra_x = dwt_data->extra_x;
    size_t extra_y = dwt_data->extra_y;
    assert(size_x <= 2047 && size_x >= 1);
    assert(size_y <= 2047 && size_y >= 1);
    assert(extra_x <= 511 && extra_x >= 0);
    assert(extra_y <= 511 && extra_y >= 0);
    /* write IMS header */
    bitio_put_bits(bio, (uint64_t) 'I', 8);
    bitio_put_bits(bio, (uint64_t) 'M', 8);
    bitio_put_bits(bio, (uint64_t) 'S', 8);
    bitio_put_bits(bio, (uint64_t) num_stages, 6);
    bitio_put_bits(bio, (uint64_t) size_x, 12);
    bitio_put_bits(bio, (uint64_t) size_y, 12);
    bitio_put_bits(bio, (uint64_t) extra_x, 10);
    bitio_put_bits(bio, (uint64_t) extra_y, 10);
    bitio_put_bit(bio, (uint8_t) 0); /* is_color = False */
    size_t offset = 128; /* 128 (114) for metadata bits */
    size_t bits0 = (trunc_bits == 0) ? (1 << 28) : trunc_bits + offset;
    assert(bits0 > 0);
    bitio_put_bits(bio, (uint64_t) bits0, 29);
    elem_t dc0 = sub_dc(dwt_data);
    assert(dc0 >= 0 && dc0 <= MAXELEM);
    uint8_t dc0_uint = (uint8_t) dc0;
    bitio_put_bits(bio, (uint64_t) dc0_uint, 8);
    /* TODO: dump DWT coeffs and compare */
    dwt2full(dwt_data);
    normalize(dwt_data);
    SPIHTState* state = spiht_encode_init(bio, dwt_data, bits0 - offset);
    spiht_encode_process(bits0 - offset, bio, dwt_data, state);
    bitio_flush(bio);
    *output_size = bitio_get_size(bio);
    bitio_free(bio);
    free_image(dwt_data);
    free_spiht_state(state);
}

void spiht_decode(uint8_t* buffer_in, size_t input_size, elem_t* buffer_out, size_t height, size_t width, size_t num_bits) {
    BitIOStreamState* bio = bitio_init(buffer_in, input_size, BITREAD);

    uint8_t magic[3];
    magic[0] = (uint8_t) bitio_get_bits(bio, 8);
    magic[1] = (uint8_t) bitio_get_bits(bio, 8);
    magic[2] = (uint8_t) bitio_get_bits(bio, 8);
    assert(magic[0] == 'I' && magic[1] == 'M' && magic[2] == 'S');
    size_t num_stages = (size_t) bitio_get_bits(bio, 6);
    size_t size_x = (size_t) bitio_get_bits(bio, 12);
    size_t size_y = (size_t) bitio_get_bits(bio, 12);
    size_t extra_x = (size_t) bitio_get_bits(bio, 10);
    size_t extra_y = (size_t) bitio_get_bits(bio, 10);
    bool_t is_color = (bool_t) bitio_get_bit(bio);
    elem_t scale = MAXELEM;

    size_t bits0 = (size_t) bitio_get_bits(bio, 29);
    size_t offset = 128; /* 128 (114) for metadata bits */
    if (num_bits > bits0) {
        fprintf(stderr, "Warning: num_bits: %d > bits0: %d, taking bits0 instead\n", num_bits, bits0);
        fprintf(stderr, "Warning: num_stages: %d, size_x: %d, size_y: %d, extra_x: %d, extra_y: %d\n", num_stages, size_x, size_y, extra_x, extra_y);
        if (bits0) num_bits = bits0; else assert(0);
    }
    num_bits -= offset;
    assert(num_bits > 0);
    uint8_t dc0_uint = (uint8_t) bitio_get_bits(bio, 8);
    elem_t dc0 = (elem_t) dc0_uint;

    DWTData* dwt_data = empty_dwtcoeff(num_stages, size_x, size_y, extra_x, extra_y);
    SPIHTState* state = spiht_decode_init(bio, dwt_data, num_bits);
    spiht_decode_process(num_bits, bio, dwt_data, state);
    idwt2full(dwt_data);
    add_dc(dwt_data, dc0);

    size_t stride = dwt_data->stride;
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            buffer_out[x + y * width] = (dwt_data->data[x + y * stride]) / scale;
        }
    }
    bitio_free(bio);
    free_image(dwt_data);
    free_spiht_state(state);
}