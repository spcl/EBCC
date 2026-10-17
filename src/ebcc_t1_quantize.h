/* Included in the build-local OpenJPEG Tier-1 implementation. */
#include "ebcc_quantize.h"

static OPJ_INT32 ebcc_t1_quantize(OPJ_FLOAT32 coefficient, OPJ_FLOAT32 step,
        OPJ_BOOL nearest) {
    if (!nearest)
        return (OPJ_INT32)opj_lrintf((coefficient / step) * (1 << T1_NMSEDEC_FRACBITS));
    double index = ebcc_nearest((double)coefficient / step);
    const OPJ_INT32 limit = INT_MAX / (1 << T1_NMSEDEC_FRACBITS);
    /* Avoid integer overflow; the outer verified residual repairs saturation. */
    index = fmax(-limit, fmin(limit, index));
    return (OPJ_INT32)index * (1 << T1_NMSEDEC_FRACBITS);
}

/* Full Tier-1 decoding yields sign(q)*(2*abs(q)+1). The usual inverse
 * quantizer multiplies by half a step; remove that nonzero midpoint first. */
static OPJ_INT32 ebcc_t1_reconstruct(OPJ_INT32 value) {
    return value - (value > 0) + (value < 0);
}
