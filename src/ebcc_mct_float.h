/* Custom Part-II MCT follows irreversible DC shift: its buffers are floats,
 * unlike upstream opj_mct_encode_custom's integer-buffer implementation. */
static OPJ_BOOL ebcc_mct_encode_float(OPJ_BYTE *coding, OPJ_SIZE_T n,
        OPJ_BYTE **data, OPJ_UINT32 components, OPJ_UINT32 is_signed)
{
    const OPJ_FLOAT32 *matrix = (const OPJ_FLOAT32 *)coding;
    OPJ_FLOAT32 *input = (OPJ_FLOAT32 *)opj_malloc(components * sizeof(*input));
    OPJ_ARG_NOT_USED(is_signed);
    if (!input) return OPJ_FALSE;
    for (OPJ_SIZE_T i = 0; i < n; ++i) {
        for (OPJ_UINT32 k = 0; k < components; ++k)
            input[k] = ((OPJ_FLOAT32 *)data[k])[i];
        for (OPJ_UINT32 j = 0; j < components; ++j) {
            double value = 0;
            for (OPJ_UINT32 k = 0; k < components; ++k)
                value += matrix[j * components + k] * (double)input[k];
            ((OPJ_FLOAT32 *)data[j])[i] = (OPJ_FLOAT32)value;
        }
    }
    opj_free(input);
    return OPJ_TRUE;
}
