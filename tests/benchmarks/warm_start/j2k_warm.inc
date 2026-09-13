/* Investigation only: appended to a temporary copy of OpenJPEG's j2k.c.
 * Limited to the benchmark's identical-input, single-tile, single-component,
 * single-layer, irreversible codestreams. This is not a public OpenJPEG API.
 */
_Thread_local opj_tcd_t *ebcc_warm_tile = NULL;

OPJ_BOOL ebcc_finish_keep(opj_codec_t *codec, opj_stream_t *stream)
{
    opj_codec_private_t *c = (opj_codec_private_t *)codec;
    return opj_j2k_write_eoc((opj_j2k_t *)c->m_codec,
                            (opj_stream_private_t *)stream, &c->m_event_mgr);
}

OPJ_BOOL ebcc_encode_again(opj_codec_t *cached, opj_codec_t *fresh,
                          opj_stream_t *stream)
{
    opj_codec_private_t *c = (opj_codec_private_t *)cached;
    opj_codec_private_t *f = (opj_codec_private_t *)fresh;
    opj_j2k_t *old = (opj_j2k_t *)c->m_codec;
    opj_j2k_t *next = (opj_j2k_t *)f->m_codec;
    opj_tcd_t *tile = old->m_tcd;
    opj_tcd_t *unused = next->m_tcd;
    opj_cp_t *saved_cp = tile->cp;
    opj_tcp_t *saved_tcp = tile->tcp;
    opj_image_t *saved_image = tile->image;
    OPJ_BOOL ok;

    assert(old->m_cp.tw * old->m_cp.th == 1);
    assert(next->m_cp.tw * next->m_cp.th == 1);
    assert(next->m_cp.tcps[0].numlayers == 1);
    assert(next->m_private_image->numcomps == 1);

    /* Fresh setup writes headers and converts CR to a tile-byte budget.
     * Reuse the completed Tier-1 code-blocks, passes and distortion values.
     * Rate allocation resets layer pass counts; Tier-2 resets tag trees.
     */
    next->m_tcd = tile;
    tile->cp = &next->m_cp;
    tile->tcp = &next->m_cp.tcps[0];
    tile->image = next->m_private_image;
    tile->cur_totnum_tp = next->m_cp.tcps[0].m_nb_tile_parts;
    tile->tcd_image->tiles->packno = 0;
    ebcc_warm_tile = tile;
    ok = opj_j2k_post_write_tile(next, (opj_stream_private_t *)stream,
                               &f->m_event_mgr);
    ebcc_warm_tile = NULL;
    tile->cp = saved_cp;
    tile->tcp = saved_tcp;
    tile->image = saved_image;
    next->m_tcd = unused;
    return ok;
}
