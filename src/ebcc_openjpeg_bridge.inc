/* Compiled inside the build-directory copy of OpenJPEG's j2k.c. These helpers
 * are private to ebcc_j2k.c, which fixes all coding parameters except target.
 * No process-global or thread-local cache is used.
 */
OPJ_BOOL ebcc_opj_finish_prepared(opj_codec_t *codec, opj_stream_t *stream)
{
    opj_codec_private_t *c = (opj_codec_private_t *)codec;
    opj_j2k_t *j2k = (opj_j2k_t *)c->m_codec;
    opj_tcd_t *prepared = j2k->m_tcd;
    OPJ_BOOL ok;

    /* Finish normally, but retain ownership of Tier-1 until codec destruction.
     * end_compress otherwise destroys it alongside the output/header buffers.
     */
    j2k->m_tcd = NULL;
    ok = opj_end_compress(codec, stream);
    j2k->m_tcd = prepared;
    if (ok) {
        opj_tcd_tilecomp_t *tile = prepared->tcd_image->tiles->comps;
        opj_image_comp_t *image = j2k->m_private_image->comps;
        /* Rate allocation and packetization only need code-blocks and image
         * metadata. Release transformed pixels, retaining neither input copy.
         */
        if (tile->ownsData) opj_aligned_free(tile->data);
        tile->data = NULL;
        tile->ownsData = OPJ_FALSE;
        opj_image_data_free(image->data);
        image->data = NULL;
    }
    return ok;
}

OPJ_BOOL ebcc_opj_encode_prepared(opj_codec_t *cached, opj_codec_t *fresh,
                                 opj_stream_t *stream)
{
    opj_codec_private_t *c = (opj_codec_private_t *)cached;
    opj_codec_private_t *f = (opj_codec_private_t *)fresh;
    opj_j2k_t *old = (opj_j2k_t *)c->m_codec;
    opj_j2k_t *next = (opj_j2k_t *)f->m_codec;
    opj_tcd_t *tile = old->m_tcd;
    opj_tcd_t *unused = next->m_tcd;
    opj_cp_t *saved_cp;
    opj_tcp_t *saved_tcp;
    opj_image_t *saved_image;
    OPJ_BOOL ok;

    /* The session API supplies these invariants; check the bridge's essential
     * assumptions even in release builds before accessing retained state.
     */
    if (!tile || old->m_cp.tw != 1 || old->m_cp.th != 1 ||
        next->m_cp.tw != 1 || next->m_cp.th != 1 ||
        next->m_cp.tcps[0].numlayers != 1 ||
        next->m_private_image->numcomps != 1 ||
        old->m_private_image->x1 != next->m_private_image->x1 ||
        old->m_private_image->y1 != next->m_private_image->y1)
        return OPJ_FALSE;

    saved_cp = tile->cp;
    saved_tcp = tile->tcp;
    saved_image = tile->image;
    /* Borrow the session-owned tile for this emission. Fresh setup handles
     * headers and CR-to-byte-budget conversion. Restore ownership on all exits.
     * Layer allocation resets pass counts; Tier-2 resets packet tag trees.
     */
    next->m_tcd = tile;
    tile->cp = &next->m_cp;
    tile->tcp = &next->m_cp.tcps[0];
    tile->image = next->m_private_image;
    tile->cur_totnum_tp = next->m_cp.tcps[0].m_nb_tile_parts;
    tile->tcd_image->tiles->packno = 0;
    tile->ebcc_reuse_tier1 = OPJ_TRUE;
    ok = opj_j2k_post_write_tile(next, (opj_stream_private_t *)stream,
                               &f->m_event_mgr);
    tile->ebcc_reuse_tier1 = OPJ_FALSE;
    tile->cp = saved_cp;
    tile->tcp = saved_tcp;
    tile->image = saved_image;
    next->m_tcd = unused;
    return ok;
}
