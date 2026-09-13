#ifndef EBCC_J2K_H
#define EBCC_J2K_H

#include "openjpeg.h"
#include <stdint.h>

#define EBCC_J2K_SAMPLE_BITS 16

/* A session owns an immutable copy of one unsigned 16-bit grayscale tile and
 * its irreversible Tier-1 encoding. Only the PSNR/rate target can change.
 * Independent sessions can be used concurrently; a session is single-owner.
 */
typedef struct ebcc_j2k_session ebcc_j2k_session_t;

ebcc_j2k_session_t *ebcc_j2k_session_create(const uint16_t *data,
                                          OPJ_UINT32 width, OPJ_UINT32 height);
void ebcc_j2k_session_destroy(ebcc_j2k_session_t *session);

/* Write one complete J2K codestream. rate is the OpenJPEG (uint16) ratio.
 * Return 1 on success, -1 on a new error, or 0 if a prior error disabled reuse.
 * After an error the caller may retry with its original input and cold encoder.
 */
int ebcc_j2k_session_encode(ebcc_j2k_session_t *session, float psnr, float rate,
                           opj_stream_write_fn write, void *output);

#endif
