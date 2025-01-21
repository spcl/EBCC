#ifndef SPIHT_F_H
#define SPIHT_F_H

#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
void spiht_encode_ims(float *buffer, size_t height, size_t width, uint8_t **out_buffer, size_t *output_size, size_t trunc_bits, int num_stages);
void spiht_decode_ims(uint8_t *buffer, size_t size, float *out_buffer, size_t height, size_t width, int num_bits);
void spiht_destroy_buffer(uint8_t *buffer);
#ifdef __cplusplus
}
#endif

#endif // SPIHT_F_H