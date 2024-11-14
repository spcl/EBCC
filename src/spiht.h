#ifndef SPIHT_F_H
#define SPIHT_F_H

#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
void spiht_encode(float *buffer, size_t height, size_t width, uint8_t **out_buffer, size_t *output_size, int num_stages, float* max_val);
void spiht_decode(uint8_t *buffer, size_t size, float *out_buffer, size_t height, size_t width, int num_bits, float max_val);
void spiht_destroy_buffer(uint8_t *buffer);
#ifdef __cplusplus
}
#endif

#endif // SPIHT_F_H