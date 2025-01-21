#ifndef SPIHT_F_H
#define SPIHT_F_H

#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
void spiht_encode_ims(float *buffer, size_t height, size_t width, uint8_t **out_buffer, size_t *output_size, int num_stages);
void spiht_encode_file(char *filename, uint8_t **out_buffer, size_t *output_size, float bit_rate);
void spiht_decode_ims(uint8_t *buffer, size_t size, float *out_buffer, size_t height, size_t width, int num_bits);
void spiht_decode_file(uint8_t *in_buffer, size_t buffer_size, char *fileout);
void spiht_destroy_buffer(uint8_t *buffer);
#ifdef __cplusplus
}
#endif

#endif // SPIHT_F_H