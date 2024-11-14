#ifndef SPIHT_H
#define SPIHT_H

#include <cstddef>
#include <cstdint>

extern "C" {
    void spiht_encode(float *buffer, size_t height, size_t width, uint8_t **out_buffer, size_t *output_size, int num_stages, float* max_val);
    void spiht_decode(uint8_t *buffer, size_t size, float *out_buffer, size_t height, size_t width, int num_bits, float max_val);
}

#endif // SPIHT_H