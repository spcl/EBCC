#include <cstdlib>
#include "Encoder.h"
#include "Decoder.h"
#include "spiht.h"


extern "C" void spiht_encode(float *buffer, size_t height, size_t width, uint8_t **out_buffer, size_t *output_size, int num_stages, float* max_val) {
    Encoder encoder;
    encoder.encode_image(buffer, height, width, out_buffer, output_size, num_stages, max_val);
}


extern "C" void spiht_decode(uint8_t *buffer, size_t size, float *out_buffer, size_t height, size_t width, int num_bits, float max_val) {
    Decoder decoder;
    decoder.decode_image(buffer, size, out_buffer, height, width, num_bits, max_val);
}

extern "C" void spiht_destroy_buffer(uint8_t *buffer) {
    free(buffer);
}