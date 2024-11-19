#include <cstdlib>
#include "Encoder.h"
#include "Decoder.h"
#include "spiht.h"


extern "C" void spiht_encode(float *buffer, size_t height, size_t width, uint8_t **out_buffer, size_t *output_size, int num_stages) {
    Encoder encoder;
    encoder.encode_image(buffer, height, width, out_buffer, output_size, num_stages);
}


extern "C" void spiht_decode(uint8_t *buffer, size_t size, float *out_buffer, size_t height, size_t width, int num_bits) {
    Decoder decoder;
    decoder.decode_image(buffer, size, out_buffer, height, width, num_bits);
}

extern "C" void spiht_destroy_buffer(uint8_t *buffer) {
    free(buffer);
}