#include <cstdlib>
#include <cstdio>
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

extern "C" void spiht_encode_file(char *filename, uint8_t **out_buffer, size_t *output_size, float bit_rate) {
    Encoder encoder;
    char* fileout = "tmp.ims";
    encoder.encode(filename, fileout, bit_rate);
    FILE *file = fopen(fileout, "rb");
    if (file == NULL) {
        *out_buffer = NULL;
        return;
    }

    fseek(file, 0, SEEK_END);
    size_t file_size = ftell(file);
    fseek(file, 0, SEEK_SET);

    *out_buffer = (uint8_t *)malloc(file_size);
    if (*out_buffer == NULL) {
        fclose(file);
        return;
    }

    fread(*out_buffer, 1, file_size, file);
    *output_size = file_size;
    fclose(file);
    remove(fileout);
}

extern "C" void spiht_decode_file(uint8_t *in_buffer, size_t buffer_size, char *fileout) {
    char* filein = "tmp.ims";
    FILE *file = fopen(filein, "wb");
    if (file == NULL) {
        return;
    }
    fwrite(in_buffer, 1, buffer_size, file);
    fclose(file);
    Decoder decoder;
    decoder.decode(filein, fileout);
    remove(filein);
}