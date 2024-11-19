#ifndef BITSTREAM_H
#define BITSTREAM_H

#include <stdio.h>
#include <stdint.h>

typedef struct {
    uint8_t *buffer;
    uint8_t bit_count; // current position in bits (0-7)
    size_t pos; // current position in bytes
    size_t size; // const, size of buffer
    uint8_t curr; // current byte
} BitStreamState;

void initialize_bitstream(BitStreamState *state, uint8_t *buffer, size_t size);
uint8_t get_bit(BitStreamState *state);
uint32_t get_bits(BitStreamState *state, uint8_t size);
void put_bit(BitStreamState *state, uint8_t bit);
void put_bits(BitStreamState *state, uint32_t bits, uint8_t size);
void flush_bits(BitStreamState *state);
void flush(BitStreamState *state);

#endif // BITSTREAM_H