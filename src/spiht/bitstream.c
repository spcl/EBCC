#include "bitstream.h"

void initialize_bitstream(BitStreamState *state, uint8_t *buffer, size_t size) {
    state->buffer = buffer;
    state->bit_count = 0;
    state->pos = 0;
    state->size = size;
    state->curr = 0;
}

uint8_t get_bit(BitStreamState *state) {
    if (state->bit_count == 0) {
        if (state->pos >= state->size) {
            return 0; // End of buffer
        }
        state->curr = state->buffer[state->pos++];
        state->bit_count = 8;
    }
    uint8_t bit = (state->curr >> (state->bit_count - 1)) & 1;
    state->bit_count--;
    return bit;
}

uint32_t get_bits(BitStreamState *state, uint8_t size) {
    uint32_t bits = 0;
    for (uint8_t i = 0; i < size; i++) {
        bits = (bits << 1) | get_bit(state);
    }
    return bits;
}

void put_bit(BitStreamState *state, uint8_t bit) {
    state->curr = (state->curr << 1) | (bit & 1);
    state->bit_count++;
    if (state->bit_count == 8) {
        if (state->pos < state->size) {
            state->buffer[state->pos++] = state->curr;
        }
        state->bit_count = 0;
        state->curr = 0;
    }
}

void put_bits(BitStreamState *state, uint32_t bits, uint8_t size) {
    for (int8_t i = size - 1; i >= 0; i--) {
        put_bit(state, (bits >> i) & 1);
    }
}

void flush_bits(BitStreamState *state) {
    while (state->bit_count > 0)
        put_bit(state, 0);
}

void flush(BitStreamState *state) {
    flush_bits(state);
}