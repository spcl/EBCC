#include "spiht_re.h"

typedef enum {
    BITREAD,
    BITWRITE,
    BITWRITEEND,
} BitIOStreamMode;

struct BitIOStreamState {
    uint8_t byte;
    uint8_t* buffer;
    size_t buffer_size;
    size_t curr;
    size_t bit_pos;
    BitIOStreamMode mode;
};
typedef struct BitIOStreamState BitIOStreamState;

BitIOStreamState* bitio_init(uint8_t* buffer, size_t buffer_size, BitIOStreamMode mode) {
    BitIOStreamState* state = (BitIOStreamState*) malloc(sizeof(BitIOStreamState));
    state->buffer = buffer;
    state->buffer_size = buffer_size;
    state->curr = 0;
    state->mode = mode;
    state->byte = 0;
    state->bit_pos = 0;
    return state;
}

void bitio_free(BitIOStreamState* state) {
    free(state);
}

void bitio_put_bit(BitIOStreamState* state, uint8_t bit) {
#ifdef DEBUG
    assert(state->mode == BITWRITE);
#endif
    state->byte = (state->byte << 1) | (bit & 1);
    state->bit_pos++;
    if (state->bit_pos == 8) {
        assert(state->curr < state->buffer_size);
        state->buffer[state->curr++] = state->byte;
        state->byte = 0;
        state->bit_pos = 0;
    }
}

void bitio_put_bits(BitIOStreamState* state, uint64_t bits, size_t num_bits) {
    assert(state->mode == BITWRITE);
    assert(num_bits <= sizeof(bits) * 8);
    for (size_t i = num_bits - 1; i >= 1 ; i--) {
        bitio_put_bit(state, (bits >> i) & 1);
    }
    bitio_put_bit(state, bits & 1);
}

uint8_t bitio_get_bit(BitIOStreamState* state) {
#ifdef DEBUG
    assert(state->mode == BITREAD);
#endif
    if (state->bit_pos == 0) {
        if (state->curr >= state->buffer_size)
            return 0;
        state->byte = state->buffer[state->curr++];
        state->bit_pos = 8;
    }
    return (state->byte >> --(state->bit_pos)) & 1;
}

uint64_t bitio_get_bits(BitIOStreamState* state, size_t num_bits) {
    assert(state->mode == BITREAD);
    uint64_t bits = 0;
    for (size_t i = 0; i < num_bits; i++)
        bits = (bits << 1) | bitio_get_bit(state);
    return bits;
}

void bitio_flush(BitIOStreamState* state) {
    assert(state->mode == BITWRITE);
    assert(state->bit_pos < 8);
    if (state->bit_pos > 0) {
        assert(state->curr < state->buffer_size);
        state->byte = state->byte << (8 - state->bit_pos);
        state->buffer[state->curr++] = state->byte;
        state->byte = 0;
        state->bit_pos = 0;
    }
    state->mode = BITWRITEEND;
}

size_t bitio_get_size(BitIOStreamState* state) {
    return state->curr;
}