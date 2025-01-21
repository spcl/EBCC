#include "spiht_re.h"

struct MaskedList {
    size_t size;
    size_t capacity;
    size_t curr;
    bool_t* mask;
    int_t* values;
};
typedef struct MaskedList MaskedList;

MaskedList* ml_init(size_t capacity) {
    MaskedList* ml = (MaskedList*) malloc(sizeof(MaskedList));
    ml->size = 0;
    ml->curr = 0;
    ml->capacity = capacity;
    ml->mask = (bool_t*) calloc(capacity, sizeof(bool_t));
    ml->values = (int_t*) calloc(capacity, sizeof(int_t));
    return ml;
}

void ml_free(MaskedList* ml) {
    free(ml->mask);
    free(ml->values);
    free(ml);
}

inline void ml_push(MaskedList* ml, int_t value) {
#ifdef DEBUG
    assert(ml->curr < ml->capacity);
#endif
    ml->values[ml->curr++] = value;
    ml->size++;
}

inline int_t ml_get(MaskedList* ml, size_t index) {
#ifdef DEBUG
    assert(index < ml->curr);
    assert(ml->mask[index] == 0);
#endif
    return ml->values[index];
}

inline void ml_remove(MaskedList* ml, size_t index) {
#ifdef DEBUG
    assert(index < ml->curr);
#endif
    ml->mask[index] = 1;
    ml->size--;
}

void ml_consolidate(MaskedList* ml) {
    size_t j = 0;
    for (size_t i = 0; i < ml->curr; i++) {
        if (!(ml->mask[i]) && i != j++) {
            ml->values[j - 1] = ml->values[i];
        } else {
            ml->mask[i] = 0;
        }
    }
    /* memset(ml->mask, 0, ml->curr * sizeof(bool_t)); */
    ml->curr = j;
    assert(j == ml->size);
}

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
    assert(state->mode == BITWRITE);
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

typedef enum {
    IMAGE,
    DWTCOEFF,
} MatrixType;

struct DWTData {
    size_t size_x;
    size_t size_y;
    size_t extra_x; /* x padding size for DWT */
    size_t extra_y; /* y padding size for DWT */
    size_t stride;
    size_t num_stages; /* Number of DWT stages*/
    elem_t* data;
    elem_t* temp;
    MatrixType type;
};
typedef struct DWTData DWTData;

struct EncoderState {
    DWTData* dwt_data;
    MaskedList* lip; /* List of insignificant pixels */
    MaskedList* lsp; /* List of significant pixels */
    MaskedList* lis; /* List of insignificant sets */
    int_t step; /* Quantization step */
};
typedef struct EncoderState EncoderState;

DWTData* load_image(elem_t* buffer, size_t height, size_t width, size_t num_stages) {
    DWTData* dwt_data = (DWTData*) malloc(sizeof(DWTData));
    size_t extra_x = 0, extra_y = 0;
    size_t size_x = width, size_y = height;
    size_t stride;
    /* TODO: support different scaling factor */
    elem_t scale = 255.0; /* Scale factor for 8-bit image data */
    while ((size_x + extra_x) % (1 << (num_stages+1)) != 0)
        extra_x++;
    while ((size_y + extra_y) % (1 << (num_stages+1)) != 0)
        extra_y++;
    stride = size_x + extra_x;

    dwt_data->size_x = size_x;
    dwt_data->size_y = size_y;
    dwt_data->num_stages = num_stages;
    dwt_data->stride = stride;
    dwt_data->extra_x = extra_x;
    dwt_data->extra_y = extra_y;
    dwt_data->data = (elem_t*) calloc((size_x+extra_x) * (size_y+extra_y), sizeof(elem_t));
    dwt_data->temp = (elem_t*) calloc((size_x+extra_x) * (size_y+extra_y), sizeof(elem_t));
    dwt_data->type = IMAGE;
    for (size_t y = 0; y < size_y; y++)
        for (size_t x = 0; x < size_x; x++)
            dwt_data->data[x + y * stride] = buffer[y * width + x] * scale;

    /* Symmetrize */
    for (size_t y = 0; y < size_y; y++)
        for (size_t x = 0; x < extra_x; x++)
            dwt_data->data[size_x + x + y * stride] = dwt_data->data[size_x - x - 1 + y * stride];
    for (size_t x = 0; x < size_x; x++)
        for (size_t y = 0; y < extra_y; y++)
            dwt_data->data[x + (size_y + y) * stride] = dwt_data->data[x + (size_y - y - 1) * stride];
    for (size_t y = size_y; y < size_y + extra_y; y++)
        for (size_t x = size_x; x < size_x + extra_x; x++)
            dwt_data->data[x + y * stride] = 0;
    return dwt_data;
}

void free_image(DWTData* dwt_data) {
    free(dwt_data->data);
    free(dwt_data->temp);
    free(dwt_data);
}


void dwt_row(DWTData* dwt_data, size_t row, size_t num_items) {
    const elem_t alpha = -1.586134342;
    const elem_t beta = -0.05298011854;
    const elem_t gamma = 0.8829110762;
    const elem_t delta = 0.44355068522;
    const elem_t xi = 1.149604398;
    elem_t* data = dwt_data->data;
    elem_t* temp = dwt_data->temp;
    size_t stride = dwt_data->stride;
    size_t row_stride = row * stride;
    for (size_t x = 0; x < num_items / 2 - 1; x++)
        temp[num_items/2 + x + row_stride] = data[2*x+1 + row_stride] + alpha*(data[2*x + row_stride] + data[2*x + 2 + row_stride]);
    temp[num_items - 1 + row_stride] = data[num_items - 1 + row_stride] + 2*alpha*data[num_items - 2 + row_stride];

    temp[0 + row_stride] = data[0 + row_stride] + beta * (temp[num_items/2 + row_stride] + temp[num_items/2 + 1 + row_stride]);
    for (size_t x = 1; x < num_items/2; x++)
        temp[x + row_stride] = data[2*x + row_stride] + beta * (temp[num_items/2 + x + row_stride] + temp[num_items/2 + x - 1 + row_stride]);

    for (size_t x = 0; x < num_items/2 - 1; x++)
        temp[num_items/2 + x + row_stride] += gamma*(temp[x + row_stride] + temp[x+1 + row_stride]);
    temp[num_items - 1 + row_stride] += gamma*(temp[num_items/2 - 1 + row_stride] + temp[num_items/2 - 2 + row_stride]);

    temp[0 + row_stride] += delta*(temp[num_items/2 + row_stride] + temp[num_items/2 + 1 + row_stride]);
    for (size_t x = 1; x < num_items/2; x++)
        temp[x + row_stride] += delta*(temp[num_items/2 + x + row_stride] + temp[num_items/2 + x - 1 + row_stride]);

    for (size_t x = 0; x < num_items/2; x++) {
        temp[x + row_stride] *= xi;
        temp[num_items/2 + x + row_stride] /= xi;
    }
}

void dwt_col(DWTData* dwt_data, size_t col, size_t num_items) {
    const elem_t alpha = -1.586134342;
    const elem_t beta = -0.05298011854;
    const elem_t gamma = 0.8829110762;
    const elem_t delta = 0.44355068522;
    const elem_t xi = 1.149604398;
    elem_t* data = dwt_data->data;
    elem_t* temp = dwt_data->temp;
    size_t stride = dwt_data->stride;
    
    assert(num_items >= 4);
    for (size_t y = 0; y < num_items/2 - 1; y++)
        data[col + (num_items/2 + y) * stride] = temp[col + (2*y+1) * stride] + alpha*(temp[col + (2*y) * stride] + temp[col + (2*y+2) * stride]);
    data[col + (num_items - 1) * stride] = temp[col + (num_items-1) * stride] + 2*alpha*temp[col + (num_items-2) * stride];

    data[col + (0) * stride] = temp[col + (0) * stride] + beta * (data[col + (num_items/2) * stride] + data[col + (num_items/2+1) * stride]);
    for (size_t y = 1; y < num_items/2; y++)
        data[col + (y) * stride] = temp[col + (2*y) * stride] + beta * (data[col + (num_items/2+y) * stride] + data[col + (num_items/2+y-1) * stride]);

    for (size_t y = 0; y < num_items/2 - 1; y++)
        data[col + (num_items/2 + y) * stride] += gamma*(data[col + (y) * stride] + data[col + (y+1) * stride]);
    data[col + (num_items - 1) * stride] += gamma*(data[col + (num_items/2 - 1) * stride] + data[col + (num_items/2-2) * stride]);

    data[col + (0) * stride] += delta*(data[col + (num_items/2) * stride] + data[col + (num_items/2 + 1) * stride]);
    for (size_t y = 1; y < num_items/2; y++)
        data[col + (y) * stride] += delta*(data[col + (num_items/2+y) * stride] + data[col + (num_items/2 + y - 1) * stride]);

    for (size_t y = 0; y < num_items/2; y++) {
        data[col + (y) * stride] *= xi;
        data[col + (num_items/2+y) * stride] /= xi;
    }

}

void dwt2(DWTData* dwt_data, size_t size_x, size_t size_y) {
    for (size_t y = 0; y < size_y; y++)
        dwt_row(dwt_data, y, size_x);
    for (size_t x = 0; x < size_x; x++)
        dwt_col(dwt_data, x, size_y);
}


void dwt2full(DWTData* dwt_data) {
    assert(dwt_data->type == IMAGE);
    size_t size_x = dwt_data->size_x + dwt_data->extra_x;
    size_t size_y = dwt_data->size_y + dwt_data->extra_y;
    for (size_t i = 0; i < dwt_data->num_stages; i++) {
        dwt2(dwt_data, size_x, size_y);
        size_x /= 2;
        size_y /= 2;
    }
    dwt_data->type = DWTCOEFF;
}

elem_t sub_dc(DWTData* dwt_data) {
    size_t size_x = dwt_data->size_x + dwt_data->extra_x;
    size_t size_y = dwt_data->size_y + dwt_data->extra_y;
    size_t stride = dwt_data->stride;
    assert(dwt_data->type == IMAGE);
    double dc = 0;
    for (size_t y = 0; y < size_y; y++)
        for (size_t x = 0; x < size_x; x++)
            dc += dwt_data->data[x + y * stride];
    dc /= (size_x * size_y);
    dc = floor(dc);
    for (size_t y = 0; y < size_y; y++)
        for (size_t x = 0; x < size_x; x++)
            dwt_data->data[x + y * stride] -= dc;
    return (elem_t) dc;
}

void normalize(DWTData* dwt_data) {
    size_t size_x = dwt_data->size_x + dwt_data->extra_x;
    size_t size_y = dwt_data->size_y + dwt_data->extra_y;
    size_t stride = dwt_data->stride;
    assert(dwt_data->type == DWTCOEFF);
    for (size_t y = 0; y < size_y; y++)
        for (size_t x = 0; x < size_x; x++) {
            if (dwt_data->data[x + y * stride] >= 0) {
                dwt_data->data[x + y * stride] = floor(dwt_data->data[x + y * stride]);
            } else {
                dwt_data->data[x + y * stride] = -floor(fabs(dwt_data->data[x + y * stride]));
            }
        }
}

void free_encoder_state(EncoderState* state) {
    ml_free(state->lip);
    ml_free(state->lsp);
    ml_free(state->lis);
    free(state);
}

EncoderState* spiht_encode_init(BitIOStreamState* bio, DWTData* dwt_data, size_t trunc_bits) {
    EncoderState* state = (EncoderState*) malloc(sizeof(EncoderState));
    elem_t max = 2.0; /* step = log2(max) is at least 0 */
    size_t size_x = dwt_data->size_x;
    size_t size_y = dwt_data->size_y;
    size_t extra_x = dwt_data->extra_x;
    size_t extra_y = dwt_data->extra_y;
    size_t stride = dwt_data->stride;
    size_t num_stages = dwt_data->num_stages;
    size_t first_stage_size_x = (size_x + extra_x) / (1 << num_stages);
    size_t first_stage_size_y = (size_y + extra_y) / (1 << num_stages);
    size_t image_size = (size_x + extra_x) * (size_y + extra_y);
    size_t capacity = (first_stage_size_x * first_stage_size_y) + trunc_bits;
    if (capacity > image_size * 8) {
        /* trunc_bits can go unlimited (2^28), cap it with image data size */
        capacity = image_size * 8;
    }
    assert (capacity <= MAXINT); /* such that x, y, set flag can be encoded into 64 bit int*/
    assert (dwt_data->type == DWTCOEFF);
    state->lip = ml_init(capacity);
    state->lsp = ml_init(capacity);
    state->lis = ml_init(capacity);
    for (size_t i = 0; i < image_size; i++) {
        elem_t absval = fabs(dwt_data->data[i]);
        if (absval > max) {
            max = absval;
        }
    }
    state->step = (int_t) floor(log(max) / log(2.0));
    assert(state->step <= MAXSTEPS);

    bitio_put_bits(bio, (uint64_t) state->step, 8);

    for (size_t y = 0; y < first_stage_size_y; y++) {
        for (size_t x = 0; x < first_stage_size_x; x++) {
            int_t pix_item = x + y * stride;
            /* set_item = (pix_item, A)  
               (pix_item, A) -> pix_item + 1
               (pix_item, B) -> - (pix_item + 1)  */
            int_t set_item = pix_item + 1; 
            ml_push(state->lip, pix_item);
            if ((x % 2 != 0) || (y % 2 != 0)) {
                ml_push(state->lis, set_item);
            }
        }
    }
    return state;
}

inline bool_t is_significant_pixel(int_t step, elem_t val) {
    /* TODO: preprocess val into integer */
    /* Version 1: return (((int_t) fabs(val)) >= (1 << step)); */
    /* Version 2: return (( fabs(val)) >= (double) (1 << step)); */
    /* Version 3: (used in imshrinker)*/
    return abs((int_t) val) >= (1 << step);
}

inline void get_successor(int_t x, int_t y, DWTData* dwt_data, int_t* sx, int_t* sy) {
    size_t size_x_pad = dwt_data->size_x + dwt_data->extra_x;
    size_t size_y_pad = dwt_data->size_y + dwt_data->extra_y;
    int_t lx = size_x_pad / (1 << dwt_data->num_stages);
    int_t ly = size_y_pad / (1 << dwt_data->num_stages);
    int_t sx_, sy_;
    if (x < lx && y < ly) {
        if (x % 2 == 1) {
            sx_ = x + lx - 1;
        } else {
            sx_ = x;
        }
        if (y % 2 == 1) {
            sy_ = y + ly - 1;
        } else {
            sy_ = y;
        }
        if (sx_ == x && sy_ == y) {
            sx_ = -1;
            sy_ = -1;
        }
    } else {
        sx_ = 2 * x;
        sy_ = 2 * y;
        if (sx_ >= size_x_pad || sy_ >= size_y_pad) {
            sx_ = -1;
            sy_ = -1;
        }
    }
    *sx = sx_;
    *sy = sy_;
}

inline bool_t is_significant_set_A(int_t step, DWTData* dwt_data, int_t pix_item, size_t count) {
    size_t stride = dwt_data->stride;
    int_t sx, sy, x, y;
    if (count > 1 && is_significant_pixel(step, dwt_data->data[pix_item])) {
        return 1;
    }
    x = pix_item % stride;
    y = pix_item / stride;
    get_successor(x, y, dwt_data, &sx, &sy);
    if (sx == -1 || sy == -1) {
        return 0;
    }
    if (is_significant_set_A(step, dwt_data, sx + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_A(step, dwt_data, sx + 1 + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_A(step, dwt_data, sx + (sy + 1) * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_A(step, dwt_data, sx + 1 + (sy + 1) * stride, count + 1)) {
        return 1;
    }
    return 0;
}

inline bool_t is_significant_set_B(int_t step, DWTData* dwt_data, int_t pix_item, size_t count) {
    size_t stride = dwt_data->stride;
    int_t sx, sy, x, y;
    if (count > 2 && is_significant_pixel(step, dwt_data->data[pix_item])) {
        return 1;
    }
    x = pix_item % stride;
    y = pix_item / stride;
    get_successor(x, y, dwt_data, &sx, &sy);
    if (sx == -1 || sy == -1) {
        return 0;
    }
    if (is_significant_set_B(step, dwt_data, sx + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_B(step, dwt_data, sx + 1 + sy * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_B(step, dwt_data, sx + (sy + 1) * stride, count + 1)) {
        return 1;
    } else if (is_significant_set_B(step, dwt_data, sx + 1 + (sy + 1) * stride, count + 1)) {
        return 1;
    }
    return 0;
}

void spiht_encode_process(size_t bits, BitIOStreamState* bio, DWTData* dwt_data, EncoderState* state) {
    size_t bit_cnt = 0;
    MaskedList* lip = state->lip;
    MaskedList* lsp = state->lsp;
    MaskedList* lis = state->lis;
    size_t stride = dwt_data->stride;
    assert(dwt_data->type == DWTCOEFF);
    /* TODO: check num bits only in the end of each loop and in put_bit */
    for (int_t step = state->step; step >= 0; step--) {
        /* Sorting pass */
        /* First process LIP */
        for (size_t i = 0; i < lip->curr; i++) {
            int_t pix_item = ml_get(lip, i);
            elem_t val = dwt_data->data[pix_item];
            bool_t sig = is_significant_pixel(step, val);
            bitio_put_bit(bio, (uint8_t) sig);
            if (++bit_cnt > bits) return;
            if (sig) {
                ml_push(lsp, pix_item);
                /* Encode the sign bit */
                bitio_put_bit(bio, (uint8_t) ((val > 0) ? 0 : 1));
                if (++bit_cnt > bits) return;
                ml_remove(lip, i);
            }
        }
        ml_consolidate(lip);

        /* now process LIS */
        for (size_t i = 0; i < lis->curr; i++) {
            int_t set_item = ml_get(lis, i);
#ifdef DEBUG
            assert(set_item != 0);
#endif
            bool_t is_set_A = set_item > 0;
            int_t sx, sy;
            if (is_set_A) {
                int_t pix_item = set_item - 1;
                int_t x = pix_item % stride;
                int_t y = pix_item / stride;
                bool_t sig = is_significant_set_A(step, dwt_data, pix_item, 1);
                bitio_put_bit(bio, (uint8_t) sig);
                if (++bit_cnt > bits) return;
                if (sig) {
                    get_successor(x, y, dwt_data, &sx, &sy);
                    /* process the four offsprings */
                    
                    for (int_t dy = 0; dy < 2; dy++) {
                        for (int_t dx = 0; dx < 2; dx++) {
                            int_t pix_item = sx + dx + (sy + dy) * stride;
                            elem_t val = dwt_data->data[pix_item];
                            sig = is_significant_pixel(step, val);
                            bitio_put_bit(bio, (uint8_t) sig);
                            if (++bit_cnt > bits) return;
                            if (sig) {
                                ml_push(lsp, pix_item);
                                bitio_put_bit(bio, (uint8_t) ((val > 0) ? 0 : 1));
                                if (++bit_cnt > bits) return;
                            } else {
                                ml_push(lip, pix_item);
                            }
                        }
                    }

                    /* test if L(i, j) != 0 */
                    get_successor(sx, sy, dwt_data, &sx, &sy);
                    if (sx != -1) {
                        /* push (sx, sy, B) in lis */
                        int_t set_item_B = - (x + y * stride + 1);
                        ml_push(lis, set_item_B);
                    }

                    ml_remove(lis, i);
                }
            } else {
                /* Set B */
                int_t pix_item = - set_item - 1;
                int_t x = pix_item % stride;
                int_t y = pix_item / stride;
                bool_t sig = is_significant_set_B(step, dwt_data, pix_item, 1);
                bitio_put_bit(bio, (uint8_t) sig);
                if (++bit_cnt > bits) return;
                if (sig) {
                    int_t set_item_A;
                    get_successor(x, y, dwt_data, &sx, &sy);
                    set_item_A = sx + sy * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + 1 + sy * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + (sy + 1) * stride + 1;
                    ml_push(lis, set_item_A);
                    set_item_A = sx + 1 + (sy + 1) * stride + 1;
                    ml_push(lis, set_item_A);
                    ml_remove(lis, i);
                }
            }
        }
        ml_consolidate(lis);

        /* Refinement pass */
        for (size_t i = 0; i < lsp->curr; i++) {
            elem_t val = dwt_data->data[ml_get(lsp, i)];
            if (is_significant_pixel(step+1, val)) {
                bitio_put_bit(bio, (uint8_t) (( ((int_t) abs((int_t) val)) >> step) & 1));
                if (++bit_cnt > bits) return;
            }
        }
        
    }
}

void encode_image(elem_t *buffer, size_t height, size_t width, uint8_t** buffer_out, size_t* output_size, size_t trunc_bits, size_t num_stages) {
    size_t buffer_size = (trunc_bits == 0) ? height * width * sizeof(elem_t) : trunc_bits / sizeof(uint8_t) + 1 ;
    if (buffer_size < 16) buffer_size = 16;
    DWTData* dwt_data = load_image(buffer, height, width, num_stages);
    *buffer_out = (uint8_t*) calloc(buffer_size, sizeof(uint8_t));
    BitIOStreamState* bio = bitio_init(*buffer_out, buffer_size, BITWRITE);
    assert(num_stages <= 32 && num_stages >= 1 );
    size_t size_x = dwt_data->size_x;
    size_t size_y = dwt_data->size_y;
    size_t extra_x = dwt_data->extra_x;
    size_t extra_y = dwt_data->extra_y;
    assert(size_x <= 2047 && size_x >= 1);
    assert(size_y <= 2047 && size_y >= 1);
    assert(extra_x <= 511 && extra_x >= 0);
    assert(extra_y <= 511 && extra_y >= 0);
    /* write IMS header */
    bitio_put_bits(bio, (uint64_t) 'I', 8);
    bitio_put_bits(bio, (uint64_t) 'M', 8);
    bitio_put_bits(bio, (uint64_t) 'S', 8);
    bitio_put_bits(bio, (uint64_t) num_stages, 6);
    bitio_put_bits(bio, (uint64_t) size_x, 12);
    bitio_put_bits(bio, (uint64_t) size_y, 12);
    bitio_put_bits(bio, (uint64_t) extra_x, 10);
    bitio_put_bits(bio, (uint64_t) extra_y, 10);
    bitio_put_bit(bio, (uint8_t) 0); /* is_color = False */
    size_t offset = 128; /* 128 (114) for metadata bits */
    size_t bits0 = (trunc_bits == 0) ? (1 << 28) : trunc_bits + offset;
    assert(bits0 > 0);
    bitio_put_bits(bio, (uint64_t) bits0, 29);
    elem_t dc0 = sub_dc(dwt_data);
    assert(dc0 >= 0 && dc0 <= 255);
    uint8_t dc0_uint = (uint8_t) dc0;
    bitio_put_bits(bio, (uint64_t) dc0_uint, 8);
    /* TODO: dump DWT coeffs and compare */
    dwt2full(dwt_data);
    normalize(dwt_data);
    EncoderState* state = spiht_encode_init(bio, dwt_data, bits0 - offset);
    spiht_encode_process(bits0 - offset, bio, dwt_data, state);
    bitio_flush(bio);
    *output_size = bitio_get_size(bio);
    bitio_free(bio);
    free_image(dwt_data);
    free_encoder_state(state);
}