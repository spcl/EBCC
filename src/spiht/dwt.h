#include "spiht_re.h"

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

DWTData* empty_dwtcoeff(size_t num_stages, size_t size_x, size_t size_y, size_t extra_x, size_t extra_y) {
    DWTData* dwt_data = (DWTData*) malloc(sizeof(DWTData));
    dwt_data->size_x = size_x;
    dwt_data->size_y = size_y;
    dwt_data->extra_x = extra_x;
    dwt_data->extra_y = extra_y;
    dwt_data->stride = size_x + extra_x;
    dwt_data->num_stages = num_stages;
    dwt_data->data = (elem_t*) calloc((size_x+extra_x) * (size_y+extra_y), sizeof(elem_t));
    dwt_data->temp = (elem_t*) calloc((size_x+extra_x) * (size_y+extra_y), sizeof(elem_t));
    dwt_data->type = DWTCOEFF;
    return dwt_data;
}

DWTData* load_image(elem_t* buffer, size_t height, size_t width, size_t num_stages) {
    DWTData* dwt_data = (DWTData*) malloc(sizeof(DWTData));
    size_t extra_x = 0, extra_y = 0;
    size_t size_x = width, size_y = height;
    size_t stride;
    /* TODO: support different scaling factor */
    elem_t scale = MAXELEM; 
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

void idwt_row(DWTData* dwt_data, size_t row, size_t num_items) {
    const elem_t alpha = -1.586134342;
    const elem_t beta = -0.05298011854;
    const elem_t gamma = 0.8829110762;
    const elem_t delta = 0.44355068522;
    const elem_t xi = 1.149604398;
    size_t stride = dwt_data->stride;
    size_t row_stride = row * stride;
    elem_t* data = dwt_data->data;
    elem_t* temp = dwt_data->temp;

    for (size_t x = 0; x < num_items/2; x++) {
        temp[x + row_stride] /= xi;
        temp[num_items/2 + x + row_stride] *= xi;
    }

    for (size_t x = 1; x < num_items/2; x++)
        temp[x + row_stride] -= delta*(temp[num_items/2 + x + row_stride] + temp[num_items/2 + x - 1 + row_stride]);
    temp[0 + row_stride] -= delta*(temp[num_items/2 + row_stride] + temp[num_items/2 + 1 + row_stride]);

    temp[num_items - 1 + row_stride] -= gamma*(temp[num_items/2 - 1 + row_stride] + temp[num_items/2 - 2 + row_stride]);
    for (size_t x = 0; x < num_items/2 - 1; x++)
        temp[num_items/2 + x + row_stride] -= gamma*(temp[x + row_stride] + temp[x+1 + row_stride]);

    for (size_t x = 1; x < num_items/2; x++)
        data[2*x + row_stride] = temp[x + row_stride] - beta * (temp[num_items/2 + x + row_stride] + temp[num_items/2 + x - 1 + row_stride]);
    data[0 + row_stride] = temp[0 + row_stride] - beta * (temp[num_items/2 + row_stride] + temp[num_items/2 + 1 + row_stride]);

    data[num_items-1 + row_stride] = temp[num_items - 1 + row_stride] - 2*alpha*data[num_items-2 + row_stride];
    for (size_t x = 0; x < num_items/2 - 1; x++)
        data[2*x+1 + row_stride] = temp[num_items/2 + x + row_stride] - alpha*(data[2*x + row_stride] + data[2*x + 2 + row_stride]);
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

void idwt_col(DWTData* dwt_data, size_t col, size_t num_items) {
    const elem_t alpha = -1.586134342;
    const elem_t beta = -0.05298011854;
    const elem_t gamma = 0.8829110762;
    const elem_t delta = 0.44355068522;
    const elem_t xi = 1.149604398;
    elem_t* data = dwt_data->data;
    elem_t* temp = dwt_data->temp;
    size_t stride = dwt_data->stride;

    for (size_t y = 0; y < num_items/2; y++) {
        data[col + (y) * stride] /= xi;
        data[col + (num_items/2+y) * stride] *= xi;
    }

    for (size_t y = 1; y < num_items/2; y++)
        data[col + (y) * stride] -= delta*(data[col + (num_items/2+y) * stride] + data[col + (num_items/2 + y - 1) * stride]);
    data[col + (0) * stride] -= delta*(data[col + (num_items/2) * stride] + data[col + (num_items/2 + 1) * stride]);

    data[col + (num_items - 1) * stride] -= gamma*(data[col + (num_items/2 - 1) * stride] + data[col + (num_items/2-2) * stride]);
    for (size_t y = 0; y < num_items/2 - 1; y++)
        data[col + (num_items/2 + y) * stride] -= gamma*(data[col + (y) * stride] + data[col + (y+1) * stride]);

    for (size_t y = 1; y < num_items/2; y++)
        temp[col + (2*y) * stride] = data[col + (y) * stride] - beta * (data[col + (num_items/2+y) * stride] + data[col + (num_items/2+y-1) * stride]);
    temp[col + (0) * stride] = data[col + (0) * stride] - beta * (data[col + (num_items/2) * stride] + data[col + (num_items/2+1) * stride]);

    temp[col + (num_items-1) * stride] = data[col + (num_items - 1) * stride] - 2*alpha*temp[col + (num_items-2) * stride];
    for (size_t y = 0; y < num_items/2 - 1; y++)
        temp[col + (2*y+1) * stride] = data[col + (num_items/2 + y) * stride] - alpha*(temp[col + (2*y) * stride] + temp[col + (2*y+2) * stride]);
}

void dwt2_(DWTData* dwt_data, size_t size_x, size_t size_y) {
    assert(size_x >= 4 && size_y >= 4);
    for (size_t y = 0; y < size_y; y++)
        dwt_row(dwt_data, y, size_x);
    for (size_t x = 0; x < size_x; x++)
        dwt_col(dwt_data, x, size_y);
}

void idwt2_(DWTData* dwt_data, size_t size_x, size_t size_y) {
    assert(size_x >= 4 && size_y >= 4);
    for (size_t x = 0; x < size_x; x++)
        idwt_col(dwt_data, x, size_y);
    for (size_t y = 0; y < size_y; y++)
        idwt_row(dwt_data, y, size_x);
}

void dwt2full(DWTData* dwt_data) {
    assert(dwt_data->type == IMAGE);
    size_t size_x = dwt_data->size_x + dwt_data->extra_x;
    size_t size_y = dwt_data->size_y + dwt_data->extra_y;
    for (size_t i = 0; i < dwt_data->num_stages; i++) {
        dwt2_(dwt_data, size_x, size_y);
        size_x /= 2;
        size_y /= 2;
    }
    dwt_data->type = DWTCOEFF;
}

void idwt2full(DWTData* dwt_data) {
    size_t num_stages = dwt_data->num_stages;
    assert(dwt_data->type == DWTCOEFF);
    assert(num_stages > 0);
    size_t size_x = (dwt_data->size_x + dwt_data->extra_x) / (1 << (num_stages - 1));
    size_t size_y = (dwt_data->size_y + dwt_data->extra_y) / (1 << (num_stages - 1));
    for (size_t i = 0; i < num_stages; i++) {
        idwt2_(dwt_data, size_x, size_y);
        size_x *= 2;
        size_y *= 2;
    }
    dwt_data->type = IMAGE;
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

void add_dc(DWTData* dwt_data, elem_t dc) {
    size_t size_x = dwt_data->size_x + dwt_data->extra_x;
    size_t size_y = dwt_data->size_y + dwt_data->extra_y;
    size_t stride = dwt_data->stride;
    assert(dwt_data->type == IMAGE);
    for (size_t y = 0; y < size_y; y++)
        for (size_t x = 0; x < size_x; x++) {
            elem_t val = dwt_data->data[x + y * stride];
            /* TODO: Maybe remove this floor? */
            val = floor(val + dc);
            if (val > MAXELEM) {
                val = MAXELEM;
            } else if (val < 0) {
                val = 0;
            }
            dwt_data->data[x + y * stride] = val;
        }
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