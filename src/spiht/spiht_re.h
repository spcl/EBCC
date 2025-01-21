#ifndef SPIHT_RE_H
#define SPIHT_RE_H
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>

/* Scale factor for 8-bit image data */
#define MAXELEM 255
#define MAXSTEPS 32
#define MAXINT 2147483646

typedef int64_t int_t;
typedef char bool_t;
typedef float elem_t;

void spiht_encode(elem_t *buffer, size_t height, size_t width, uint8_t** buffer_out, size_t* output_size, size_t trunc_bits, size_t num_stages);
void spiht_decode(uint8_t* buffer_in, size_t input_size, elem_t* buffer_out, size_t height, size_t width, size_t num_bits);

#endif /* SPIHT_RE_H */ 