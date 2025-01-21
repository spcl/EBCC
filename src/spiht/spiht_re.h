
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <assert.h>

#define MAXSTEPS 32
#define MAXINT 2147483646

typedef int64_t int_t;
typedef char bool_t;
typedef float elem_t;

void encode_image(elem_t *buffer, size_t height, size_t width, uint8_t** buffer_out, size_t* output_size, size_t trunc_bits, size_t num_stages);