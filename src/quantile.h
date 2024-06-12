#include <stdlib.h>
#include <math.h>
#include <float.h>

void heapify(double* heap, size_t idx, size_t size) {
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;

    if (left < size && heap[left] < heap[smallest])
        smallest = left;
    if (right < size && heap[right] < heap[smallest])
        smallest = right;

    if (smallest != idx) {
        double tmp = heap[idx];
        heap[idx] = heap[smallest];
        heap[smallest] = tmp;
        heapify(heap, smallest, size);
    }
}

void zero_out_function(double *data, size_t length, double quantile) {
    for (size_t i = 0; i < length; ++i) {
        if (fabs(data[i]) <= quantile) {
            data[i] = 0;
        }
    }
}

// zeroes out all of the values below the q-th quantile
// greater q is greater deletion
double zero_out_quantile(double *data, size_t length, double q_ratio) {
    size_t quantile_index = ((double) length * q_ratio);
    size_t heap_size = length - quantile_index;
    double *heap = (double *) calloc(heap_size, sizeof(double));

    for (size_t i = 0; i < length; ++i) {
        double a = fabs(data[i]);
        if (a > heap[0]) {
            heap[0] = a;
            heapify(heap, 0, heap_size);
        }
    }

    double quantile = heap[0];
    free(heap);

    zero_out_function(data, length, quantile);

    return quantile;
}
