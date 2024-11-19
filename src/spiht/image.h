#ifndef IMAGE_H
#define IMAGE_H

#include <stddef.h>

void loadPGM(const char *filename, size_t *width, size_t *height, unsigned int *maxVal, void **data);

#endif // IMAGE_H