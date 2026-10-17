#ifndef EBCC_QUANTIZE_H
#define EBCC_QUANTIZE_H
#include <math.h>
/* Round to nearest, ties to even, independently of the host rounding mode. */
static inline double ebcc_nearest(double value) {
    double magnitude = fabs(value), index = floor(magnitude);
    double fraction = magnitude - index;
    if (fraction > 0.5 || (fraction == 0.5 && fmod(index, 2.0) != 0)) ++index;
    return copysign(index, value);
}
#endif
