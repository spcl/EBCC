#ifndef HDF5_STUB_H
#define HDF5_STUB_H

// HDF5 version 1.14.6-3.fc42
// Extracted HDF5 symbols required by h5z_j2k.c

#include <stddef.h>
#include <stdint.h>

/* Filter types and constants */
typedef int H5Z_filter_t;

/* Filter flags */
/**
 * reverse direction; read
 */
#define H5Z_FLAG_REVERSE 0x0100

/* Filter class version */
#define H5Z_CLASS_T_VERS 1

/* Filter function type */
typedef size_t (*H5Z_func_t)(unsigned int flags, size_t cd_nelmts, 
                             const unsigned int cd_values[], size_t nbytes, 
                             size_t *buf_size, void **buf);

/* Filter class structure */
typedef struct H5Z_class2_t {
    int                  version;         /* Version number of the H5Z_class_t struct */
    H5Z_filter_t         id;              /* Filter ID number */
    unsigned             encoder_present; /* Does this filter have an encoder? */
    unsigned             decoder_present; /* Does this filter have a decoder? */
    const char          *name;            /* Filter name for debugging */
    void                *can_apply;       /* The "can apply" callback */
    void                *set_local;       /* The "set local" callback */
    H5Z_func_t           filter;          /* The actual filter function */
} H5Z_class2_t;

/* Plugin type enumeration */
typedef enum H5PL_type_t {
    H5PL_TYPE_ERROR  = -1, /* Error */
    H5PL_TYPE_FILTER = 0,  /* Filter */
    H5PL_TYPE_VOL    = 1,  /* VOL connector */
    H5PL_TYPE_VFD    = 2,  /* VFD */
    H5PL_TYPE_NONE   = 3   /* Sentinel: This must be last! */
} H5PL_type_t;

#endif /* HDF5_STUB_H */
