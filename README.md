# How to install
```
cd src
make setup
make
```

# User guide
Our compression algorithm is implemented as an HDF5 filter. The filter expects a chunk size equal to a single 2D “frame”, and can scale to any number of dimensions greater or equan than 2.

Its base functionality compresses the data using JPEG2000. The user provides a "base" compression ratio for this functionality.
The user can also enable a compression of the residual in order to improve accuracy. The residual is the difference between the original frame and the uncompressed frame. The residual is compressed independently of the base frame, and is summed to the uncompressed base frame when decompressing. We support three modes of operation:
1. NONE: there is no residual
2. SPARSIFICATION_FACTOR: the residual is wavelet encoded and sparsified with the desired sparsification factor. The sparsification factor roughly corresponds to the compression ratio of the residual
3. MAX_ERROR: the residual is wavelet encoded and sparsified. The sparsification factor is found through an iterative process that tries out several factors and selects the largest one that keeps the max error below the selected threshold

The input parameters depend on the chosen mode of operation. As HDF5 filters support integer parameters only, a translation from float and double types to integer representation is required. The script `params_helper.py` can be used to generate these parameters.

