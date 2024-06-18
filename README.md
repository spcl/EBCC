# How to install
```
cd src
# download and install dependency libraries (openjpeg and wavelib)
# HDF5 and NetCDF4 library (with headers) have to be installed manually
make setup
# compile the HDF5 filter for compression
make
# the compiled filter is stored in `src/libh5z_j2k.so`
```

# User guide
Our compression algorithm is implemented as an HDF5 filter. The filter expects a chunk size equal to a single 2D “frame”, and can scale to any number of dimensions greater or equan than 2.

Its base functionality compresses the data using JPEG2000. The user provides a "base" compression ratio for this functionality.
The user can also enable a compression of the residual in order to improve accuracy. The residual is the difference between the original frame and the uncompressed frame. The residual is compressed independently of the base frame, and is summed to the uncompressed base frame when decompressing. We support three modes of operation:
1. `NONE`: there is no residual
2. `SPARSIFICATION_FACTOR`: the residual is wavelet encoded and sparsified with the desired sparsification factor. The sparsification factor roughly corresponds to the compression ratio of the residual
3. `MAX_ERROR`: the residual is wavelet encoded and sparsified. The sparsification factor is found through an iterative process that tries out several factors and selects the largest one that keeps the max error below the selected threshold

The input parameters depend on the chosen mode of operation. As HDF5 filters support integer parameters only, a translation from float and double types to integer representation is required. We provide a python wrapper `JP2SPWV_Filter` in `filter_wrapper.py` to simply the process. Here is an example of using the filter (also in `test.py`):

```python
import os
import h5py
import numpy as np
# specify the path of compiled HDF5 plugin libh5z_j2k.so
os.environ["HDF5_PLUGIN_PATH"] = os.path.join(os.getcwd(), 'src')

# generate data
chunk_shape = (721, 1440)
data = np.zeros((10, *chunk_shape))

# initialize the wrapper
jp2spwv_filter = JP2SPWV_Filter(
    base_cr=100, # base compression ratio
    height=chunk_shape[0],  # height of each 2D data chunk
    width=chunk_shape[1],  # width of each 2D data chunk
    data_dim=len(data.shape), # data dimension, required to specify the HDF5 chunk shape
    residual_opt=("max_error_target", 1.0)) # specify the max error target to be 1.0
    # other possible residual_opt can be
    # `("quantile_error_target", xxx)` : max_error does not exceed the specified quantile value calculated from the compression error with only base compression method
    # `("fixed_sparsification", xxx)`: specify a fixed sparsification ratio for the sparse wavelet compression

# create hdf5 storage
f = h5py.File('test.hdf5', 'a') 
f.create_dataset('compressed', shape=data.shape, dtype=data.dtype, **jp2spwv_filter)

# compress data
f['compressed'][:] = data
# compressed data now stored in `test.hdf5`

```


