import os
import pyproj # Avoids HDF5 plugin error
from ebcc import EBCC_FILTER_DIR

os.environ["HDF5_PLUGIN_PATH"] = EBCC_FILTER_DIR
import h5py
import numpy as np
from ebcc.filter_wrapper import EBCC_Filter
import xarray

f = h5py.File(f'geopotential_pl_small_ebcc.nc', 'w')
orig = xarray.open_dataset(f'geopotential_pl_small.nc')
data = orig['z'].data

ebcc_filter = EBCC_Filter(
    base_cr=30, # base compression ratio
    height=data.shape[-2],  # height of each 2D data chunk
    width=data.shape[-1],  # width of each 2D data chunk
    data_dim=len(data.shape), # data dimension, required to specify the HDF5 chunk shape
    residual_opt=("max_error_target", 10.0)) # specify the maximum absolute error target
    # other possible residual_opt can be
    # `("max_error_target", xxx)` : the max_error does not exceed the specified value

ebcc_filter = dict(ebcc_filter)
dtype = ebcc_filter.pop('dtype')
data_comp = f.create_dataset('z', shape=data.shape, dtype=dtype, **ebcc_filter)
f['z'][...] = data[...]
#for name,val in data.attrs.items():
#    data_comp.attrs[name] = val
#f.create_dataset('time', data=orig['time'][...], dtype=orig['time'].dtype)
#f.create_dataset('longitude', data=orig['longitude'][...])
#f.create_dataset('latitude', data=orig['latitude'][...])
#f.create_dataset('level', data=orig['level'][...])

f.close()

original_size = os.path.getsize(f'geopotential_pl_small.nc')
compressed_size = os.path.getsize(f'geopotential_pl_small_ebcc.nc')

print(f'EBCC: achieved compression ratio of {original_size/compressed_size}')
