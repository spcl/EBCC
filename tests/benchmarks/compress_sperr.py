import os
import h5py
import hdf5plugin
import numpy as np
import xarray

f = h5py.File('geopotential_pl_small_sperr.nc', 'w')
orig = xarray.open_dataset('geopotential_pl_small.nc')
data = orig['z'].data

data_comp = f.create_dataset('z', shape=data.shape, dtype=np.float32, chunks = (1, 1, data.shape[-2], data.shape[-1]), **hdf5plugin.Sperr(absolute=10.0))
f['z'][...] = data[...]

f.close()

original_size = os.path.getsize('geopotential_pl_small.nc')
compressed_size = os.path.getsize('geopotential_pl_small_sperr.nc')

print(f'SPERR: achieved compression ratio of {original_size/compressed_size}')

