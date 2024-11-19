import os
base_folder = "/home/huanglangwen/Documents/compression-filter"
#os.environ["HDF5_PLUGIN_PATH"] = str(os.path.join(base_folder, 'src/build/lib')) + ":" + str(plugins_folder)
os.environ["HDF5_PLUGIN_PATH"] = str(os.path.join(base_folder, 'src/build/lib'))
import h5py
import hdf5plugin
#plugins_folder = hdf5plugin.PLUGIN_PATH
import numpy as np
import xarray


f_orig = xarray.open_dataset(f'geopotential_pl_small.nc')
file_sz = f'geopotential_pl_small_sz.nc'
file_jp2 = f'geopotential_pl_small_jp2.nc'
f_sz = h5py.File(file_sz, 'r')
f_jp2 = h5py.File(file_jp2, 'r')

data_orig = f_orig['z'].data
data_sz = f_sz['z'][...]
data_jp2 = f_jp2['z'][...]

rmse_sz = np.sqrt(np.square(data_sz - data_orig).mean())
rmse_jp2 = np.sqrt(np.square(data_jp2 - data_orig).mean())
max_sz = np.max(np.abs(data_sz-data_orig))
max_jp2 = np.max(np.abs(data_jp2-data_orig))

size_sz = os.path.getsize(file_sz)
size_jp2 = os.path.getsize(file_jp2)

print(" \tSZ3\tJP2")
print(f"Size\t{size_sz/1024/1024:.2f}MB\t{size_jp2/1024/1024:.2f}MB")
print(f"RMSE\t{rmse_sz:.3f}\t{rmse_jp2:.3f}")
print(f"MAX\t{max_sz:.3f}\t{max_jp2:.3f}")
