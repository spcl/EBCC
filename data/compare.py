import os
import pyproj # Avoids HDF5 plugin error
base_folder = "/home/huanglangwen/Documents/compression-filter"
#os.environ["HDF5_PLUGIN_PATH"] = str(os.path.join(base_folder, 'src/build/lib')) + ":" + str(plugins_folder)
os.environ["HDF5_PLUGIN_PATH"] = str(os.path.join(base_folder, 'src/build/lib'))
import h5py
import hdf5plugin
#plugins_folder = hdf5plugin.PLUGIN_PATH
import numpy as np
import xarray

def stats(data_orig, filename):
    f = h5py.File(filename, 'r')
    data = f['z'][...]
    rmse = np.sqrt(np.square(data - data_orig).mean())
    max_diff = np.max(np.abs(data - data_orig))
    size = os.path.getsize(filename) / 1024 / 1024
    return rmse, max_diff, size


f_orig = xarray.open_dataset(f'geopotential_pl_small.nc')
data_orig = f_orig['z'].data
rmse_jp2, max_jp2, size_jp2 = stats(data_orig, f'geopotential_pl_small_jp2.nc')
rmse_sperr, max_sperr, size_sperr = stats(data_orig, f'geopotential_pl_small_sperr.nc')
rmse_sz, max_sz, size_sz = stats(data_orig, f'geopotential_pl_small_sz.nc')
rmse_sz3, max_sz3, size_sz3 = stats(data_orig, f'geopotential_pl_small_sz3.nc')
rmse_sz3_2, max_sz3_2, size_sz3_2 = stats(data_orig, f'geopotential_pl_small_sz3_2.nc')

print(" \tSZ\tSZ3\tSZ3\tJP2SPWV\tSPERR")
print(f"Size\t{size_sz:.2f}MB\t{size_sz3:.2f}MB\t{size_sz3_2:.2f}MB\t{size_jp2:.2f}MB\t{size_sperr:.2f}MB")
print(f"RMSE\t{rmse_sz:.3f}\t{rmse_sz3:.3f}\t{rmse_sz3_2:.3f}\t{rmse_jp2:.3f}\t{rmse_sperr:.3f}")
print(f"MAX\t{max_sz:.3f}\t{max_sz3:.3f}\t{max_sz3_2:.3f}\t{max_jp2:.3f}\t{max_sperr:.3f}")
