import os
os.environ["HDF5_PLUGIN_PATH"] = os.path.join(os.getcwd(), 'src')

import h5py
import numpy as np
import matplotlib.pyplot as plt
from filter_wrapper import JP2SPWV_Filter

try:
    os.remove('test.hdf5')
except:
    pass

f = h5py.File('test.hdf5', 'a')

data = np.load('test_data.npy')

jp2spwv_filter = JP2SPWV_Filter(
    base_cr=100, # base compression ratio
    height=data.shape[0],  # height of each 2D data chunk
    width=data.shape[1],  # width of each 2D data chunk
    data_dim=len(data.shape), # data dimension, required to specify the HDF5 chunk shape
    residual_opt=("max_error_target", 1.0),# specify the max error target to be 1.0
    filter_path=os.path.join(os.path.dirname(__file__), 'src')) # directory to the compiled HDF5 filter plugin
    # other possible residual_opt can be
    # `("quantile_error_target", xxx)` : max_error does not exceed the specified quantile value calculated from the compression error with only base compression method
    # `("fixed_sparsification", xxx)`: specify a fixed sparsification ratio for the sparse wavelet compression

print(dict(jp2spwv_filter))
f.create_dataset('compressed', shape=data.shape, dtype=data.dtype, **jp2spwv_filter)

f['compressed'][:] = data

fig, (ax1, ax2) = plt.subplots(1, 2, layout='constrained')
ax1.imshow(data)
ax2.imshow(f['compressed'][:])

fig.savefig('test.pdf', bbox_inches='tight')

original_size = os.path.getsize('test_data.npy')
compressed_size = os.path.getsize('test.hdf5')

print(f'achieved compression ratio of {original_size/compressed_size}')
