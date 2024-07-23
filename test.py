import os
current_folder = os.path.dirname(os.path.realpath(__file__))
os.environ["HDF5_PLUGIN_PATH"] = os.path.join(current_folder, 'src/build/lib')

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from filter_wrapper import JP2SPWV_Filter

try:
    os.remove(f'{current_folder}/test.hdf5')
except:
    pass

f = h5py.File(f'{current_folder}/test.hdf5', 'a')

data = np.load(f'{current_folder}/test_data.npy')

jp2spwv_filter = JP2SPWV_Filter(
    base_cr=100, # base compression ratio
    height=data.shape[0],  # height of each 2D data chunk
    width=data.shape[1],  # width of each 2D data chunk
    data_dim=len(data.shape), # data dimension, required to specify the HDF5 chunk shape
    residual_opt=("relative_error_target", 0.0019),# specify the relative error target to be 0.0019
    filter_path=os.path.join(current_folder, 'src')) # directory to the compiled HDF5 filter plugin
    # other possible residual_opt can be
    # `("max_error_target", xxx)` : the max_error does not exceed the specified value
    # `("quantile_target", xxx)` : specifies the quantile used to sparsify the wavelet transformed residual
    # `("fixed_sparsification", xxx)`: specify a fixed sparsification ratio for the sparse wavelet compression

print(dict(jp2spwv_filter))
f.create_dataset('compressed', shape=data.shape, dtype=data.dtype, **jp2spwv_filter)

f['compressed'][:] = data
uncompressed = f['compressed'][:]

# check if error target is correctly enforced
max_error = np.max(np.abs(data - uncompressed) / np.abs(data))
print('achieved max relative error:', max_error)

fig, (ax1, ax2) = plt.subplots(1, 2, layout='constrained')
ax1.imshow(data)
ax2.imshow(uncompressed)

fig.savefig(f'{current_folder}/test.pdf', bbox_inches='tight')

original_size = os.path.getsize(f'{current_folder}/test_data.npy')
compressed_size = os.path.getsize(f'{current_folder}/test.hdf5')

print(f'achieved compression ratio of {original_size/compressed_size}')
