import os
from ebcc import EBCC_FILTER_DIR
os.environ["HDF5_PLUGIN_PATH"] = EBCC_FILTER_DIR

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from ebcc.filter_wrapper import EBCC_Filter

current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(current_dir, '../data')

try:
    os.remove(f'{data_dir}/test.hdf5')
except:
    pass

f = h5py.File(f'{data_dir}/test.hdf5', 'a')

data = np.load(f'{data_dir}/test_data.npy')
#data = np.zeros((721, 1440)) + 99

ebcc_filter = EBCC_Filter(
    base_cr=100, # base compression ratio
    height=data.shape[0],  # height of each 2D data chunk
    width=data.shape[1],  # width of each 2D data chunk
    data_dim=len(data.shape), # data dimension, required to specify the HDF5 chunk shape
    residual_opt=("relative_error_target", 0.009),# specify the relative error target to be 0.0019
    )
    # other possible residual_opt can be
    # `None` : no residual compression
    # `("max_error_target", xxx)` : the max_error does not exceed the specified value
    # `("quantile_target", xxx)` : specifies the quantile used to sparsify the wavelet transformed residual
    # `("fixed_sparsification", xxx)`: specify a fixed sparsification ratio for the sparse wavelet compression

print(dict(ebcc_filter))
f.create_dataset('compressed', shape=data.shape,  **ebcc_filter)

f['compressed'][:] = data
uncompressed = f['compressed'][:]

# check if error target is correctly enforced
#max_error = np.max(np.abs(data - uncompressed) / np.abs(data))
data_range = (np.max(data) - np.min(data))
max_error = np.max(np.abs(data - uncompressed))
if data_range > 0:
    rel_error = max_error / data_range
    print('achieved max relative error:', rel_error)
else:
    print('achieved max absolute error:', max_error)

fig, (ax1, ax2) = plt.subplots(1, 2, layout='constrained')
ax1.imshow(data)
ax2.imshow(uncompressed)

fig.savefig(f'{data_dir}/test.pdf', bbox_inches='tight')

original_size = os.path.getsize(f'{data_dir}/test_data.npy')
compressed_size = os.path.getsize(f'{data_dir}/test.hdf5')

print(f'achieved compression ratio of {original_size/compressed_size}')
