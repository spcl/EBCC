import os
os.environ["HDF5_PLUGIN_PATH"] = os.path.join(os.getcwd(), 'src')

import h5py
import numpy as np
import matplotlib.pyplot as plt

try:
    os.remove('test.hdf5')
except:
    pass

f = h5py.File('test.hdf5', 'a')

data = np.load('test_data.npy')

# python params_helper.py -b 100 -e 1
opts = (721,1440,1120403456,2,1065353216)

f.create_dataset('compressed', shape=data.shape, dtype=data.dtype, chunks=data.shape, compression=308, compression_opts=opts)

f['compressed'][:] = data

fig, (ax1, ax2) = plt.subplots(1, 2, layout='constrained')
ax1.imshow(data)
ax2.imshow(f['compressed'][:])

fig.savefig('test.pdf', bbox_inches='tight')

original_size = os.path.getsize('test_data.npy')
compressed_size = os.path.getsize('test.hdf5')

print(f'achieved compression ratio of {original_size/compressed_size}')