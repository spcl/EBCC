import zarr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from zarr_any_numcodecs import AnyNumcodecsArrayBytesCodec
from ebcc.zarr_filter import EBCCZarrFilter
from ebcc.filter_wrapper import EBCC_Filter
import os
import shutil

current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(current_dir, '../data')

# Remove existing zarr store if it exists
try:
    shutil.rmtree(f'{data_dir}/test.zarr')
except:
    pass

width = 32
height = 32
data = np.zeros((1, 1, height, width), dtype=np.float32)
for i in range(height):
    for j in range(width):
        data[0, 0, i, j] = i + j
atol = 1e-2
ebcc_filter = EBCC_Filter(base_cr=2, 
                        height=height, 
                        width=width, 
                        residual_opt=("max_error_target", atol))

print(dict(ebcc_filter))

zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)

encoded = zarr_filter.encode(data)

decoded = zarr_filter.decode(encoded)
decoded = decoded.reshape(data.shape)

assert np.allclose(data, decoded, atol=atol), "Decoded data does not match original data"

z = zarr.create_array(
    store=f"{data_dir}/test.zarr",
    data=data,
    chunks=(1, 1, height, width),
    serializer=AnyNumcodecsArrayBytesCodec(zarr_filter),
    compressors=None,
    overwrite=True,
    )

decoded_zarr = z[:]
assert np.allclose(data, decoded_zarr, atol=atol), "Decoded Zarr data does not match original data"

# check if error target is correctly enforced
data_2d = data[0, 0]  # extract 2D data for visualization and error calculation
decoded_2d = decoded_zarr[0, 0]

data_range = (np.max(data_2d) - np.min(data_2d))
max_error = np.max(np.abs(data_2d - decoded_2d))
if data_range > 0:
    rel_error = max_error / data_range
    print('achieved max relative error:', rel_error)
else:
    print('achieved max absolute error:', max_error)

fig, (ax1, ax2) = plt.subplots(1, 2, layout='constrained')
ax1.imshow(data_2d)
ax1.set_title('Original')
ax2.imshow(decoded_2d)
ax2.set_title('Decompressed')

fig.savefig(f'{data_dir}/test_zarr.pdf', bbox_inches='tight')

# Calculate compression ratio (approximate since Zarr uses directory structure)
original_size = data.nbytes
compressed_size = sum(
    os.path.getsize(os.path.join(dirpath, filename))
    for dirpath, dirnames, filenames in os.walk(f'{data_dir}/test.zarr')
    for filename in filenames
)

print(f'achieved compression ratio of {original_size/compressed_size}')