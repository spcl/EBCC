import zarr
import numpy as np
from zarr_any_numcodecs import AnyNumcodecsArrayBytesCodec
from ebcc.zarr_filter import EBCCZarrFilter
from ebcc.filter_wrapper import EBCC_Filter
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(current_dir, '../data')

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

zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)

encoded = zarr_filter.encode(data)

decoded = zarr_filter.decode(encoded)
decoded = decoded.reshape(data.shape)

assert np.allclose(data, decoded, atol=atol), "Decoded data does not match original data"

z = zarr.create_array(
    store=f"{data_dir}/test.zarr",
    data=data,
    chunks="auto",
    serializer=AnyNumcodecsArrayBytesCodec(zarr_filter),
    compressors=None,
    overwrite=True,
    )

decoded_zarr = z[:]
assert np.allclose(data, decoded_zarr, atol=atol), "Decoded Zarr data does not match original data"