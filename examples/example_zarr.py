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

width = 32
height = 32
data = np.zeros((1, 1, height, width), dtype=np.float32)
for i in range(height):
    for j in range(width):
        data[0, 0, i, j] = i + j
atol = 1e-2

for compressor in ["j2k", "jxl"]:
    print(f"\n=== Testing base compressor: {compressor} ===")
    store_path = f'{data_dir}/test_{compressor}.zarr'

    try:
        shutil.rmtree(store_path)
    except Exception:
        pass

    if compressor == "j2k":
        ebcc_filter = EBCC_Filter(
            base_cr=2,
            height=height,
            width=width,
            residual_opt=("max_error_target", atol),
            base_compressor="j2k",
        )
    else:
        ebcc_filter = EBCC_Filter(
            base_distance=1.0,
            height=height,
            width=width,
            residual_opt=("max_error_target", atol),
            base_compressor="jxl",
        )

    print(dict(ebcc_filter))

    zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)

    z = zarr.create_array(
        store=store_path,
        data=data,
        chunks=(1, 1, height, width),
        serializer=AnyNumcodecsArrayBytesCodec(zarr_filter),
        compressors=None,
        overwrite=True,
    )

    decoded_zarr = z[:]
    assert np.allclose(data, decoded_zarr, atol=atol), \
        f"Decoded Zarr data does not match original data for {compressor}"

    data_2d = data[0, 0]
    decoded_2d = decoded_zarr[0, 0]

    data_range = np.max(data_2d) - np.min(data_2d)
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
    fig.savefig(f'{data_dir}/test_zarr_{compressor}.pdf', bbox_inches='tight')
    plt.close(fig)

    original_size = data.nbytes
    compressed_size = sum(
        os.path.getsize(os.path.join(dirpath, filename))
        for dirpath, dirnames, filenames in os.walk(store_path)
        for filename in filenames
    )
    print(f'achieved compression ratio of {original_size / compressed_size}')
