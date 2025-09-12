import os
from ebcc import EBCC_FILTER_DIR
os.environ["HDF5_PLUGIN_PATH"] = EBCC_FILTER_DIR

import h5py
import numpy as np
from ebcc.filter_wrapper import EBCC_Filter
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='EBCC HDF5 compression example')
    parser.add_argument('--data', default=None, help='Path to input .npy data file (default: data/test_data.npy relative to this script)')
    parser.add_argument('--output-hdf5', default=None, help='Output HDF5 file path (default: data/test.hdf5 relative to this script)')
    parser.add_argument('--output-pdf', default=None, help='Output PDF figure path (default: data/test.pdf relative to this script)')
    parser.add_argument('--base-cr', type=int, default=100, help='Base compression ratio (default: 100)')
    parser.add_argument('--residual-mode', choices=['none', 'relative_error_target', 'max_error_target'], default='relative_error_target', help='Residual compression mode')
    parser.add_argument('--residual-target', type=float, default=0.009, help='Residual error target value')
    parser.add_argument('--plot', action='store_true', help='Enable plotting and PDF output')
    return parser.parse_args()


def main():
    args = parse_args()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, '../data')

    npy_path = args.data or f'{data_dir}/test_data.npy'
    hdf5_path = args.output_hdf5 or f'{data_dir}/test.hdf5'
    pdf_path = args.output_pdf or f'{data_dir}/test.pdf'

    try:
        os.remove(hdf5_path)
    except Exception:
        pass

    f = h5py.File(hdf5_path, 'a')

    data = np.load(npy_path)
    # data = np.zeros((721, 1440)) + 99

    if args.residual_mode == 'none':
        residual_opt = None
    else:
        residual_opt = (args.residual_mode, args.residual_target)

    ebcc_filter = EBCC_Filter(
        base_cr=args.base_cr,  # base compression ratio
        height=data.shape[0],  # height of each 2D data chunk
        width=data.shape[1],  # width of each 2D data chunk
        data_dim=len(data.shape),  # data dimension, required to specify the HDF5 chunk shape
        residual_opt=residual_opt,  # residual compression option
    )

    print(dict(ebcc_filter))
    f.create_dataset('compressed', shape=data.shape, **ebcc_filter)

    f['compressed'][:] = data
    uncompressed = f['compressed'][:]

    # check if error target is correctly enforced
    data_range = (np.max(data) - np.min(data))
    max_error = np.max(np.abs(data - uncompressed))
    if data_range > 0:
        rel_error = max_error / data_range
        print('achieved max relative error:', rel_error)
    else:
        print('achieved max absolute error:', max_error)

    if args.plot:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
        fig, (ax1, ax2) = plt.subplots(1, 2, layout='constrained')
        ax1.imshow(data)
        ax2.imshow(uncompressed)
        fig.savefig(pdf_path, bbox_inches='tight')

    original_size = os.path.getsize(npy_path)
    compressed_size = os.path.getsize(hdf5_path)

    print(f'achieved compression ratio of {original_size / compressed_size}')


if __name__ == '__main__':
    main()
