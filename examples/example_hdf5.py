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
    parser.add_argument('--output-hdf5', default=None, help='Output HDF5 file path (default: data/test_<compressor>.hdf5 relative to this script)')
    parser.add_argument('--output-pdf', default=None, help='Output PDF figure path (default: data/test_<compressor>.pdf relative to this script)')
    parser.add_argument('--base-cr', type=int, default=100, help='Base compression ratio for J2K (default: 100)')
    parser.add_argument('--base-distance', type=float, default=3.0, help='Base distance for JXL (default: 1.0)')
    parser.add_argument('--base-compressor', choices=['j2k', 'jxl', 'both'], default='both',
                        help='Base compressor backend to use (default: both)')
    parser.add_argument('--residual-mode', choices=['none', 'relative_error_target', 'max_error_target'], default='relative_error_target', help='Residual compression mode')
    parser.add_argument('--residual-target', type=float, default=0.009, help='Residual error target value')
    parser.add_argument('--plot', action='store_true', help='Enable plotting and PDF output')
    return parser.parse_args()


def run_compression(compressor, data, npy_path, hdf5_path, pdf_path, args):
    print(f"\n=== Testing base compressor: {compressor} ===")

    try:
        os.remove(hdf5_path)
    except Exception:
        pass

    if args.residual_mode == 'none':
        residual_opt = None
    else:
        residual_opt = (args.residual_mode, args.residual_target)

    if compressor == "j2k":
        ebcc_filter = EBCC_Filter(
            base_cr=args.base_cr,
            height=data.shape[0],
            width=data.shape[1],
            data_dim=len(data.shape),
            residual_opt=residual_opt,
            base_compressor="j2k",
        )
    else:
        ebcc_filter = EBCC_Filter(
            base_distance=args.base_distance,
            height=data.shape[0],
            width=data.shape[1],
            data_dim=len(data.shape),
            residual_opt=residual_opt,
            base_compressor="jxl",
        )

    print(dict(ebcc_filter))

    f = h5py.File(hdf5_path, 'a')
    f.create_dataset('compressed', shape=data.shape, **ebcc_filter)
    f['compressed'][:] = data
    uncompressed = f['compressed'][:]
    f.close()

    data_range = np.max(data) - np.min(data)
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
        plt.close(fig)

    original_size = os.path.getsize(npy_path)
    compressed_size = os.path.getsize(hdf5_path)
    print(f'achieved compression ratio of {original_size / compressed_size}')


def main():
    args = parse_args()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, '../data')

    npy_path = args.data or f'{data_dir}/test_data.npy'
    data = np.load(npy_path)

    compressors = ["j2k", "jxl"] if args.base_compressor == "both" else [args.base_compressor]

    for compressor in compressors:
        hdf5_path = args.output_hdf5 or f'{data_dir}/test_{compressor}.hdf5'
        pdf_path = args.output_pdf or f'{data_dir}/test_{compressor}.pdf'
        run_compression(compressor, data, npy_path, hdf5_path, pdf_path, args)


if __name__ == '__main__':
    main()
