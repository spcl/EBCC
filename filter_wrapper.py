import argparse
import struct
import os
from typing import Tuple
from collections.abc import Mapping

def double_to_uint32(f):
    packed = struct.pack('d', f)
    return struct.unpack('II', packed)

def float_to_uint32(f):
    packed = struct.pack('f', f)
    return struct.unpack('I', packed)[0]

class JP2SPWV_Filter(Mapping):
    FILTER_ID = 308

    def __init__(self, base_cr: float, height: int, width: int, residual_opt: Tuple[str, float], data_dim: int = 2, filter_path: str = None):
        assert height > 0 and width > 0
        base_cr = float(base_cr)
        hdf_filter_opts = [int(height), int(width), float_to_uint32(base_cr)]
        self.base_cr = base_cr
        self.height = height
        self.width = width
        self.residual_opt = residual_opt
        residual_type_str, residual_opt_val = residual_opt
        self.data_dim = data_dim
        residual_opt_val = float(residual_opt_val)
        if residual_type_str == "quantile_target":
            residual_type = 1
            hdf_filter_opts.extend([residual_type, float_to_uint32(residual_opt_val)])
        elif residual_type_str == "max_error_target":
            residual_type = 2
            hdf_filter_opts.extend([residual_type, float_to_uint32(residual_opt_val)])
        elif residual_type_str == "relative_error_target":
            residual_type = 3
            hdf_filter_opts.extend([residual_type, float_to_uint32(residual_opt_val)])
        elif residual_type_str == "fixed_sparsification":
            residual_type = 4
            q_a, q_b = double_to_uint32(residual_opt_val)
            hdf_filter_opts.extend([residual_type, q_a, q_b])
        else:
            print(f""""Unknown residual_type {residual_type_str}, has to be one of 'quantile_target', 
                  'max_error_target', 'relative_error_target' or 'fixed_sparsification""")

        self.hdf_filter_opts = tuple(hdf_filter_opts)
        self.chunks = (*[1 for _ in range(self.data_dim - 2)], height, width)
        if filter_path is None:
            filter_path = os.path.join(os.path.dirname(__file__), 'src')
        os.environ["HDF5_PLUGIN_PATH"] = filter_path

    # https://github.com/silx-kit/hdf5plugin/blob/main/src/hdf5plugin/_filters.py
    @property
    def _kwargs(self):
        return {
            'chunks': self.chunks,
            'compression': self.FILTER_ID,
            'compression_opts': self.hdf_filter_opts
        }

    def __hash__(self):
        return hash((self.FILTER_ID, self.hdf_filter_opts))

    def __len__(self):
        return len(self._kwargs)

    def __iter__(self):
        return iter(self._kwargs)

    def __getitem__(self, item):
        return self._kwargs[item]

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-b', '--base_cr', type=str, default=200, help='base compression ratio')
    parser.add_argument('-h', '--height', type=int, default=721, help='height of the data slice or size of latitude dim')
    parser.add_argument('-w', '--width', type=int, default=1440, help='width of the data slice or size of longitude dim')
    parser.add_argument('-m', '--max_error_target', default=None, type=float)
    parser.add_argument('-r', '--relative_error_target', default=None, type=float)
    parser.add_argument('-q', '--quantile_target', default=None, type=float)
    parser.add_argument('-s', '--fixed_sparsification', default=None, type=float)

    args = parser.parse_args()

    residual_type = 0

    base_cr = float(args.base_cr)

    if args.quantile_target:
        residual_opt_val = float(args.quantile_target)
        residual_type = "quantile_target"
    elif args.max_error_target:
        residual_opt_val = float(args.max_error_target)
        residual_type = "max_error_target"
    elif args.relative_error_target:
        residual_opt_val = float(args.relative_error_target)
        residual_type = "relative_error_target"
    elif args.sparsification:
        residual_opt_val = float(args.fixed_sparsification)
        residual_type = "fixed_sparsification"

    jp2spwv_filter = JP2SPWV_Filter(base_cr=args.base_cr, 
                                    height=args.height, 
                                    width=args.width, 
                                    residual_opt=(residual_type, residual_opt_val))


    print(dict(jp2spwv_filter))


