import argparse
import struct
import os
import sys
from typing import Tuple
from collections.abc import Mapping

def double_to_uint32(f):
    packed = struct.pack('d', f)
    return struct.unpack('II', packed)

def float_to_uint32(f):
    packed = struct.pack('f', f)
    return struct.unpack('I', packed)[0]

class EBCC_Filter(Mapping):
    FILTER_ID = 308

    def __init__(self, base_cr: float, height: int, width: int, residual_opt: Tuple[str, float], data_dim: int = 2):
        assert height > 0 and width > 0
        base_cr = float(base_cr)
        hdf_filter_opts = [int(height), int(width), float_to_uint32(base_cr)]
        self.base_cr = base_cr
        self.height = height
        self.width = width
        self.residual_opt = residual_opt
        if residual_opt is None:
            self.residual_opt = ("none", 0)
        residual_type_str, residual_opt_val = self.residual_opt
        self.data_dim = data_dim
        residual_opt_val = float(residual_opt_val)
        if residual_type_str == "none":
            residual_type = 0
            hdf_filter_opts.append(residual_type)
        elif residual_type_str == "quantile_target":
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

    # https://github.com/silx-kit/hdf5plugin/blob/main/src/hdf5plugin/_filters.py
    @property
    def _kwargs(self):
        return {
            'dtype': 'float32',
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
    parser.add_argument('-H', '--height', type=int, default=721, help='height of the data slice or size of latitude dim')
    parser.add_argument('-W', '--width', type=int, default=1440, help='width of the data slice or size of longitude dim')
    parser.add_argument('-m', '--max_error_target', default=None, type=float, help='max error target')
    parser.add_argument('-r', '--relative_error_target', default=None, type=float, help='relative error target')
    parser.add_argument('-q', '--quantile_target', default=None, type=float, help='[DEPRECATED!] quantile target')
    parser.add_argument('-s', '--fixed_sparsification', default=None, type=float, help='[DEPRECATED!] fixed sparsification')
    parser.add_argument('--help-cdo', action='store_true', help='print CDO help')

    args = parser.parse_args()

    residual_type = 0

    base_cr = float(args.base_cr)

    if args.quantile_target:
        residual_opt_val = float(args.quantile_target)
        residual_type = "quantile_target"
        print("Quantial target option is deprecated!", file=sys.stderr)
        exit(1)
    elif args.max_error_target:
        residual_opt_val = float(args.max_error_target)
        residual_type = "max_error_target"
    elif args.relative_error_target:
        residual_opt_val = float(args.relative_error_target)
        residual_type = "relative_error_target"
    elif args.fixed_sparsification:
        residual_opt_val = float(args.fixed_sparsification)
        residual_type = "fixed_sparsification"
        print("Fixed sparsification option is deprecated!", file=sys.stderr)
        exit(1)
    else:
        print('Using default settings: relative error target of 0.01', file=sys.stderr)
        residual_opt_val = 0.01
        residual_type = "relative_error_target"

    ebcc_filter = EBCC_Filter(base_cr=args.base_cr, 
                                    height=args.height, 
                                    width=args.width, 
                                    residual_opt=(residual_type, residual_opt_val))
    
    print("======Configuration======", file=sys.stderr)
    print(f"Base compression ratio: {args.base_cr}", file=sys.stderr)
    print(f"HeightxWidth: {args.height}x{args.width}", file=sys.stderr)
    print(f"Residual option: {residual_type}, {residual_opt_val}", file=sys.stderr)

    opts = ",".join([str(opt) for opt in ebcc_filter.hdf_filter_opts])
    opts = f"{EBCC_Filter.FILTER_ID},{opts}"

    if args.help_cdo:
        print(f"Compression using cdo: cdo -b F32 -f nc4 --filter {opts} copy original.nc compressed.nc")
        print(f"Make sure to check chunksize of original.nc divides the tile size {args.height}x{args.width}")

    print(opts)
