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

    def __init__(self, base_cr: float, height: int, width: int, residual_opt: Tuple[str, float], data_dim: int = 2, chunks: Tuple[int, ...] = None):
        assert height > 0 and width > 0
        base_cr = float(base_cr)
        hdf_filter_opts = [int(height), int(width), float_to_uint32(base_cr)]
        self.base_cr = base_cr
        self.height = height
        self.width = width
        self.residual_opt = residual_opt
        if residual_opt is None:
            self.residual_opt = ("fixed_cr", 0)
        residual_type_str, residual_opt_val = self.residual_opt
        self.data_dim = data_dim
        residual_opt_val = float(residual_opt_val)
        if residual_type_str == "fixed_cr":
            residual_type = 0
            hdf_filter_opts.append(residual_type)
        elif residual_type_str == "max_error_target":
            residual_type = 1
            hdf_filter_opts.extend([residual_type, float_to_uint32(residual_opt_val)])
        elif residual_type_str == "relative_error_target":
            residual_type = 2
            hdf_filter_opts.extend([residual_type, float_to_uint32(residual_opt_val)])
        else:
            print(f""""Unknown residual_type {residual_type_str}, has to be one of 
                  'fixed_cr', 'max_error_target' or 'relative_error_target'""")

        self.hdf_filter_opts = tuple(hdf_filter_opts)
        if chunks is None:
            self.chunks = (*[1 for _ in range(self.data_dim - 2)], height, width)
        else:
            self.chunks = tuple(int(n) for n in chunks)
            if len(self.chunks) != data_dim or self.chunks[-2:] != (height, width) or any(n <= 0 for n in self.chunks):
                raise ValueError("chunks must match data_dim and end in (height, width)")
            from math import prod
            if prod(self.chunks[:-2]) > 127:
                raise ValueError("MCT supports at most 127 leading components per chunk")

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

    parser.add_argument('-b', '--base_cr', type=float, default=None, help='target compression ratio; selects FIXED_CR mode, ignored with a warning if -m or -r is supplied')
    parser.add_argument('-H', '--height', type=int, default=721, help='height of the data slice or size of latitude dim')
    parser.add_argument('-W', '--width', type=int, default=1440, help='width of the data slice or size of longitude dim')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('-m', '--max_error_target', default=None, type=float, help='max error target')
    mode.add_argument('-r', '--relative_error_target', default=None, type=float, help='relative error target')
    parser.add_argument('--help-cdo', action='store_true', help='print CDO help')

    args = parser.parse_args()

    base_cr = 200.0 if args.base_cr is None else args.base_cr

    if args.max_error_target is not None:
        residual_opt_val = float(args.max_error_target)
        residual_type = "max_error_target"
    elif args.relative_error_target is not None:
        residual_opt_val = float(args.relative_error_target)
        residual_type = "relative_error_target"
    elif args.base_cr is not None:
        residual_opt_val = 0
        residual_type = "fixed_cr"
    else:
        print('Using default settings: relative error target of 0.01', file=sys.stderr)
        residual_opt_val = 0.01
        residual_type = "relative_error_target"

    if args.base_cr is not None and residual_type != "fixed_cr":
        print("Warning: --base_cr is ignored when an error bound is supplied.", file=sys.stderr)

    ebcc_filter = EBCC_Filter(base_cr=base_cr,
                                    height=args.height, 
                                    width=args.width, 
                                    residual_opt=(residual_type, residual_opt_val))
    
    print("======Configuration======", file=sys.stderr)
    if residual_type == "fixed_cr":
        print(f"base_cr: {base_cr}", file=sys.stderr)
    print(f"HeightxWidth: {args.height}x{args.width}", file=sys.stderr)
    print(f"Residual option: {residual_type}, {residual_opt_val}", file=sys.stderr)

    opts = ",".join([str(opt) for opt in ebcc_filter.hdf_filter_opts])
    opts = f"{EBCC_Filter.FILTER_ID},{opts}"

    if args.help_cdo:
        print(f"Compression using cdo: cdo -b F32 -f nc4 --filter {opts} copy original.nc compressed.nc")
        print(f"Make sure to check chunksize of original.nc divides the tile size {args.height}x{args.width}")

    print(opts)
