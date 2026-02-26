import argparse
import struct
import sys
from typing import Optional, Tuple
from collections.abc import Mapping


def float_to_uint32(f):
    packed = struct.pack('f', f)
    return struct.unpack('I', packed)[0]


class EBCC_Filter(Mapping):
    FILTER_ID = 308
    FILTER_ID_J2K = 308
    FILTER_ID_JXL = 309

    def __init__(
        self,
        base_param: Optional[float] = None,
        base_cr: Optional[float] = None,
        height: int = 721,
        width: int = 1440,
        residual_opt: Optional[Tuple[str, float]] = None,
        data_dim: int = 2,
        base_compressor: str = "j2k",
        base_distance: Optional[float] = None,
    ):
        assert height > 0 and width > 0

        base_compressor = base_compressor.lower()
        if base_compressor not in {"j2k", "jxl"}:
            raise ValueError("base_compressor must be 'j2k' or 'jxl'")

        if residual_opt is None:
            residual_opt = ("none", 0.0)

        if base_compressor == "j2k":
            if base_distance is not None:
                raise ValueError("base_distance can only be used when base_compressor='jxl'")
            if base_cr is None:
                if base_param is not None:
                    base_cr = float(base_param)
                else:
                    base_cr = 200.0
            elif base_param is not None and float(base_param) != float(base_cr):
                raise ValueError("base_param and base_cr are both set but have different values")
            base_param = float(base_cr)
            self.base_cr = base_param
            self.base_distance = None
            self.filter_id = self.FILTER_ID_J2K
        else:
            if base_cr is not None:
                raise ValueError("base_cr cannot be used when base_compressor='jxl'; use base_distance")
            if base_distance is None:
                if base_param is not None:
                    base_distance = float(base_param)
                else:
                    base_distance = 1.0
            elif base_param is not None and float(base_param) != float(base_distance):
                raise ValueError("base_param and base_distance are both set but have different values")
            base_param = float(base_distance)
            self.base_cr = None
            self.base_distance = base_param
            self.filter_id = self.FILTER_ID_JXL

        hdf_filter_opts = [int(height), int(width), float_to_uint32(base_param)]

        self.base_param = base_param
        self.base_compressor = base_compressor
        self.height = int(height)
        self.width = int(width)
        self.residual_opt = residual_opt
        self.data_dim = data_dim

        residual_type_str, residual_opt_val = residual_opt
        residual_opt_val = float(residual_opt_val)
        if residual_type_str == "none":
            hdf_filter_opts.append(0)
        elif residual_type_str == "max_error_target":
            hdf_filter_opts.extend([1, float_to_uint32(residual_opt_val)])
        elif residual_type_str == "relative_error_target":
            hdf_filter_opts.extend([2, float_to_uint32(residual_opt_val)])
        else:
            raise ValueError(
                "Unknown residual_type; must be one of 'none', 'max_error_target', 'relative_error_target'"
            )

        self.hdf_filter_opts = tuple(hdf_filter_opts)
        self.chunks = (*[1 for _ in range(self.data_dim - 2)], self.height, self.width)

    @property
    def _kwargs(self):
        return {
            'dtype': 'float32',
            'chunks': self.chunks,
            'compression': self.filter_id,
            'compression_opts': self.hdf_filter_opts,
        }

    def __hash__(self):
        return hash((self.filter_id, self.hdf_filter_opts))

    def __len__(self):
        return len(self._kwargs)

    def __iter__(self):
        return iter(self._kwargs)

    def __getitem__(self, item):
        return self._kwargs[item]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', '--base_param', type=float, default=None, help='base parameter (CR for J2K, distance for JXL)')
    parser.add_argument('-b', '--base_cr', type=float, default=None, help='base compression ratio (J2K)')
    parser.add_argument('--base_distance', type=float, default=None, help='base distance (JXL)')
    parser.add_argument('--base_compressor', choices=['j2k', 'jxl'], default='j2k', help='base compressor backend')
    parser.add_argument('-H', '--height', type=int, default=721, help='height of the data slice or size of latitude dim')
    parser.add_argument('-W', '--width', type=int, default=1440, help='width of the data slice or size of longitude dim')
    parser.add_argument('-m', '--max_error_target', default=None, type=float, help='max error target')
    parser.add_argument('-r', '--relative_error_target', default=None, type=float, help='relative error target')
    parser.add_argument('--help-cdo', action='store_true', help='print CDO help')

    args = parser.parse_args()

    if args.max_error_target is not None:
        residual_opt = ("max_error_target", float(args.max_error_target))
    elif args.relative_error_target is not None:
        residual_opt = ("relative_error_target", float(args.relative_error_target))
    else:
        print('Using default settings: relative error target of 0.01', file=sys.stderr)
        residual_opt = ("relative_error_target", 0.01)

    ebcc_filter = EBCC_Filter(
        base_param=args.base_param,
        base_cr=args.base_cr,
        base_distance=args.base_distance,
        base_compressor=args.base_compressor,
        height=args.height,
        width=args.width,
        residual_opt=residual_opt,
    )

    print("======Configuration======", file=sys.stderr)
    print(f"Base compressor: {ebcc_filter.base_compressor}", file=sys.stderr)
    print(f"Base parameter: {ebcc_filter.base_param}", file=sys.stderr)
    print(f"HeightxWidth: {args.height}x{args.width}", file=sys.stderr)
    print(f"Residual option: {residual_opt[0]}, {residual_opt[1]}", file=sys.stderr)

    opts = ",".join([str(opt) for opt in ebcc_filter.hdf_filter_opts])
    opts = f"{ebcc_filter.filter_id},{opts}"

    if args.help_cdo:
        print(f"Compression using cdo: cdo -b F32 -f nc4 --filter {opts} copy original.nc compressed.nc")
        print(f"Make sure to check chunksize of original.nc divides the tile size {args.height}x{args.width}")

    print(opts)
