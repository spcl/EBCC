import argparse
import struct
import sys

def double_to_uint32(f):
    packed = struct.pack('d', f)
    return struct.unpack('II', packed)

def float_to_uint32(f):
    packed = struct.pack('f', f)
    return struct.unpack('I', packed)[0]

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--base', help='base compression ratio', required=True)
parser.add_argument('-g', '--grid', help='comma separated grid height and width value (e.g. 721,1440)')
parser.add_argument('-e', '--max_error_target')
parser.add_argument('-q', '--quantile')
parser.add_argument('-s', '--sparsification')

args = parser.parse_args()

height, width = 721, 1440
base_cr = 0
residual_type = 0

base_cr = float(args.base)

if args.grid:
    try:
        height, width = args.grid.split(',')
    except:
        print('-g --grid should be used with comma separated values (e.g. 721,1440)', file=sys.stderr)
        sys.exit(-1)


if args.quantile:
    quantile = float(args.quantile)
    residual_type = 3
elif args.max_error_target:
    max_error_target = float(args.max_error_target)
    residual_type = 2
elif args.sparsification:
    sparsification = float(args.sparsification)
    residual_type = 1


out_string = f'{height},{width},{float_to_uint32(base_cr)},{residual_type}'

if residual_type == 1:
    out_string += f',{float_to_uint32(sparsification)}'
elif residual_type == 2:
    out_string += f',{float_to_uint32(max_error_target)}'
elif residual_type == 3:
    q_a, q_b = double_to_uint32(quantile)
    out_string += f',{q_a},{q_b}'

print(out_string)
