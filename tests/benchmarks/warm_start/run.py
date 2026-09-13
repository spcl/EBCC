"""Measure isolated warm-start prototypes and require byte-identical output."""

import argparse
import ctypes as c
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
import xarray as xr

from ebcc.zarr_filter import CodecConfigT


def load_probe(path):
    lib = c.CDLL(str(path.resolve()))
    lib.probe_create.argtypes = [c.c_void_p, c.c_uint, c.c_uint]
    lib.probe_create.restype = c.c_void_p
    lib.probe_encode.argtypes = [c.c_void_p, c.c_float, c.c_float, c.c_int]
    lib.probe_encode.restype = c.c_int
    lib.probe_data.argtypes = [c.c_void_p]
    lib.probe_data.restype = c.c_void_p
    lib.probe_size.argtypes = [c.c_void_p]
    lib.probe_size.restype = c.c_size_t
    lib.probe_destroy.argtypes = [c.c_void_p]
    lib.probe_destroy.restype = None
    return lib


def load_codec(path):
    lib = c.CDLL(str(path.resolve()))
    lib.ebcc_encode.argtypes = [c.c_void_p, c.POINTER(CodecConfigT), c.POINTER(c.c_void_p)]
    lib.ebcc_encode.restype = c.c_size_t
    lib.ebcc_decode.argtypes = [c.c_void_p, c.c_size_t, c.POINTER(c.c_void_p)]
    lib.ebcc_decode.restype = c.c_size_t
    lib.free_buffer.argtypes = [c.c_void_p]
    lib.free_buffer.restype = None
    return lib


def clocked(fn):
    wall, cpu = time.perf_counter(), time.process_time()
    result = fn()
    return result, time.perf_counter() - wall, time.process_time() - cpu


def record(output, **row):
    with output.open("a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)


def micro(lib, data, meta, args):
    # Non-monotone requests and mode switches exercise cache reuse in both
    # directions. Rates are OpenJPEG rates, half EBCC's float32-relative CR.
    requests = [("psnr", q) for q in [96, 94, 95, 78, 58]]
    requests += [("cr", r) for r in [15, 50, 500]]
    requests += [("psnr", q) for q in [96, 58]]
    scaled = ((data - data.min()) / np.ptp(data) * np.iinfo(np.uint16).max).astype(np.uint16)
    for repeat in range(args.repeats):
        sessions = [lib.probe_create(scaled.ctypes.data, data.shape[1], data.shape[0]) for _ in range(2)]
        assert all(sessions)
        try:
            for index, (mode, value) in enumerate(requests):
                outputs = {}
                metrics = {}
                for warm in ([0, 1] if repeat % 2 == 0 else [1, 0]):
                    session = sessions[warm]
                    psnr, rate = (value, 0) if mode == "psnr" else (0, value)
                    ok, wall, cpu = clocked(lambda: lib.probe_encode(session, psnr, rate, warm))
                    assert ok
                    outputs[warm] = c.string_at(lib.probe_data(session), lib.probe_size(session))
                    metrics[warm] = dict(wall=wall, cpu=cpu)
                assert outputs[0] == outputs[1], (meta, mode, value, index)
                record(args.output, kind="j2k", **meta, repeat=repeat, index=index,
                       mode=mode, value=value, bytes=len(outputs[0]), identical=True,
                       cold=metrics[0], warm=metrics[1],
                       sha256=hashlib.sha256(outputs[0]).hexdigest())
        finally:
            for session in sessions:
                lib.probe_destroy(session)


def codec(lib, data, meta, args):
    config = CodecConfigT()
    config.dims[:] = (1, *data.shape)
    config.base_cr = 30
    config.residual_compression_type = 2  # RELATIVE_ERROR
    for repeat in range(args.repeats):
        for mode in ["cr", "psnr"]:
            os.environ["EBCC_BASE_SEARCH"] = mode
            for target in [0.0001, 0.001, 0.01]:
                config.error = target
                outputs, decoded_outputs, metrics = {}, {}, {}
                for warm in ([0, 1] if repeat % 2 == 0 else [1, 0]):
                    os.environ["EBCC_WARM_START"] = str(warm)
                    out = c.c_void_p()
                    size, wall, cpu = clocked(lambda: lib.ebcc_encode(data.ctypes.data, c.byref(config), c.byref(out)))
                    assert size and out.value
                    try:
                        outputs[warm] = c.string_at(out, size)
                        decoded = c.c_void_p()
                        count, read_wall, read_cpu = clocked(lambda: lib.ebcc_decode(out, size, c.byref(decoded)))
                        assert count == data.size and decoded.value
                        try:
                            values = np.ctypeslib.as_array(c.cast(decoded, c.POINTER(c.c_float)), shape=(data.size,)).copy().reshape(data.shape)
                        finally:
                            lib.free_buffer(decoded)
                    finally:
                        lib.free_buffer(out)
                    error = float(np.max(np.abs(data - values)))
                    assert error <= target * float(np.ptp(data)), (meta, mode, target, warm, error)
                    decoded_outputs[warm] = values
                    metrics[warm] = dict(wall=wall, cpu=cpu, decode_wall=read_wall,
                                         decode_cpu=read_cpu, max_error=error)
                assert outputs[0] == outputs[1], (meta, mode, target)
                assert np.array_equal(decoded_outputs[0], decoded_outputs[1])
                record(args.output, kind="ebcc", **meta, repeat=repeat, mode=mode,
                       target=target, scale=args.scale, bytes=len(outputs[0]),
                       input_bytes=data.nbytes, identical=True, passed=True,
                       cold=metrics[0], warm=metrics[1],
                       sha256=hashlib.sha256(outputs[0]).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--j2k-library", type=Path)
    parser.add_argument("--skip-micro", action="store_true",
                        help="Compare full EBCC cold/warm paths in the integrated library")
    parser.add_argument("--ebcc-library", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=Path("data/era5_pl_sample.nc"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--variables", nargs="+", default=["z", "t", "q", "u", "v", "w"])
    parser.add_argument("--levels", nargs="+", type=float, default=[500])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--scale", type=float, default=6)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a new output file to avoid mixing experiments")
    if args.repeats < 1:
        parser.error("repeats must be positive")
    if not args.skip_micro and args.j2k_library is None:
        parser.error("--j2k-library is required unless --skip-micro is used")
    os.environ["EBCC_LOG_LEVEL"] = "3"
    os.environ["EBCC_PSNR_ERROR_SCALE"] = str(args.scale)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    j2k = None if args.skip_micro else load_probe(args.j2k_library)
    ebcc = load_codec(args.ebcc_library)
    with xr.open_dataset(args.input) as dataset:
        for variable in args.variables:
            for level in args.levels:
                data = np.ascontiguousarray(dataset[variable].isel(valid_time=0).sel(pressure_level=level).values, dtype=np.float32)
                assert np.isfinite(data).all() and np.ptp(data) > 0
                meta = dict(variable=variable, pressure_level=level, time_index=0, shape=data.shape)
                if j2k is not None:
                    micro(j2k, data, meta, args)
                codec(ebcc, data, meta, args)


if __name__ == "__main__":
    main()
