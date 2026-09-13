"""Compare CR and experimental PSNR searches on identical float32 input.

Run with EBCC_LOG_LEVEL=0 and redirect stderr to retain iteration traces.
Write timing includes HDF5 creation, compression, and close; read timing
includes opening and decompressing the resulting file.
"""

import argparse
from contextlib import contextmanager, nullcontext
import json
import os
from pathlib import Path
import random
import re
import sys
import time

from ebcc import EBCC_FILTER_DIR

os.environ["HDF5_PLUGIN_PATH"] = EBCC_FILTER_DIR

import h5py
import numpy as np
import xarray as xr
from ebcc.filter_wrapper import EBCC_Filter


@contextmanager
def capture_stderr(path):
    """Capture native codec traces as well as Python stderr."""
    sys.stderr.flush()
    saved = os.dup(2)
    try:
        with path.open("w") as log:
            os.dup2(log.fileno(), 2)
            yield
    finally:
        sys.stderr.flush()
        os.dup2(saved, 2)
        os.close(saved)


def trace_counts(path):
    counts = dict(base_trials=0, base_searches=0, final_reencodes=0,
                  residual_trials=0, residual_searches=0, residual_kept=0,
                  lossless_fallbacks=0, warm_preparations=0, warm_reuses=0)
    for line in path.read_text().splitlines():
        match = re.search(r"ebcc_(?:codec|j2k)\.c:\d+: (.*)", line)
        if not match:
            continue
        message = match[1]
        if message == "Base warm start: prepare":
            counts["warm_preparations"] += 1
        elif message == "Base warm start: reuse":
            counts["warm_reuses"] += 1
        elif message.startswith("PSNR search start:"):
            counts["base_searches"] += 1
        elif message.startswith("PSNR trial:"):
            counts["base_trials"] += 1
        elif message.startswith("current_cr:") and "cr_lo:" not in message:
            counts["base_searches"] += 1
            counts["final_reencodes"] += 1
        elif message.startswith(("cr_lo:", "cr_hi:")) and "jp2_length:" in message:
            counts["base_trials"] += 1
        elif message.startswith("current_cr:") and "cr_lo:" in message:
            counts["base_trials"] += 1
        elif "exceeds 1000, while base compression still feasible" in message:
            counts["final_reencodes"] -= 1
        elif message.startswith("trunc_lo:"):
            key = "residual_searches" if "optimization_error_target:" in message else "residual_trials"
            counts[key] += 1
        elif message.startswith("coeffs_size:"):
            counts["residual_kept"] += int(re.search(r"coeffs_size: (\d+)", message)[1]) > 0
        elif "PSNR search found no feasible" in message:
            counts["lossless_fallbacks"] += 1
    counts["base_encode_calls"] = (counts["base_trials"] + counts["base_searches"]
                                   + counts["final_reencodes"] + counts["lossless_fallbacks"])
    return counts


def run_case(data, variable, target, mode, scale, repeat, args):
    os.environ["EBCC_BASE_SEARCH"] = mode
    os.environ["EBCC_PSNR_ERROR_SCALE"] = str(scale or args.scale)
    tag = f"{variable}-{mode}-rel{target}-scale{scale}-repeat{repeat}"
    print(f"BEGIN {tag}", file=sys.stderr, flush=True)
    output = args.output_dir / f"{tag}.h5"
    trace = args.output_dir / f"{tag}.log"
    options = EBCC_Filter(base_cr=30, height=data.shape[-2],
                          width=data.shape[-1], data_dim=data.ndim,
                          residual_opt=("relative_error_target", target))
    with capture_stderr(trace) if args.trace else nullcontext():
        start = time.perf_counter()
        cpu_start = time.process_time()
        with h5py.File(output, "w") as f:
            f.create_dataset(variable, shape=data.shape, **options)
            f[variable][...] = data
        write_seconds = time.perf_counter() - start
        write_cpu_seconds = time.process_time() - cpu_start
        start = time.perf_counter()
        cpu_start = time.process_time()
        with h5py.File(output, "r") as f:
            decoded = f[variable][...]
        read_seconds = time.perf_counter() - start
        read_cpu_seconds = time.process_time() - cpu_start
    ranges = np.ptp(data, axis=(-2, -1))
    errors = np.max(np.abs(data - decoded), axis=(-2, -1))
    relative = np.divide(errors, ranges, out=np.zeros_like(errors), where=ranges != 0)
    relative = np.where((ranges == 0) & (errors != 0), np.inf, relative)
    passed = bool(np.all(errors <= target * ranges))
    result = dict(variable=variable, mode=mode, target=target, scale=scale, repeat=repeat,
                  time_indices=args.time_indices, pressure_levels=args.pressure_levels,
                  shape=data.shape,
                  input_bytes=data.nbytes, file_bytes=output.stat().st_size,
                  write_seconds=write_seconds, read_seconds=read_seconds,
                  write_cpu_seconds=write_cpu_seconds, read_cpu_seconds=read_cpu_seconds,
                  write_mb_s=data.nbytes / 1e6 / write_seconds,
                  total_seconds=write_seconds + read_seconds,
                  compression_ratio=data.nbytes / output.stat().st_size,
                  max_chunk_relative_error=float(relative.max()), passed=passed)
    if args.trace:
        result.update(trace_counts(trace))
    with (args.output_dir / "results.jsonl").open("a") as f:
        f.write(json.dumps(result) + "\n")
    print(json.dumps(result), flush=True)
    print(f"END {tag}", file=sys.stderr, flush=True)
    if not passed and not args.keep_going:
        raise AssertionError(f"Error bound violated for {tag}: {result}")
    return passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/geopotential_pl_small.nc"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--modes", nargs="+", choices=("cr", "psnr"), default=["cr", "psnr"])
    parser.add_argument("--targets", nargs="+", type=float, default=[0.01, 0.1])
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--scale", type=float, default=8.0)
    parser.add_argument("--scales", nargs="+", type=float, help="Sweep multiple initial-guess scales")
    parser.add_argument("--variables", nargs="+", default=["z"])
    parser.add_argument("--time-indices", nargs="+", type=int, help="Select indices along the first non-spatial dimension")
    parser.add_argument("--pressure-levels", nargs="+", type=float, help="Select pressure_level coordinates in hPa")
    parser.add_argument("--randomize", action="store_true", help="Shuffle mode/scale order with a fixed seed")
    parser.add_argument("--trace", action="store_true", help="Capture TRACE logs per case and count search trials")
    parser.add_argument("--keep-going", action="store_true", help="Record failed bounds and finish other configurations")
    args = parser.parse_args()
    if args.repeats < 1 or not np.isfinite(args.scale) or args.scale <= 0:
        parser.error("repeats and scale must be positive")
    if any(not np.isfinite(target) or target <= 0 for target in args.targets):
        parser.error("targets must be finite and positive")
    scales = args.scales if args.scales is not None else [args.scale]
    if any(not np.isfinite(scale) or scale <= 0 for scale in scales):
        parser.error("scales must be finite and positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.trace:
        os.environ["EBCC_LOG_LEVEL"] = "0"
    rng = random.Random(0)
    all_passed = True
    with xr.open_dataset(args.input) as dataset:
        for variable in args.variables:
            array = dataset[variable]
            if args.pressure_levels is not None:
                if "pressure_level" not in array.dims:
                    parser.error(f"{variable} has no pressure_level dimension")
                array = array.sel(pressure_level=args.pressure_levels)
            if array.ndim < 2:
                parser.error(f"{variable} has fewer than two dimensions")
            if args.time_indices is not None:
                if array.ndim < 3:
                    parser.error("time selection requires a non-spatial dimension")
                array = array.isel({array.dims[0]: args.time_indices})
            data = array.values.astype(np.float32)
            if not np.isfinite(data).all():
                parser.error(f"{variable} contains nonfinite values")
            for repeat in range(args.repeats):
                for target in args.targets:
                    configs = [(mode, scale) for mode in args.modes
                               for scale in (scales if mode == "psnr" else [None])]
                    if args.randomize:
                        rng.shuffle(configs)
                    elif repeat % 2:
                        configs.reverse()
                    for mode, scale in configs:
                        passed = run_case(data, variable, target, mode, scale, repeat, args)
                        all_passed = passed and all_passed
    if not all_passed:
        raise SystemExit("One or more configurations violated the error bound; see results.jsonl")


if __name__ == "__main__":
    main()
