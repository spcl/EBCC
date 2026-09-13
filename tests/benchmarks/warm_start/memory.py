"""Sample RSS during repeated compression; run each setting in a fresh process."""

import argparse
import ctypes as c
import gc
import json
import os
from pathlib import Path
import threading

import numpy as np
import psutil
import xarray as xr
from ebcc.zarr_filter import CodecConfigT
from run import load_codec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--warm", choices=["0", "1"], required=True)
    args = parser.parse_args()
    os.environ.update(EBCC_WARM_START=args.warm, EBCC_BASE_SEARCH="psnr",
                      EBCC_PSNR_ERROR_SCALE="6", EBCC_LOG_LEVEL="3")
    with xr.open_dataset("data/era5_pl_sample.nc") as dataset:
        data = np.ascontiguousarray(dataset.q.isel(valid_time=0).sel(pressure_level=500).values)
    lib = load_codec(args.library)
    config = CodecConfigT()
    config.dims[:] = (1, *data.shape)
    config.base_cr = 30
    config.residual_compression_type = 2
    config.error = 0.0001
    gc.collect()
    process = psutil.Process()
    baseline = process.memory_info().rss
    peak = [baseline]
    stop = threading.Event()

    def sample():
        while not stop.wait(0.005):
            peak[0] = max(peak[0], process.memory_info().rss)

    sampler = threading.Thread(target=sample)
    sampler.start()
    try:
        for _ in range(3):
            output = c.c_void_p()
            size = lib.ebcc_encode(data.ctypes.data, c.byref(config), c.byref(output))
            assert size and output.value
            lib.free_buffer(output)
    finally:
        stop.set()
        sampler.join()
    print(json.dumps(dict(warm=int(args.warm), baseline_rss=baseline,
                          peak_rss=peak[0], increase_rss=peak[0] - baseline,
                          end_rss=process.memory_info().rss,
                          variable="q", pressure_level=500, time_index=0,
                          target=0.0001, scale=6, mode="psnr", encodes=3,
                          sampling_seconds=0.005)))


if __name__ == "__main__":
    main()
