# PSNR base-search experiment

Measured on 2026-09-13 against the CR search at source baseline
`91b4b52f76ef93c1275b56fe8af6ea0fec123291`, using one rebuilt library for both
modes. This is an opt-in experiment: CR remains the default.

## Search

Set `EBCC_BASE_SEARCH=psnr` to use OpenJPEG's fixed-quality allocation
(`cp_fixed_quality`, `tcp_distoratio`) instead of rate allocation for the
base-layer feedback searches. The HDF5 options and encoded format are unchanged.

The initial guess assumes `RMSE = absolute_error_target / scale`:

```text
initial_PSNR = 20 log10(chunk_range * scale / absolute_error_target)
```

This initial experiment used `EBCC_PSNR_ERROR_SCALE=16`. The subsequent
[surface-variable sweep](SFC_SCALE_SWEEP.md) selected 8 as the default.
The existing 1% optimization slack is
already included in `absolute_error_target`. For relative targets of 1% and
10%, scale 16 gives initial guesses of approximately 64.17 dB and 44.17 dB.

Each trial decodes the candidate and measures RMSE, maximum error, and the
fraction within the requested error bound. RMSE alone does not establish the
pointwise bound: acceptance still uses the existing error-quantile target.
The observed tail/RMSE ratio predicts the next PSNR step until feasible and
infeasible candidates bracket the solution; bisection then refines to 0.25 dB.
The search is bounded to 24 trials and PSNR values from 1 to 120 dB. It retains
the smallest feasible payload and its decoded data to avoid a final re-encode.
If none is feasible, it uses the existing lossless-JP2 fallback semantics.

The initial CR encode and later pure-base consistency encode are retained.
The residual search and pure-base candidate selection are unchanged. Thus
results include both the change of search variable and removal of final
feasible re-encodes; this is not an ablation isolating those two effects.

## Input and timing scope

- `data/geopotential_pl_small.nc`, variable `z`, shape `(2, 11, 721, 1440)`.
- Float32 input: 91,365,120 bytes, in 22 chunks of `(1, 1, 721, 1440)`.
- `base_cr=30`; chunk-range-relative targets 0.01 and 0.1.
- macOS, Python 3.11.14, Release CMake build; TRACE logging enabled for all
  comparisons. No parallel benchmark runs were launched.
- One dataset execution per configuration. These are exploratory single-run
  timings, not statistically established speedups; scale 16 was selected using
  this dataset and still needs validation on other variables/data.
- Write time includes HDF5 file/dataset creation, compression, writing and
  close. Input loading and float32 conversion are excluded. Read/decompression
  is timed separately. MB/s uses 1,000,000 bytes. These are not codec-only rates.

## Results

| Target | Search | Scale | Base trials | JPEG2000 encode calls* | Residual iterations | Write time (s) | Write MB/s | Compression ratio | Max chunk-relative error |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1% | CR | — | 594 | 678 | 97 | 92.5515 | 0.9872 | 346.759 | 0.00986228 |
| 1% | PSNR | 4 | 310 | 354 | 99 | 47.8992 | 1.9074 | 349.363 | 0.00990433 |
| 1% | PSNR | 16 | 223 | 267 | 98 | 36.2379 | 2.5213 | 347.334 | 0.00985489 |
| 10% | CR | — | 319 | 368 | 0 | 49.9173 | 1.8303 | 1668.647 | 0.09885735 |
| 10% | PSNR | 4 | 370 | 414 | 0 | 53.9034 | 1.6950 | 3205.119 | 0.09885735 |
| 10% | PSNR | 16 | 307 | 351 | 0 | 45.1157 | 2.0251 | 3177.254 | 0.09860717 |

All six configurations passed the exact per-chunk relative error assertion.
Read/decompression times were 0.12–0.19 seconds.

Validation after selecting scale 16: 22 tests passed with PSNR enabled
(NetCDF/HDF5, C API, and both original performance-benchmark cases with timing
repetitions disabled), including the new offset-field/sparse-outlier tests.
The 20 NetCDF/HDF5 and C API tests also passed with CR selected (the three
PSNR-specific tests explicitly override that selection).

At scale 16, base trials fell 62.5% at 1% error and 3.8% at 10% error.
Observed write-time speedups were 2.55x and 1.11x respectively. Residual
iterations were effectively unchanged. At 10%, PSNR also explored stronger
compression beyond the CR search's upper stopping threshold, so the outputs
do not have matching compression ratios.

*Counts come from TRACE messages and source control flow. CR trials are its
bracketing/bisection loop executions; PSNR trials include its first quality
candidate. Encode-call totals also include the initial and consistency encodes
(44 per target) and, for CR, final feasible re-encodes (40 at 1%, 5 at 10%).
No PSNR lossless fallback occurred in these dataset comparisons. Encode-call
counts are inferred, not direct function-entry profiler counts. Residual
iterations count truncation-search trials, not internal SPIHT bit-plane steps.

## Reproduce

Build the library using the project's CMake configuration and copy the resulting
plugin into the installed `ebcc` package before running. Both modes must load
the same rebuilt plugin. The scripts assume `ebcc` is already importable.

```bash
# Same-process, same-library comparison; retain traces and JSONL results.
EBCC_LOG_LEVEL=0 .venv/bin/python tests/benchmarks/compare_base_search.py \
  --scale 16 --output-dir /private/tmp/ebcc-base-search-compare \
  > /private/tmp/ebcc-base-search-results.log \
  2> /private/tmp/ebcc-base-search-trace.log

# Original benchmark and regression tests with PSNR enabled.
EBCC_BASE_SEARCH=psnr EBCC_LOG_LEVEL=3 .venv/bin/python -m pytest \
  tests/test_netcdf.py tests/test_c_api.py \
  tests/benchmarks/test_compression_benchmarks.py::TestCompressionBenchmarks::test_ebcc_compression_performance \
  --benchmark-disable -q
```

The comparison driver accepts `--repeats`, alternates mode order between
repetitions, and appends records to `results.jsonl` in the requested output
directory. Use a fresh output directory for each experiment.

Trace counting: within each `BEGIN`/`END` segment, count `PSNR trial:` for PSNR;
for CR count `cr_lo:` / `cr_hi:` lines containing `jp2_length:` and
`current_cr:` lines containing `cr_lo:`. Residual iterations are `trunc_lo:`
lines without `optimization_error_target:`. Initialization and summary lines
are excluded.

## Remaining limitations

The max-error/RMSE relationship depends on the field's error distribution;
scale 16 is an empirical seed, not an error guarantee. The decoded quantile
check and residual correction remain necessary. The existing behavior that
can accept an out-of-bound pure-base fallback for unattainable absolute
targets is not fixed by this experiment.

Before changing the default, repeat timings and test other variables, sparse
outliers, scales and chunk sizes. Reusing the second base search is a separate
optimization to evaluate after this comparison.
