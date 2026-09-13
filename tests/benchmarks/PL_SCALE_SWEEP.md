# ERA5 pressure-level PSNR scale sweep

## Recommendation

For this pressure-level sample, use **6 as a single initial-guess scale** (`RMSE ≈ max_error / 6`). Scales 6 and 8 have almost identical total base-trial counts across screening and validation (1181 and 1184), but scale 6 used less compression CPU time in both phases. This differs from the surface-data recommendation of 8 for the feedback-loop objective. The codec default remains 8; this experiment does not change it.

For a bound-dependent choice, use **6 at 0.01% and 0.1%, and 8 at 1%**. At 0.1%, base counts are essentially tied: screening slightly favors 8 (186 versus 187), validation slightly favors 6 (190 versus 191); scale 6 is recommended there for its lower observed validation time and residual work.

For range-relative bound `r`, the seed is `20 * log10(scale / (0.99 * r))`, including the existing 1% optimization slack. Set `EBCC_PSNR_ERROR_SCALE` with `EBCC_BASE_SEARCH=psnr`. The scale controls the starting guess; decoded error checks enforce acceptance.

## Method and scope

- Input: `data/era5_pl_sample.nc`; experiment date: 2026-09-13; source commit: `399ade6`, with only pressure-level selection added to the benchmark driver.
- All six variables: geopotential (`z`), temperature (`t`), specific humidity (`q`), zonal wind (`u`), meridional wind (`v`), and pressure vertical velocity (`w`).
- Screening: 2024-03-01 at 1000, 500, and 100 hPa; scales 4, 6, 8, 16, 32; bounds 0.0001, 0.001, 0.01. This gives 90 configurations.
- Validation: 2024-03-02 at 850, 250, and 10 hPa; CR and PSNR scales 6 and 8; the same three bounds. This gives 54 configurations.
- Each configuration contains three full-resolution 721 × 1440 float32 chunks, one per pressure level, totaling 12,458,880 input bytes. No spatial cropping or downsampling. All selected values were finite.
- Total: **144 configurations / 432 spatial chunks**. This samples six of the file's 37 pressure levels, with different levels on the two days; it is not an exhaustive all-level study or a time-only holdout.
- Cases ran sequentially, once per configuration, with deterministic shuffled mode/scale order, TRACE capture, `base_cr=30`, and unchanged codec settings. Wall timing includes HDF5 creation, compression, writing and close; CPU timing measures process CPU over that interval. Input loading and decompression are excluded. MB/s is decimal.
- Every chunk was checked against exactly `relative_target * chunk_range`, without an added error allowance. **All 144 configurations passed**, with no PSNR lossless fallback. Worst observed error was 99.6192% of the requested bound.
- Timings are single executions, not confidence intervals or statistically established rankings. TRACE overhead is included in every case.

## Screening

Each row covers all six variables and three bounds (54 chunks).

| Scale | Base trials at 0.01% | At 0.1% | At 1% | Total base trials | Residual iterations | Write CPU (s) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 221 | 253 | 276 | 750 | 605 | 146.025 |
| 6 | 154 | 187 | 243 | 584 | 658 | 116.248 |
| 8 | 200 | 186 | 197 | 583 | 559 | 121.819 |
| 16 | 231 | 258 | 235 | 724 | 617 | 148.786 |
| 32 | 263 | 236 | 261 | 760 | 594 | 154.674 |

## Validation

Each row covers identical validation data and bounds: 224,259,840 input bytes across 54 chunks.

| Search | Scale | Base trials | Inferred base encode calls | Residual iterations | Write wall (s) | Write CPU (s) | Write MB/s | Compression ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CR | — | 821 | 1037 | 325 | 159.049 | 155.045 | 1.410 | 15.725× |
| PSNR | 6 | 597 | 705 | 573 | 111.832 | 109.022 | 2.005 | 16.192× |
| PSNR | 8 | 601 | 709 | 627 | 122.023 | 118.396 | 1.838 | 16.222× |
| PSNR per-bound | 6 / 6 / 8 | 556 | 664 | 587 | 106.268 | 103.601 | 2.110 | 16.200× |

Scale 6 reduced base trials by 27.3% versus CR, with 1.42× observed write throughput and 29.7% less CPU time. It was 1.09× faster than scale 8. PSNR used more residual-search iterations than CR, so the speedup does not imply both feedback loops improved. The per-bound row combines already measured cases; it is not an additional adaptive-policy run.

### Per-bound validation

| Error bound | Selected scale | CR base trials | Selected base trials | CR write (s) | Selected write (s) | Selected MB/s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.01% | 6 | 252 | 157 | 61.846 | 39.614 | 1.887 |
| 0.1% | 6 | 230 | 190 | 46.972 | 34.181 | 2.187 |
| 1% | 8 | 339 | 209 | 50.232 | 32.473 | 2.302 |

### Per-variable comparison

Base counts below combine both phases (six spatial slices × three bounds per scale). CPU times are validation only.

| Variable | Base trials, scale 6 | Base trials, scale 8 | Validation CPU, scale 6 (s) | Validation CPU, scale 8 (s) | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| z | 233 | 209 | 13.874 | 13.361 | 8 favors both trials and time |
| t | 200 | 191 | 18.967 | 19.099 | 8 favors trials; timing nearly tied |
| q | 192 | 197 | 20.880 | 23.783 | 6 favors both combined trials and validation time |
| u | 191 | 199 | 17.875 | 18.791 | 6 favors both trials and time |
| v | 199 | 192 | 17.613 | 20.294 | 8 favors combined trials; 6 was faster in validation |
| w | 166 | 196 | 19.814 | 23.069 | 6 favors both trials and time |

## Reproduce

The driver assumes `ebcc` is importable and uses its shipped `EBCC_FILTER_DIR`.

```bash
.venv/bin/python tests/benchmarks/compare_base_search.py \
  --input data/era5_pl_sample.nc --variables z t q u v w \
  --time-indices 0 --pressure-levels 1000 500 100 \
  --targets 0.0001 0.001 0.01 --modes psnr --scales 4 6 8 16 32 \
  --trace --randomize --keep-going --output-dir /private/tmp/ebcc-pl-screen-repro

.venv/bin/python tests/benchmarks/compare_base_search.py \
  --input data/era5_pl_sample.nc --variables z t q u v w \
  --time-indices 1 --pressure-levels 850 250 10 \
  --targets 0.0001 0.001 0.01 --modes cr psnr --scales 6 8 \
  --trace --randomize --keep-going --output-dir /private/tmp/ebcc-pl-validation-repro
```

Use new output folders. [Raw results](results/pl_scale_sweep.jsonl) include phase, selected time and pressure levels, timings, errors, sizes, and trace counters. Per-case traces and compressed files from this run are in `/private/tmp/ebcc-pl-scale-screen` and `/private/tmp/ebcc-pl-scale-validation`. Counter definitions follow [the base-search experiment](BASE_SEARCH.md); residual counts are truncation trials, not internal SPIHT bit-plane steps.

SHA-256 fingerprints:

```text
data/era5_pl_sample.nc
2eda14bea9d671a2de97dad00878b2f4e0cb671aece8746a28d64244385e5ef4

ebcc/libh5z_ebcc.cpython-311-darwin.so
a95f8fee3be472574ac2f7fa7176189cfe98c2682eb12dd8edaf005a6dc4a61a
```
