# ERA5 surface-variable PSNR scale sweep

## Recommendation

Use **8 as a general initial-guess scale** (`RMSE ≈ max_error / 8`). It minimizes
combined base-search trials across the two tested days among the finalists,
with substantially less residual work than scale 6. It replaces the initial
geopotential-only choice of 16 as the opt-in PSNR search's default. CR remains
the default search mode.

If selecting a scale per error bound, **6 at 0.01%**, and **8 at 0.1% and 1%**
performed better in both screening and validation. These are empirical seeds,
not error guarantees; the decoded error checks remain necessary.

| Relative error bound | Fraction passed to EBCC | Suggested scale | Initial PSNR with 1% optimization slack |
| --- | ---: | ---: | ---: |
| 0.01% | 0.0001 | 6 | 95.65 dB |
| 0.1% | 0.001 | 8 | 78.15 dB |
| 1% | 0.01 | 8 | 58.15 dB |

For a range-relative bound `r`, the seed is
`20 * log10(scale / (0.99 * r))`. The scale can be overridden through
`EBCC_PSNR_ERROR_SCALE`. The C implementation does not automatically switch
scale by bound: its default is the single compromise value 8.

## Method

Experiment date: 2026-09-13. Input: `data/era5_sfc_sample.nc`.

- Six variables: `u10`, `v10`, `t2m`, `msl`, `tisr`, and `tp` (10 m wind
  components, 2 m temperature, sea-level pressure, incident solar radiation,
  and precipitation).
- Both time slices: 2024-03-01 and 2024-03-02. All values were finite.
- Each trial used one full-resolution `721 × 1440` float32 spatial slice:
  4,152,960 input bytes, one HDF5 chunk. No spatial cropping/downsampling.
- `base_cr=30`, relative error mode, original residual algorithm and pure-base
  comparison, 0.25 dB PSNR search tolerance, and the existing 1% optimization
  slack. The filter binary was identical throughout the sweep.
- First-day screening: scales 4, 8, 16, 32, 64 for all six variables and three
  bounds (90 configurations). Refinement added 6 and 12 (36 configurations).
- Second-day validation: CR plus PSNR scales 6, 8, 16 across all variables and
  bounds (72 configurations). Thus **198 configurations** were executed.
- Mode/scale order was shuffled deterministically within each variable/bound
  group. Cases ran sequentially, with native TRACE capture and separate
  write/compression and read/decompression timers. Validation also recorded
  process CPU time.
- Every result checked maximum error per chunk against exactly
  `relative_target * chunk_range`, without the old 50% allowance.

**All 198 configurations passed.** No PSNR lossless fallback occurred.

After changing the default to 8 and rebuilding the plugin, **22 integration
and benchmark tests passed** with PSNR enabled, including NetCDF/HDF5, the
C API, offset-field/sparse-outlier regression cases, and both original
geopotential performance-benchmark cases (timing repetitions disabled).

## First-day scale selection

Each entry sums base-search trials over six variables. All candidates have
the same 18 variable/bound combinations. Residual iterations are the total
truncation-search iterations over those combinations.

| Scale | Base trials at 0.01% | At 0.1% | At 1% | Total base trials | Residual iterations |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 78 | 83 | 92 | 253 | 139 |
| 6 | 55 | 67 | 82 | 204 | 222 |
| 8 | 73 | 67 | 65 | 205 | 144 |
| 12 | 83 | 70 | 70 | 223 | 140 |
| 16 | 81 | 67 | 76 | 224 | 141 |
| 32 | 88 | 70 | 85 | 243 | 140 |
| 64 | 109 | 77 | 90 | 276 | 125 |

Scales 6 and 8 were nearly tied on base trials, but 6 needed 54% more residual
iterations. Scale 16 was retained as a control because it was the previous
default. First-day wall timings shifted between coarse and refinement phases;
they are retained in the raw data but not used to claim a speed advantage
between those phases.

## Second-day validation

Each row sums the same 18 variable/bound combinations: 74,753,280 input bytes.
Write timing includes HDF5 creation, compression, file writing and close;
it excludes input loading and read/decompression. CPU timing includes CPU
consumed by the process in that write interval. MB/s uses 1,000,000 bytes.

| Search | Scale | Base trials | Inferred base encode calls | Residual iterations | Write wall (s) | Write CPU (s) | Write MB/s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CR | — | 284 | 354 | 109 | 56.864 | 55.570 | 1.315 |
| PSNR | 6 | 221 | 257 | 180 | 43.597 | 42.542 | 1.715 |
| PSNR | 8 | 212 | 248 | 163 | 45.337 | 44.043 | 1.649 |
| PSNR | 16 | 257 | 293 | 197 | 57.503 | 56.049 | 1.300 |

Scale 8 reduced base trials by **25.4% versus CR** and **17.5% versus scale 16**.
Its observed write speedups were **1.25× versus CR** and **1.27× versus scale 16**.
Across both days, the finalist base-trial totals were 425 (scale 6),
**417 (scale 8)**, and 481 (scale 16). Corresponding residual totals were
402, **307**, and 338.

Scale 6 was about **4% faster in wall time** than 8 in validation, despite more
base and residual trials. Fewer iterations are not identical to less CPU work:
candidate sizes, PSNR levels, residual paths, and system conditions matter.
Choose 8 for the general feedback-loop reduction objective; 6 remains a useful
throughput alternative, especially at the tightest bound. These are single
executions per configuration, not statistically established runtime rankings.

Aggregate compression ratios (sum of input bytes / sum of output file bytes)
were 12.781× for CR, 13.312× for PSNR/6, **13.263× for PSNR/8**, and 13.284× for
PSNR/16. PSNR/8's worst observed error was 99.0361% of the requested bound.

### Per-bound choices

The following sums each selected bound/scale over all six second-day
variables. They are assembled from the validation cases, not an extra run
with an automatic scale-selection policy.

| Error bound | Scale | CR base trials → PSNR | CR write seconds → PSNR |
| --- | ---: | ---: | ---: |
| 0.01% | 6 | 106 → 63 | 23.462 → 15.787 |
| 0.1% | 8 | 96 → 65 | 17.883 → 13.400 |
| 1% | 8 | 82 → 66 | 15.519 → 10.746 |

Radiation behaves differently: in the first-day coarse sweep it favored scale
32, and among validation finalists it favored 16 for base-trial count. This
is another reason not to treat the chosen constant as universal.

## Reproduce and inspect

The driver assumes `ebcc` is importable and a plugin built from the experimental
PSNR source is installed in `ebcc.EBCC_FILTER_DIR`.

```bash
# First-day coarse screen
.venv/bin/python tests/benchmarks/compare_base_search.py \
  --input data/era5_sfc_sample.nc \
  --variables u10 v10 t2m msl tisr tp --time-indices 0 \
  --targets 0.0001 0.001 0.01 --modes psnr --scales 4 8 16 32 64 \
  --trace --randomize --keep-going --output-dir /private/tmp/ebcc-sfc-screen

# Repeat on day 0 with --scales 6 12 for refinement, using a new output folder.

# Second-day validation against CR and the previous scale
.venv/bin/python tests/benchmarks/compare_base_search.py \
  --input data/era5_sfc_sample.nc \
  --variables u10 v10 t2m msl tisr tp --time-indices 1 \
  --targets 0.0001 0.001 0.01 --modes cr psnr --scales 6 8 16 \
  --trace --randomize --keep-going --output-dir /private/tmp/ebcc-sfc-validation
```

`--trace` writes a separate log per case and adds iteration counters to
`results.jsonl`. `--keep-going` records all error-bound failures and exits
unsuccessfully after the sweep if any occurred. `--repeats` is available for
follow-up timing studies. Use a new output directory for each experiment.

[Raw results for all 198 configurations](results/sfc_scale_sweep.jsonl) include
the phase, variable, selected time slice, target, scale, timings, compression
ratio, error, and trace counters. Only validation records have CPU timings;
that instrumentation was added after screening. Counter definitions match
the [base-search experiment](BASE_SEARCH.md): CR counts its loop trials,
PSNR counts quality trials, and encode-call totals include inferred initial
and final encodes. Residual counts exclude SPIHT's internal bit-plane steps.

SHA-256 fingerprints during the sweep (before changing only the default scale
from 16 to 8; all sweep runs explicitly specified their scale):

```text
data/era5_sfc_sample.nc
751e24fcc11f7ab42fdbaf1c2fa4e74e9d58a392612af699f2706f5f971b454f

ebcc/libh5z_ebcc.cpython-311-darwin.so
a65f9e0699adcb8663bc3b162d14eb468dfc0e66eac6bd415b00ae6e822e4720
```

This covers two days from one ERA5 sample, not all weather regimes or chunk
sizes. The known permissive pure-base fallback for unattainable absolute
bounds remains outside this scale-tuning experiment.
