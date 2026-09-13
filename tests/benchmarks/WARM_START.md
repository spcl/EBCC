# Reusing OpenJPEG work across EBCC search trials

## Finding

**Warm start is feasible and substantially reduces measured compression time.** Retaining JPEG2000 Tier-1 code-blocks makes subsequent quality/rate requests about **36× faster** in the isolated encoder probe. Including the first full preparation, a sequence of ten requests is **8.07× faster**. In complete EBCC compression, the measured speedup is **3.00× with CR search** and **2.69× with PSNR search**. All compared outputs were byte-identical.

This report records the initial isolated prototype. Warm start has since been [integrated into EBCC](WARM_START_IMPLEMENTATION.md). The original experiment left the shipped library, `src/ebcc_codec.c`, OpenJPEG submodule, and defaults unchanged. Its historical prototype lives under [warm_start/](warm_start/) and builds temporary source copies pinned to the original codec commit.

## What can be reused

Current EBCC creates a new encoder for every trial. OpenJPEG then performs DC level shifting, the wavelet transform, Tier-1 code-block encoding, rate allocation, and Tier-2 packet writing. For fixed input and coding parameters, changing only PSNR or CR requires the latter two stages.

The retained state includes the code-block bitstreams, coding-pass byte lengths and distortion reductions, tile distortion, and supporting tile structures. Each request supplies a new target, reruns rate allocation, and writes a normal standalone single-layer codestream. There is no approximation or reuse of a previously decoded candidate.

In the prototype, the first lossy request does the full transform and Tier-1 encoding. All later lossy requests in the same `ebcc_encode` call reuse that state, including the pure-base comparison and switches between CR and PSNR. Lossless requests retain the existing cold path because the reversible wavelet transform differs. The cache is released before the quantized input is freed.

**Search iteration counts are unchanged.** One expensive lossy transform/Tier-1 preparation replaces repeated preparations for the chunk; rate allocation, packet writing, decoding, and error checks still run for every candidate. Residual compression also remains unchanged. This explains why full EBCC speedup is much smaller than the repeated-encode speedup.

## Prototype implementation

- `prepare.py` copies the pinned OpenJPEG core sources into a supplied temporary directory. It never edits the submodule.
- In the temporary `tcd.c`, a thread-local pointer identifies the cached tile and skips only DC shift, MCT, DWT, and Tier-1 preparation for that tile. Rate allocation and Tier-2 still execute.
- `j2k_warm.inc` appends two experimental helpers to temporary `j2k.c`: finish the first stream without destroying its tile state, and emit a new stream using that state.
- A fresh target encoder still generates headers and converts CR into a tile byte budget. The cached tile is temporarily attached to that encoder for allocation and packetization, then detached before cleanup. Layer pass counts and packet tag trees are reset by the existing allocation/Tier-2 code; the tile packet counter is explicitly reset.
- The temporary EBCC build routes eligible calls through this cache when `EBCC_WARM_START=1`. Cold and warm measurements use the same temporary library. The shipped library did not recognize this variable at the time of the experiment; the integrated implementation now supports it.
- Eligibility is deliberately narrow: the same quantized input, one tile, one unsigned 16-bit grayscale component, one quality layer, default JPEG2000 coding parameters, and the irreversible transform. Multi-tile and reversible/lossless calls use the original path.

OpenJPEG's current public interface supports targets at encoder setup and predefined quality layers, but does not expose a prepared-Tier-1 reuse operation. `opj_end_compress` destroys the tile coder. Production integration therefore needs a maintained OpenJPEG extension or patch with an explicit prepared-input lifetime. Merely calling encode again on a completed public codec is insufficient. See the pinned [tile encoder](https://github.com/uclouvain/openjpeg/blob/210a8a5690d0da66f02d49420d7176a21ef409dc/src/lib/openjp2/tcd.c), [codec lifecycle](https://github.com/uclouvain/openjpeg/blob/210a8a5690d0da66f02d49420d7176a21ef409dc/src/lib/openjp2/j2k.c), and [public encoder parameters](https://www.openjpeg.org/doxygen/structopj__cparameters.html).

## Measurement method

- Date: 2026-09-13. EBCC source commit: `399ade6`. OpenJPEG commit: `210a8a5690d0da66f02d49420d7176a21ef409dc`.
- Input: all six `z`, `t`, `q`, `u`, `v`, `w` variables from `data/era5_pl_sample.nc`, 2024-03-01, 500 hPa. Each input is a full 721 × 1440 float32 spatial slice, 4,152,960 bytes.
- JPEG2000 probe: each variable uses a ten-request sequence, twice: PSNR 96, 94, 95, 78, 58 dB; OpenJPEG CR 15, 50, 500; PSNR 96, 58 dB. This tests decreasing/increasing quality, repeated targets, and CR/PSNR switching. OpenJPEG CR is half EBCC's float32-relative CR.
- Complete EBCC probe: each variable × three range-relative targets (0.0001, 0.001, 0.01) × CR/PSNR search × two repetitions. `base_cr=30`, PSNR scale 6, and the existing 1% optimization slack.
- Cold/warm execution order alternates by repetition. Cases execute sequentially. Native compression wall time and process CPU time are measured; Python byte copying/comparison and final round-trip decoding are outside that timer. EBCC's internal feedback-loop decoding is inside it. Logging is WARN for both paths.
- These are **C API compression times**, without HDF5/disk overhead, and should not be directly compared to the earlier HDF5 throughput numbers. Two repeats provide an initial performance estimate, not confidence intervals.

## JPEG2000 encode timings

| Requests | Paired measurements | Cold mean (ms) | Warm mean (ms) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| All requests, including first preparation | 120 | 120.461 | 14.923 | 8.07× |
| Subsequent requests only | 108 | 120.157 | 3.317 | 36.22× |
| Subsequent PSNR requests | 72 | 120.051 | 3.516 | 34.14× |
| Subsequent CR requests | 36 | 120.370 | 2.919 | 41.24× |

## Complete EBCC compression

Each search row covers 36 paired measurements: six variables × three bounds × two repetitions.

| Search | Paired measurements | Cold mean/chunk (ms) | Warm mean/chunk (ms) | Speedup | Cold MB/s | Warm MB/s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CR | 36 | 3051.341 | 1016.997 | 3.00× | 1.361 | 4.084 |
| PSNR | 36 | 2166.267 | 805.594 | 2.69× | 1.917 | 5.155 |

### By error bound

| Search / bound | Paired measurements | Cold mean/chunk (ms) | Warm mean/chunk (ms) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| CR / 0.01% | 12 | 3589.104 | 1540.516 | 2.33× |
| CR / 0.1% | 12 | 2647.393 | 902.695 | 2.93× |
| CR / 1% | 12 | 2917.527 | 607.779 | 4.80× |
| PSNR / 0.01% | 12 | 2123.255 | 1116.351 | 1.90× |
| PSNR / 0.1% | 12 | 2148.461 | 774.829 | 2.77× |
| PSNR / 1% | 12 | 2227.086 | 525.601 | 4.24× |

### By variable, both search modes and all bounds

| Variable | Paired measurements | Cold mean/chunk (ms) | Warm mean/chunk (ms) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| z | 12 | 1940.175 | 555.939 | 3.49× |
| t | 12 | 2290.560 | 761.603 | 3.01× |
| q | 12 | 3403.330 | 1262.605 | 2.70× |
| u | 12 | 2638.789 | 943.636 | 2.80× |
| v | 12 | 2499.687 | 843.061 | 2.97× |
| w | 12 | 2880.285 | 1100.927 | 2.62× |

## Correctness checks

- **120/120 JPEG2000 paired measurements had identical codestream bytes.** This includes 108 subsequent requests that reused previously prepared state.
- **72/72 complete EBCC paired measurements had identical compressed bytes and decoded float arrays.** Both paths passed the exact per-chunk relative error bound in all 144 full-codec executions.
- The existing **10 C API regression tests passed with warm start enabled in each search mode** (20 test executions). These include chunking/padding, relative bounds, plain EBCC payloads, and lossless fallback regressions. Tests loaded the temporary `libebcc_probe` through the fixture's `EBCC_FILTER_PATH`, leaving the installed library untouched.
- [Raw paired results](results/warm_start.jsonl) contain timings, encoded sizes, per-case output SHA-256 hashes, and measured maximum errors. Byte equality also preserves compression ratio and residual payload decisions.

## Production integration requirements

Expose a prepared encoder object that owns the immutable input-dependent state and releases it at chunk completion. Provide a target-setting operation that resets mutable layer/packet state and produces a standalone codestream. Validate that dimensions, precision, signedness, transform, quantization, resolution levels, code-block settings, ROI, and component transform settings match before reuse. Give reversible/lossless coding a separate prepared state or keep its cold fallback.

The current probe retains the first encoder's tile/image allocations until chunk completion, so retained memory increases; peak memory and parallel chunk scaling were not measured. A production implementation can discard transform/input buffers once no longer needed while retaining code-block streams and pass statistics. It should use explicit ownership instead of the probe's private-structure swapping and add concurrency, broader input/settings, allocation-failure, and memory-sanitizer tests.

A practical next step is a small maintained OpenJPEG extension exposing prepare / emit-at-target / destroy, used by one EBCC chunk across both base searches. The measured gain is large enough to justify that work. Predefined multi-quality layers are another public-API option, but require choosing targets in advance; that alternative was not benchmarked here.

## Reproduce

Use the workdir Python environment. The build needs headers and the zstd static library from an existing EBCC build; adjust the platform-specific build directory below.

```bash
.venv/bin/python tests/benchmarks/warm_start/prepare.py \
  --work-dir /private/tmp/ebcc-warm-start-repro \
  --config-dir build/temp.macosx-11.0-arm64-cpython-311/ebcc.libh5z_ebcc/openjpeg/src/lib/openjp2 \
  --zstd-library build/temp.macosx-11.0-arm64-cpython-311/ebcc.libh5z_ebcc/openjpeg/bin/libzstd.a

.venv/bin/python tests/benchmarks/warm_start/run.py \
  --j2k-library /private/tmp/ebcc-warm-start-repro/build/libwarm_probe.dylib \
  --ebcc-library /private/tmp/ebcc-warm-start-repro/build/libebcc_probe.dylib \
  --output /private/tmp/ebcc-warm-start-repro/results.jsonl
```

The runner assumes `ebcc` is importable; it does not insert system paths. Use a fresh results filename. Linux builds produce `.so` library names; this probe was validated on macOS arm64.

Experiment binary SHA-256:

```text
libebcc_probe.dylib
68510a78618484aeb17b94deb38fdb7672c9cb80cbcb3f9e727df6e82505d9c1

libwarm_probe.dylib
038aa0900818c4937b32a78a0f8a2b9410aebb951a859803919ba1ea4916717e
```
