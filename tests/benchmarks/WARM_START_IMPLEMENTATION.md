# Integrated JPEG2000 warm start

Warm start is enabled by default for single-tile EBCC error-bound compression.
Set `EBCC_WARM_START=0` to compare with the original cold path. It applies to
both CR and PSNR searches, including their final pure-base comparison. Search
targets, iteration logic, error checks, compressed format and decoding remain
unchanged. Constant fields, single-encode `NONE` mode, multi-tile encodes and
reversible/lossless fallback use their existing paths.

## Ownership and build

[ebcc_j2k.c](../../src/ebcc_j2k.c) implements an opaque prepared-input session.
Creation copies one unsigned 16-bit grayscale tile. Encoding accepts only a
new PSNR or rate, so the input, dimensions, precision, component count, wavelet,
resolution count, code-block settings and layer count cannot change within a
session. `ebcc_encode` owns this session and passes it explicitly through both
base searches. Destruction occurs before freeing the chunk's quantized input;
there is no global or thread-local cache.

The first lossy encode performs DC shift, DWT and Tier-1 code-block encoding.
It retains code-block bitstreams, pass lengths/distortions and tile metadata.
Transformed image/component pixel buffers are freed immediately after this
preparation. Subsequent requests create target-specific headers, rerun rate
allocation and Tier-2 packetization, and emit independent J2K codestreams.
Existing rate allocation and packetization reset layer pass counts and tag
trees; the bridge explicitly resets the packet counter and restores borrowed
tile ownership on both success and failure.

A creation failure selects cold encoding. An encode failure releases the
prepared codec, disables reuse for the rest of that call, discards any partial
output and retries through the cold path. A session failure cannot leak cached
state into another chunk. Lossless requests bypass the lossy cache because
their transform differs.

[OpenJPEGWarmStart.cmake](../../src/cmake/OpenJPEGWarmStart.cmake) copies the
bundled core source/header files into the build tree and extends only that
copy. It adds an instance-local reuse flag to the private tile coder, guards
the input-dependent preparation stages, and includes
[ebcc_openjpeg_bridge.inc](../../src/ebcc_openjpeg_bridge.inc) in the J2K
translation unit. All core translation units use the same extended header.
The upstream submodule stays pristine. Hook checks fail configuration if an
upstream update removes or duplicates the expected insertion points.

The bundled OpenJPEG and bridge symbols are hidden in the EBCC shared library
to prevent binding the extended private state to a different OpenJPEG loaded
by the host. CMake static consumers link the bundled, extended OpenJPEG target;
a system OpenJPEG is not a substitute for that private dependency.

## Validation

- 48 NetCDF/HDF5, C API, Zarr and original performance-benchmark tests passed
  with warm start enabled in each search mode: 96 test executions. Benchmark
  timing repetitions were disabled in these integration runs.
- New Python cases compare exact cold/warm bytes at all three relative bounds,
  cover multi-tile bypass and constant fields, exercise concurrent independent
  chunks, and assert exactly one preparation per HDF5-sized spatial chunk.
- The native session regression passed under AddressSanitizer. It tests an
  owned input snapshot, non-monotone PSNR/CR requests, repeated targets,
  first-encode failure, post-preparation write failure, disabled-session reuse,
  and destruction after failure. It does not rely on assertions enabled only
  in debug builds.
- Symbol inspection of the final plugin found no exported `opj_*` or
  `ebcc_opj_*` entry points.

## Measurements

The paired C API comparison uses all six pressure-level variables (`z`, `t`,
`q`, `u`, `v`, `w`) at 500 hPa on 2024-03-01, full 721 × 1440 spatial chunks,
relative bounds 0.0001, 0.001, 0.01, both search modes, and two repetitions.
PSNR scale is 6; `base_cr=30`. Cold/warm order alternates by repetition.
Compression timing includes EBCC's internal decoding and residual work,
but excludes final round-trip decoding and HDF5/disk overhead. Logging is WARN.
These are small-sample performance measurements, not confidence intervals.

**All 72 paired comparisons produced identical compressed bytes and decoded arrays, and both paths met the error bound in all 144 encodes.** Their output SHA-256 values also match all 72 corresponding pre-integration reference records.

| Search | Cold mean/chunk (s) | Warm mean/chunk (s) | Speedup | Cold MB/s | Warm MB/s |
| --- | ---: | ---: | ---: | ---: | ---: |
| CR | 3.013 | 1.000 | 3.01× | 1.378 | 4.153 |
| PSNR | 2.134 | 0.788 | 2.71× | 1.946 | 5.273 |

[Raw paired measurements](results/warm_start_integrated.jsonl) include wall/CPU timings, output sizes, errors and hashes. This measures the final build with symbol visibility isolation; the PSNR scale default remains 8, while this comparison explicitly uses scale 6.

### Trace counts per chunk

Each row is one 721 × 1440 chunk at 500 hPa, time index 0, relative bound 0.1%, scale 6. Base iterations include both searches; emit counts also include initial and final encodes.

| Variable | Search | Base iterations | Base emits | Full preparations | Reuses | Residual iterations |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| z | CR | 18 | 22 | 1 | 21 | 14 |
| z | PSNR | 14 | 16 | 1 | 15 | 0 |
| q | CR | 14 | 18 | 1 | 17 | 0 |
| q | PSNR | 12 | 14 | 1 | 13 | 17 |

[Trace-derived records](results/warm_start_integrated_trace.jsonl) retain timings and error checks. Native logs are in `/private/tmp/ebcc-warm-final-trace`.

### Sampled process memory

Fresh processes each compressed the full-resolution humidity slice at 500 hPa three times with PSNR, scale 6, and relative bound 0.01%. RSS was sampled every 5 ms after loading input and the library.

| Setting | Baseline RSS (MiB) | Peak RSS (MiB) | Increase over baseline (MiB) |
| --- | ---: | ---: | ---: |
| Cold | 204.22 | 421.92 | 217.70 |
| Warm | 205.42 | 422.38 | 216.95 |

There was no material peak-RSS increase in this check. RSS includes the process, allocator retention and residual compressor; it is not a direct measurement of cache bytes and can miss peaks between samples. Parallel peak memory and other inputs were not measured. [Raw memory samples](results/warm_start_integrated_memory.jsonl).

Final plugin SHA-256 (macOS arm64):

```text
5294f8624fb83a38cd6a9bacde03717c76054776c9eb0c774e10459a3ccc07f4
```

## Reproduce

Build/install EBCC normally with its bundled OpenJPEG. The Python commands
assume `ebcc` is importable and use the workdir environment.

```bash
# Replace the library filename for the current platform.
.venv/bin/python tests/benchmarks/warm_start/run.py --skip-micro \
  --ebcc-library ebcc/libh5z_ebcc.cpython-311-darwin.so \
  --output /private/tmp/ebcc-warm-integrated-repro.jsonl

# Native lifecycle/failure checks, with AddressSanitizer where supported.
cmake -S src -B /private/tmp/ebcc-warm-asan-repro \
  -DEBCC_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS='-fsanitize=address -fno-omit-frame-pointer'
cmake --build /private/tmp/ebcc-warm-asan-repro --target test_j2k_session -j 4
ctest --test-dir /private/tmp/ebcc-warm-asan-repro --output-on-failure

# Run each memory setting in a fresh process.
.venv/bin/python tests/benchmarks/warm_start/memory.py \
  --library ebcc/libh5z_ebcc.cpython-311-darwin.so --warm 0
.venv/bin/python tests/benchmarks/warm_start/memory.py \
  --library ebcc/libh5z_ebcc.cpython-311-darwin.so --warm 1
```

TRACE logs distinguish `Base warm start: prepare` from `Base warm start: reuse`.
`compare_base_search.py --trace` records both counts alongside base and residual
iterations. Rate allocation and feedback decoding still occur on every trial;
the count of expensive transform/Tier-1 preparations falls to one per eligible
chunk.

The [earlier prototype report](WARM_START.md) remains a historical comparison.
Its build script pins the pre-integration codec source so it stays reproducible.
