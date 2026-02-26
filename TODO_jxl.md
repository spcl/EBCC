# TODO: Add `libjxl` as an Additional Base Compressor

## Goal
- Add JPEG XL (`libjxl`) as a second base compressor backend alongside current JPEG2000 (`openjpeg`), without breaking existing users, data, or filter configurations.

## Non-goals
- Do not remove or regress current JPEG2000 path.
- Do not change EBCC residual compression semantics.
- Do not break legacy decode support.
- Windows support is out of scope for the initial implementation (CI workflow already has Windows commented out). Can be revisited later.

## Current State (for context)
- Base compression is hardcoded to JPEG2000 in `src/ebcc_codec.c` (`j2k_encode_internal`, `j2k_decode_internal`, `error_bound_j2k_compression`).
- HDF5 filter options (`cd_values` array) currently encode:
  - `[height, width, base_cr, residual_type, residual_param?]`
  - `base_cr` and `residual_param` are IEEE 754 float bit-patterns reinterpreted as `uint32` (via `memcpy`/`float_to_uint32`), not truncated integers.
  - `cd_nelmts` is 4 (when `residual_type == NONE`) or 5 (when `MAX_ERROR` or `RELATIVE_ERROR`).
- Serialized EBCC header is a fixed 48-byte struct (`ebcc_header_t`) enforced by `static_assert`:
  - `magic[4]` ("EBCC"), `version` (1), `flags` (1 byte), `reserved` (2 bytes), then payload metadata fields.
  - The 2-byte `reserved` field is currently unused and zero-filled.
- Legacy decode path (`ebcc_decode_legacy`) handles pre-header data without the magic/version prefix.

## High-Level Design
- Introduce a base compressor abstraction in C as a **function-pointer dispatch table**:
  ```c
  typedef struct {
      size_t (*encode)(uint16_t *data, int h, int w, float param, codec_data_buffer_t *out);
      size_t (*decode)(codec_data_buffer_t *stream, int h, int w, uint16_t **out);
      float  (*error_bound_search)(uint16_t *scaled_data, size_t *image_dims, size_t *tile_dims,
                                   float param, codec_data_buffer_t *buf, float **decoded,
                                   float minval, float maxval, float *data, size_t tot_size,
                                   float error_target, double base_quantile_target);
      const char *name;        // for logging
      uint8_t id;              // for serialization
      const char *param_name;  // "cr" or "distance"
      float param_min;
      float param_max;
  } base_codec_backend_t;
  ```
  - `encode` writes into a caller-owned `codec_data_buffer_t` (cleared/rewound by caller before each call); returns the number of bytes written.
  - `decode` reads from a `codec_data_buffer_t` whose `offset` is already rewound to the start of the compressed payload.
  - `param_name` / `param_min` / `param_max` expose the backend-specific quality knob so generic error-bound search code can validate and log it without hardcoding backend assumptions.
- Keep JPEG2000 implementation as the first backend (`id=0`).
- Add JPEG XL implementation as a second backend (`id=1`).
- Record selected base compressor in the serialized payload (header v2), while keeping v1/legacy decode support.

## Phase 0: Feasibility Spike
- [ ] Write a standalone C test program (outside EBCC) that:
  - Encodes a uint16 grayscale buffer with `libjxl` at various `distance` values.
  - Decodes it back and measures max/mean error in uint16 space.
  - Measures compressed size to understand the `distance` → compression ratio relationship.
- [ ] Verify `libjxl` API supports 16-bit single-channel (grayscale) encoding without issues.
  - `JxlPixelFormat` with `data_type = JXL_TYPE_UINT16`, `num_channels = 1`.
- [ ] Measure encode/decode speed relative to OpenJPEG at comparable compression ratios.
  - This matters because the error-bound binary search does multiple encode-decode cycles.
  - If JXL is slower per iteration, consider using lower `effort` values during the search phase.
- [ ] Document findings: quality-vs-distance curve, speed characteristics, any API quirks.

**Rationale:** De-risks the entire effort. If libjxl's 16-bit grayscale support is awkward or the distance-to-ratio mapping is unpredictable, we need to know before touching any EBCC code.

## Phase 1: Build and Dependency Strategy
- [ ] Choose dependency strategy (must decide before any code changes):
  - **Option A: vendored/submodule `libjxl`** — consistent builds, larger repo, more complex CMake.
  - **Option B: system `libjxl` via `find_package`/pkg-config** — simpler repo, but CI must install it.
- [ ] Make JXL support **optional at compile time**:
  - Add CMake option: `-DENABLE_JXL=ON/OFF` (default OFF initially, flip to ON once stable).
  - When disabled, `BASE_COMPRESSOR_JXL` returns a clear runtime error: "JXL support not compiled".
  - Conditionally compile JXL source files and link JXL libraries only when enabled.
- [ ] Update CMake (`src/CMakeLists.txt`):
  - Add includes/libs for JXL behind the `ENABLE_JXL` guard.
  - Link `h5z_ebcc` and `ebcc` targets with JXL when enabled.
  - Add clear configure-time message: "JXL support: ON/OFF".
- [ ] Update `setup.py`:
  - Respect `ENABLE_JXL` environment variable or CMake cache.
  - Ensure build works in wheel jobs with chosen dependency strategy.
  - Keep submodule init logic consistent if vendored.
- [ ] Update `.gitmodules` if vendoring is chosen.

## Phase 2: Header v2 Design and Serialization
- [ ] Repurpose the `reserved` field in `ebcc_header_t`:
  - Split `uint16_t reserved` into `uint8_t base_compressor` + `uint8_t reserved`.
  - `base_compressor = 0` → J2K, `base_compressor = 1` → JXL.
  - This maintains the 48-byte header size (static assert preserved).
  - **Key insight:** v1 headers have `reserved = 0x0000`, so the new `base_compressor` byte is already `0` (J2K) in all existing data. No special v1-vs-v2 detection needed beyond checking version.
- [ ] Bump `EBCC_HEADER_VERSION` to `2`.
- [ ] Decode behavior:
  - v2 payload: read `base_compressor` field, dispatch decoder accordingly.
  - v1 payload: `base_compressor` byte is `0` (from old `reserved`), so J2K is selected automatically. Alternatively, ignore the field for `version == 1` and hardcode J2K.
  - Legacy payload (no magic): keep `ebcc_decode_legacy` unchanged.
- [ ] Add validation: reject unknown `base_compressor` values with a clear error message.

## Phase 3: Codec Abstraction Refactor
- [ ] Define `base_codec_backend_t` dispatch table as described in High-Level Design.
- [ ] Create backend instances:
  ```c
  static const base_codec_backend_t j2k_backend = { j2k_encode, j2k_decode, error_bound_j2k_compression, "j2k", 0, "cr",        1.0f, 1000.0f };
  static const base_codec_backend_t jxl_backend = { jxl_encode, jxl_decode, error_bound_jxl_compression, "jxl", 1, "distance",  0.0f,   25.0f };
  ```
- [ ] Refactor `ebcc_encode` / `ebcc_decode` to select backend from config/header:
  - Replace direct `j2k_encode_internal()` / `j2k_decode_internal()` calls with `backend->encode()` / `backend->decode()`.
- [ ] Refactor error-bound search — keep separate per-backend implementations, exposed via the dispatch table's `error_bound_search` pointer:
  - **Keep `error_bound_j2k_compression`** unchanged: dynamically expands the `[cr_lo, cr_hi]` bracket (halving/doubling) because compression ratio is unbounded above, then bisects.
  - **Add `error_bound_jxl_compression`**: trivial binary search directly over `[dist_lo=0, dist_hi=25]` — bounds are known upfront, no bracket-expansion phase needed. Higher distance → more compression → more error; search for the highest distance that satisfies the error bound.
  - Add `emulate_jxl_compression` (mirrors `emulate_j2k_compression`): one encode–decode cycle at a given distance, returns the error quantile.
  - Call the selected backend's `error_bound_search` via `backend->error_bound_search(...)` in the generic encode path.
- [ ] Rename remaining `j2k_*` references in the generic encode/decode flow to `base_codec_*`.
- [ ] Update log messages to include selected backend name (e.g., `"[ebcc] using backend: jxl"`).

## Phase 4: API and Config Plumbing
- [ ] Add base compressor enum and rename quality parameter in codec config (`src/ebcc_codec.h`):
  ```c
  typedef enum { BASE_COMPRESSOR_J2K = 0, BASE_COMPRESSOR_JXL = 1 } base_compressor_t;
  ```
  - Rename `codec_config_t.base_cr` → `codec_config_t.base_param`. It holds compression ratio for J2K and `distance` for JXL; the backend interprets it.
  - Add `base_compressor_t base_compressor` field to `codec_config_t`, defaulting to `BASE_COMPRESSOR_J2K`.
- [ ] Register separate HDF5 filter IDs per backend (`src/h5z_ebcc.c`):
  - `H5Z_FILTER_EBCC_J2K = 308` — existing filter, no change to `cd_values` layout or `populate_config`.
  - `H5Z_FILTER_EBCC_JXL = 309` — new filter, identical `cd_values` layout (4 or 5 elements, same as J2K).
  - Each backend has its own `H5Z_class2_t` struct and a thin filter wrapper that sets `config.base_compressor` before calling `ebcc_encode`. `populate_config` is shared and unchanged.
  - `ebcc_decode` requires no caller-supplied backend — it reads the backend from the serialized header.
  - **`populate_config` is not modified.** The `cd_values` parameter layout remains `[height, width, base_param_bits, residual_type, error_bits?]` for both backends; `base_param_bits` is the IEEE 754 bit-pattern of the quality parameter (compression ratio for J2K, distance for JXL).
- [ ] Update Python wrapper (`ebcc/filter_wrapper.py`):
  - Add optional `base_compressor: str = "j2k"` parameter to `EBCC_Filter.__init__`.
  - J2K constructor keeps `base_cr: float`; JXL constructor uses `base_distance: float` instead — both are serialized identically into `cd_values[2]` as float bit-patterns.
  - Select `FILTER_ID` based on backend: `308` for `"j2k"`, `309` for `"jxl"`. `hdf_filter_opts` layout is unchanged.
  - When `base_compressor == "j2k"`, behavior is identical to the current code (full backward compatibility).
- [ ] Update Zarr wrapper (`ebcc/zarr_filter.py`):
  - Add `base_compressor` field to the ctypes `CodecConfigT` struct to match the updated C struct.
  - Use separate `codec_id` values per backend: `'ebcc_j2k'` (existing) and `'ebcc_jxl'` (new). The `arglist` layout is unchanged.
  - Update `get_config()` / `from_config()` to include `base_compressor`, defaulting to `"j2k"` for old configs that lack the field.
  - Test that old Zarr-encoded data (without `base_compressor` in config) still decodes correctly.

## Phase 5: Implement `libjxl` Backend
- [ ] Add JXL encode path for scaled 16-bit grayscale frames/tiles:
  - Use `JxlEncoderFrameSettings` with `distance` as the quality control.
  - Read `config->base_param` as `base_distance` and pass it directly to the JXL encoder — no conversion or binary search over `distance` needed.
- [ ] Add JXL decode path to reconstruct uint16 output, then apply existing min/max float rescaling.
- [ ] Implement `emulate_jxl_compression` and `error_bound_jxl_compression`:
  - `emulate_jxl_compression`: one JXL encode–decode cycle at a given `distance`, returns the error quantile.
  - `error_bound_jxl_compression`: binary search over `[0.0, 25.0]` — no bracket-expansion needed since distance is bounded on both ends. Finds the highest `distance` (most compression) satisfying the error bound.
- [ ] Validate functional equivalence in EBCC flow:
  - Works with `NONE`, `MAX_ERROR`, `RELATIVE_ERROR` residual modes.
  - Works with const-field optimization path (min == max → no base codec payload).
- [ ] Consider JXL `effort` parameter:
  - During error-bound binary search (which does many encode-decode iterations), use a lower `effort` (e.g., 3-5) for speed. For the final encode, optionally use higher `effort`.
  - Expose `effort` as an optional config parameter or environment variable.

## Phase 6: Tests
- [ ] Add parameterized tests for base compressor in:
  - `tests/test_netcdf.py` — add `@pytest.mark.parametrize("compressor", ["j2k", "jxl"])`.
  - `tests/test_zarr.py` — same parametrization.
- [ ] Required coverage for each backend (`j2k`, `jxl`):
  - Round-trip shape/dtype correctness.
  - Error-bound compliance (`max_error_target`, `relative_error_target`).
  - Constant field behavior (min == max).
  - Multi-dimensional data behavior.
- [ ] Add decode compatibility tests:
  - Save reference v1/J2K-encoded payloads as test fixtures.
  - Verify they decode correctly after v2 code changes.
  - Test that v2/JXL payloads fail gracefully when JXL is not compiled in.
- [ ] Keep benchmark suite usable; optionally add backend comparison dimension.

## Phase 7: CI and Release
- [ ] Update workflow to build and test with JXL on release targets (`.github/workflows/build.yml`):
  - linux amd64
  - linux arm64
  - macOS arm64
- [ ] Ensure wheel artifact includes all runtime dependencies needed by backend selection.
- [ ] Keep HDF5 example test step passing with default settings.
- [ ] Add a CI job that builds with `-DENABLE_JXL=OFF` to verify the optional-compilation path works.

## Phase 8: Docs and Migration Notes
- [ ] Update `README.md` user guide:
  - New `base_compressor` option.
  - Default remains JPEG2000.
  - Backward compatibility statement (old filter opts and old payloads still decode).
- [ ] Add examples for:
  - Python: `EBCC_Filter(..., base_compressor="jxl", base_distance=1.0)` — selects filter ID 309, same `cd_values` layout.
  - Zarr: `EBCCZarrFilter(..., base_compressor="jxl")` — uses `codec_id = 'ebcc_jxl'`, same `arglist`.
  - CDO: filter option layout is unchanged; only the leading filter ID (308 vs 309) differs.
- [ ] Document:
  - Build flags: `-DENABLE_JXL=ON/OFF`.
  - Any environment variables introduced for JXL (e.g., JXL effort).
  - Performance characteristics: expected speed/quality differences vs J2K.
  - The `base_distance` parameter range and effect on compression ratio / quality.

## Acceptance Criteria
- [ ] Existing J2K tests pass unchanged (or with explicit backend param defaulting to J2K).
- [ ] New JXL tests pass across supported platforms (linux amd64/arm64, macOS arm64).
- [ ] Old encoded data (legacy + v1 header) decode successfully with the new code.
- [ ] New v2 encoded data decodes correctly and enforces error bounds.
- [ ] Build with `-DENABLE_JXL=OFF` compiles and runs all J2K tests.
- [ ] Build pipeline can produce release artifacts with JXL support.

## Risks and Mitigations
- **Risk:** packaging native JXL dependencies across platforms is fragile.
  - Mitigation: decide dependency strategy early; validate CI matrix before deep code refactor. Make JXL optional at compile time so broken JXL builds don't block J2K users.
- **Risk:** silent compatibility break from a changed filter option layout — **eliminated by the plugin-ID approach.** `cd_values` layout is identical for both backends (4 or 5 elements); the backend is inferred entirely from the HDF5 filter ID (308 vs 309). Old filter configurations and old encoded data are unaffected.
- **Risk:** JXL encode/decode speed may be slower than J2K, amplified by binary search iterations.
  - Mitigation: use lower JXL `effort` during error-bound search iterations. Benchmark in Phase 0 spike.
- **Risk:** `libjxl` 16-bit grayscale support may have quirks or limitations.
  - Mitigation: Phase 0 feasibility spike validates this before any EBCC code changes.

## Suggested Execution Order
1. Feasibility spike: standalone libjxl test (Phase 0)
2. Dependency strategy decision + optional CMake flag (Phase 1)
3. Header v2 design: repurpose `reserved` byte (Phase 2)
4. Codec abstraction: dispatch table, refactor J2K into first backend (Phase 3)
5. JXL backend implementation (Phase 5)
6. Config plumbing: filter options, Python/Zarr wrappers (Phase 4)
7. Serialization + decode compatibility (Phase 2, finalize)
8. Tests and CI hardening (Phases 6/7)
9. Docs and release notes (Phase 8)
