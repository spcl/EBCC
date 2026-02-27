#!/usr/bin/env python3
import ctypes
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


_LOCAL_TMP_ROOT = Path(os.environ.get("TMPDIR", "/tmp")) / "ebcc_matplotlib"
_LOCAL_CACHE_DIR = _LOCAL_TMP_ROOT / "cache"
_LOCAL_MPL_DIR = _LOCAL_TMP_ROOT / "mplconfig"
_LOCAL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
_LOCAL_MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(_LOCAL_CACHE_DIR))
os.environ.setdefault("MPLCONFIGDIR", str(_LOCAL_MPL_DIR))

import matplotlib.pyplot as plt


DISTANCES = [0.01, 0.1, 0.5, 1.0, 2.0, 4.0]
BASE_COMPRESSOR_J2K = 0
BASE_COMPRESSOR_JXL = 1
RESIDUAL_NONE = 0


class CodecConfigT(ctypes.Structure):
    _fields_ = [
        ("dims", ctypes.c_size_t * 3),
        ("base_param", ctypes.c_float),
        ("base_compressor", ctypes.c_int),
        ("residual_compression_type", ctypes.c_int),
        ("residual_cr", ctypes.c_float),
        ("error", ctypes.c_float),
    ]


def float_to_u32(value: float) -> int:
    return int(np.array([value], dtype=np.float32).view(np.uint32)[0])


def find_codec_lib(root: Path) -> Path:
    candidates = [
        root / "ebcc" / "libh5z_ebcc_jxl.dylib",
        root / "src" / "build" / "lib" / "libh5z_ebcc_jxl.dylib",
        root / "ebcc" / "libh5z_ebcc.cpython-311-darwin.so",
        root / "libh5z_ebcc.cpython-311-darwin.so",
    ]
    for path in candidates:
        if path.exists():
            return path

    extra = sorted((root / "ebcc").glob("libh5z_ebcc*"))
    if extra:
        return extra[0]
    raise FileNotFoundError("Could not find EBCC shared library under project root")


def load_codec_lib(lib_path: Path) -> ctypes.CDLL:
    lib = ctypes.CDLL(str(lib_path))
    lib.populate_config.argtypes = [
        ctypes.POINTER(CodecConfigT),
        ctypes.c_size_t,
        np.ctypeslib.ndpointer(ctypes.c_uint, flags="C_CONTIGUOUS"),
        ctypes.c_size_t,
    ]
    lib.ebcc_encode.restype = ctypes.c_size_t
    lib.ebcc_encode.argtypes = [
        np.ctypeslib.ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
        ctypes.POINTER(CodecConfigT),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
    ]
    lib.ebcc_decode.restype = ctypes.c_size_t
    lib.ebcc_decode.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.POINTER(ctypes.c_float)),
    ]
    lib.free_buffer.argtypes = [ctypes.c_void_p]
    return lib


def build_cd_values(height: int, width: int, base_param: float) -> np.ndarray:
    return np.array(
        [height, width, float_to_u32(base_param), RESIDUAL_NONE],
        dtype=np.uint32,
    )


def encode(lib: ctypes.CDLL, data_flat: np.ndarray, cd_values: np.ndarray, compressor: int) -> bytes:
    config = CodecConfigT()
    out_buffer = ctypes.POINTER(ctypes.c_ubyte)()
    lib.populate_config(ctypes.byref(config), len(cd_values), cd_values, data_flat.nbytes)
    config.base_compressor = compressor
    encoded_size = lib.ebcc_encode(data_flat, ctypes.byref(config), ctypes.byref(out_buffer))
    out = ctypes.cast(out_buffer, ctypes.POINTER(ctypes.c_ubyte * encoded_size)).contents
    encoded = bytes(out)
    lib.free_buffer(out_buffer)
    return encoded


def decode(lib: ctypes.CDLL, encoded: bytes) -> np.ndarray:
    in_array_type = ctypes.c_ubyte * len(encoded)
    in_array = in_array_type.from_buffer_copy(encoded)
    out_buffer = ctypes.POINTER(ctypes.c_float)()
    decoded_size = lib.ebcc_decode(
        ctypes.cast(in_array, ctypes.c_void_p),
        len(encoded),
        ctypes.byref(out_buffer),
    )
    decoded = np.ctypeslib.as_array(out_buffer, shape=(decoded_size,)).copy()
    lib.free_buffer(out_buffer)
    return decoded


def float_to_uint16_range_scaled(data: np.ndarray) -> tuple[np.ndarray, float, float]:
    minval = float(np.min(data))
    maxval = float(np.max(data))
    rng = maxval - minval
    if rng <= 0.0:
        return np.zeros_like(data, dtype=np.uint16), minval, maxval
    scaled = np.rint((data - minval) * (65535.0 / rng)).astype(np.uint16)
    return scaled, minval, maxval


def uint16_to_float_range_scaled(data_u16: np.ndarray, minval: float, maxval: float) -> np.ndarray:
    rng = maxval - minval
    if rng <= 0.0:
        return np.full(data_u16.shape, minval, dtype=np.float32)
    return (data_u16.astype(np.float32) * (rng / 65535.0) + minval).astype(np.float32)


def normalize_to_unit_interval(data: np.ndarray) -> tuple[np.ndarray, float, float]:
    minval = float(np.min(data))
    maxval = float(np.max(data))
    rng = maxval - minval
    if rng <= 0.0:
        return np.zeros_like(data, dtype=np.float32), minval, maxval
    norm = ((data - minval) / rng).astype(np.float32)
    return norm, minval, maxval


def denormalize_from_unit_interval(data_norm: np.ndarray, minval: float, maxval: float) -> np.ndarray:
    rng = maxval - minval
    if rng <= 0.0:
        return np.full(data_norm.shape, minval, dtype=np.float32)
    return (data_norm.astype(np.float32) * rng + minval).astype(np.float32)


def write_pgm_u16(path: Path, data_u16: np.ndarray) -> None:
    if data_u16.ndim != 2:
        raise ValueError(f"PGM expects 2D array, got shape {data_u16.shape}")
    h, w = data_u16.shape
    with open(path, "wb") as f:
        f.write(f"P5\n{w} {h}\n65535\n".encode("ascii"))
        f.write(data_u16.astype(">u2", copy=False).tobytes())


def read_pgm_u16(path: Path) -> np.ndarray:
    with open(path, "rb") as f:
        magic = f.readline().strip()
        if magic != b"P5":
            raise ValueError(f"Unsupported PGM magic in {path}: {magic!r}")

        def _read_non_comment_line() -> bytes:
            while True:
                line = f.readline()
                if not line:
                    raise ValueError(f"Unexpected EOF while reading {path}")
                s = line.strip()
                if s and not s.startswith(b"#"):
                    return s

        dims = _read_non_comment_line()
        w, h = map(int, dims.split())
        maxval = int(_read_non_comment_line())
        if maxval != 65535:
            raise ValueError(f"Expected 16-bit PGM maxval=65535, got {maxval}")
        data = f.read()
    out = np.frombuffer(data, dtype=">u2")
    if out.size != h * w:
        raise ValueError(f"PGM payload size mismatch in {path}: {out.size} vs {h*w}")
    return out.reshape(h, w).copy()


def write_pfm_gray(path: Path, data: np.ndarray) -> None:
    if data.ndim != 2:
        raise ValueError(f"PFM expects 2D array, got shape {data.shape}")
    h, w = data.shape
    with open(path, "wb") as f:
        f.write(b"Pf\n")
        f.write(f"{w} {h}\n".encode("ascii"))
        # Negative scale means little-endian float payload.
        f.write(b"-1.0\n")
        f.write(data.astype("<f4", copy=False).tobytes())


def read_pfm_gray(path: Path) -> np.ndarray:
    with open(path, "rb") as f:
        magic = f.readline().strip()
        if magic != b"Pf":
            raise ValueError(f"Unsupported PFM magic in {path}: {magic!r}")
        dims = f.readline().strip()
        w, h = map(int, dims.split())
        scale = float(f.readline().strip())
        payload = f.read()
    dtype = ">f4" if scale > 0 else "<f4"
    out = np.frombuffer(payload, dtype=dtype)
    if out.size != h * w:
        raise ValueError(f"PFM payload size mismatch in {path}: {out.size} vs {h*w}")
    return out.reshape(h, w).astype(np.float32, copy=True)


def run_cjxl_roundtrip(input_path: Path, decoded_path: Path, distance: float) -> int:
    cjxl_bin = shutil.which("cjxl")
    djxl_bin = shutil.which("djxl")
    if cjxl_bin is None or djxl_bin is None:
        raise RuntimeError("Could not find cjxl/djxl in PATH")

    encoded_path = decoded_path.with_suffix(".jxl")
    subprocess.run(
        [cjxl_bin, str(input_path), str(encoded_path), "--distance", str(distance), "--quiet", "--effort", "10", "-x", "color_space=Gra_D65_Rel_Lin"], #"--disable_perceptual_optimizations" # cause crash
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        [djxl_bin, str(encoded_path), str(decoded_path), "--quiet"],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return encoded_path.stat().st_size


def match_j2k_to_target_size(
    lib: ctypes.CDLL,
    data_flat: np.ndarray,
    height: int,
    width: int,
    target_size_bytes: int,
) -> tuple[float, bytes]:
    cache: dict[float, bytes] = {}

    def encode_j2k_for_cr(base_cr: float) -> bytes:
        base_cr = float(max(1.0, min(1000.0, base_cr)))
        key = round(base_cr, 6)
        if key not in cache:
            cd_values = build_cd_values(height, width, base_cr)
            cache[key] = encode(lib, data_flat, cd_values, BASE_COMPRESSOR_J2K)
        return cache[key]

    lo, hi = 1.0, 1000.0
    enc_lo = encode_j2k_for_cr(lo)
    enc_hi = encode_j2k_for_cr(hi)
    size_lo, size_hi = len(enc_lo), len(enc_hi)

    if target_size_bytes >= size_lo:
        return lo, enc_lo
    if target_size_bytes <= size_hi:
        return hi, enc_hi

    best_cr, best_enc = lo, enc_lo
    best_gap = abs(size_lo - target_size_bytes)

    for _ in range(20):
        mid = 0.5 * (lo + hi)
        enc_mid = encode_j2k_for_cr(mid)
        size_mid = len(enc_mid)
        gap = abs(size_mid - target_size_bytes)
        if gap < best_gap:
            best_gap = gap
            best_cr = mid
            best_enc = enc_mid
        if size_mid > target_size_bytes:
            lo = mid
        else:
            hi = mid

    # Check final boundaries as well.
    for cr in (lo, hi):
        enc = encode_j2k_for_cr(cr)
        gap = abs(len(enc) - target_size_bytes)
        if gap < best_gap:
            best_gap = gap
            best_cr = cr
            best_enc = enc

    return best_cr, best_enc


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    data_path = root / "data" / "test_data.npy"
    out_path = root / "data" / "jxl_compression_histogram.png"

    original = np.load(data_path).astype(np.float32, copy=False)
    if original.ndim != 2:
        raise ValueError(f"Expected 2D array in {data_path}, got shape {original.shape}")

    height, width = map(int, original.shape)
    original_flat = np.ascontiguousarray(original.ravel())
    data_range = float(np.max(original_flat) - np.min(original_flat))

    lib_path = find_codec_lib(root)
    lib = load_codec_lib(lib_path)
    print(f"Using codec library: {lib_path}")

    top_ebcc_jxl_errors = {}
    top_ebcc_jxl_sizes = {}
    top_ebcc_jxl_cr = {}
    for distance in DISTANCES:
        cd_values = build_cd_values(height, width, distance)
        encoded = encode(lib, original_flat, cd_values, BASE_COMPRESSOR_JXL)
        decoded = decode(lib, encoded)
        if decoded.size != original_flat.size:
            raise RuntimeError(
                f"Decoded size mismatch at distance={distance}: "
                f"{decoded.size} vs {original_flat.size}"
            )
        err = decoded - original_flat
        top_ebcc_jxl_errors[distance] = err
        top_ebcc_jxl_sizes[distance] = len(encoded)
        top_ebcc_jxl_cr[distance] = original_flat.nbytes / float(len(encoded))
        max_abs = float(np.max(np.abs(err)))
        mean_abs = float(np.mean(np.abs(err)))
        print(
            f"distance={distance:>4g}  encoded_bytes={len(encoded):>8d}  "
            f"max_abs_err={max_abs:.6g}  mean_abs_err={mean_abs:.6g}"
        )

    top_cjxl_pgm_errors = {}
    top_cjxl_pgm_sizes = {}
    top_cjxl_pgm_cr = {}
    top_cjxl_pfm_errors = {}
    top_cjxl_pfm_sizes = {}
    top_cjxl_pfm_cr = {}

    pgm_u16, pgm_min, pgm_max = float_to_uint16_range_scaled(original)
    pfm_norm, pfm_min, pfm_max = normalize_to_unit_interval(original)
    pgm_source_bytes = int(pgm_u16.nbytes)
    pfm_source_bytes = int(original.nbytes)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        pgm_input = tmp / "input.pgm"
        pfm_input = tmp / "input.pfm"
        write_pgm_u16(pgm_input, pgm_u16)
        write_pfm_gray(pfm_input, pfm_norm)

        for distance in DISTANCES:
            pgm_decoded_path = tmp / f"decoded_pgm_d{distance:g}.pgm"
            pgm_encoded_size = run_cjxl_roundtrip(pgm_input, pgm_decoded_path, distance)
            pgm_decoded_u16 = read_pgm_u16(pgm_decoded_path)
            pgm_decoded_f32 = uint16_to_float_range_scaled(pgm_decoded_u16, pgm_min, pgm_max)
            pgm_err = np.ascontiguousarray(pgm_decoded_f32.ravel() - original_flat)
            top_cjxl_pgm_errors[distance] = pgm_err
            top_cjxl_pgm_sizes[distance] = pgm_encoded_size
            top_cjxl_pgm_cr[distance] = pgm_source_bytes / float(pgm_encoded_size)

            pfm_decoded_path = tmp / f"decoded_pfm_d{distance:g}.pfm"
            pfm_encoded_size = run_cjxl_roundtrip(pfm_input, pfm_decoded_path, distance)
            pfm_decoded_norm = read_pfm_gray(pfm_decoded_path)
            pfm_decoded_f32 = denormalize_from_unit_interval(pfm_decoded_norm, pfm_min, pfm_max)
            pfm_err = np.ascontiguousarray(pfm_decoded_f32.ravel() - original_flat)
            top_cjxl_pfm_errors[distance] = pfm_err
            top_cjxl_pfm_sizes[distance] = pfm_encoded_size
            top_cjxl_pfm_cr[distance] = pfm_source_bytes / float(pfm_encoded_size)

            print(
                f"distance={distance:>4g}  cjxl_pgm_cr={top_cjxl_pgm_cr[distance]:8.3f}  "
                f"cjxl_pfm_cr={top_cjxl_pfm_cr[distance]:8.3f}"
            )

    def build_j2k_counterpart(
        target_sizes: dict[float, int],
        ratio_numerator_bytes: int,
        name: str,
    ) -> tuple[dict[float, np.ndarray], dict[float, float]]:
        errors = {}
        ratios = {}
        for distance in DISTANCES:
            target_size = target_sizes[distance]
            matched_cr_param, encoded_j2k = match_j2k_to_target_size(
                lib, original_flat, height, width, target_size
            )
            decoded_j2k = decode(lib, encoded_j2k)
            if decoded_j2k.size != original_flat.size:
                raise RuntimeError(
                    f"Decoded J2K size mismatch at distance={distance} for {name}: "
                    f"{decoded_j2k.size} vs {original_flat.size}"
                )
            errors[distance] = decoded_j2k - original_flat
            ratios[distance] = ratio_numerator_bytes / float(len(encoded_j2k))
            print(
                f"{name} distance={distance:>4g}  target_bytes={target_size:8d}  "
                f"matched_j2k_cr_param={matched_cr_param:8.3f}  "
                f"achieved_cr={ratios[distance]:8.3f}"
            )
        return errors, ratios

    bot_j2k_for_ebcc_jxl_errors, bot_j2k_for_ebcc_jxl_cr = build_j2k_counterpart(
        top_ebcc_jxl_sizes,
        original_flat.nbytes,
        "j2k_vs_ebcc_jxl",
    )
    bot_j2k_for_cjxl_pgm_errors, bot_j2k_for_cjxl_pgm_cr = build_j2k_counterpart(
        top_cjxl_pgm_sizes,
        pgm_source_bytes,
        "j2k_vs_cjxl_pgm",
    )
    bot_j2k_for_cjxl_pfm_errors, bot_j2k_for_cjxl_pfm_cr = build_j2k_counterpart(
        top_cjxl_pfm_sizes,
        pfm_source_bytes,
        "j2k_vs_cjxl_pfm",
    )

    global_min = min(
        min(float(np.min(err)) for err in top_ebcc_jxl_errors.values()),
        min(float(np.min(err)) for err in top_cjxl_pgm_errors.values()),
        min(float(np.min(err)) for err in top_cjxl_pfm_errors.values()),
        min(float(np.min(err)) for err in bot_j2k_for_ebcc_jxl_errors.values()),
        min(float(np.min(err)) for err in bot_j2k_for_cjxl_pgm_errors.values()),
        min(float(np.min(err)) for err in bot_j2k_for_cjxl_pfm_errors.values()),
    )
    global_max = max(
        max(float(np.max(err)) for err in top_ebcc_jxl_errors.values()),
        max(float(np.max(err)) for err in top_cjxl_pgm_errors.values()),
        max(float(np.max(err)) for err in top_cjxl_pfm_errors.values()),
        max(float(np.max(err)) for err in bot_j2k_for_ebcc_jxl_errors.values()),
        max(float(np.max(err)) for err in bot_j2k_for_cjxl_pgm_errors.values()),
        max(float(np.max(err)) for err in bot_j2k_for_cjxl_pfm_errors.values()),
    )
    bins = np.linspace(global_min, global_max, 120)

    fig, axes = plt.subplots(2, 3, figsize=(26, 12), sharex=True, sharey=True)
    ax_top_ebcc_jxl = axes[0, 0]
    ax_top_cjxl_pgm = axes[0, 1]
    ax_top_cjxl_pfm = axes[0, 2]
    ax_bot_j2k_ebcc_jxl = axes[1, 0]
    ax_bot_j2k_cjxl_pgm = axes[1, 1]
    ax_bot_j2k_cjxl_pfm = axes[1, 2]
    cmap = plt.get_cmap("tab10")

    for idx, distance in enumerate(DISTANCES):
        ax_top_ebcc_jxl.hist(
            top_ebcc_jxl_errors[distance],
            bins=bins,
            histtype="step",
            linewidth=1.8,
            density=True,
            color=cmap(idx),
            label=f"d={distance:g}, CR={top_ebcc_jxl_cr[distance]:.2f}",
        )
        ax_top_cjxl_pgm.hist(
            top_cjxl_pgm_errors[distance],
            bins=bins,
            histtype="step",
            linewidth=1.8,
            density=True,
            color=cmap(idx),
            label=f"d={distance:g}, CR={top_cjxl_pgm_cr[distance]:.2f}",
        )
        ax_top_cjxl_pfm.hist(
            top_cjxl_pfm_errors[distance],
            bins=bins,
            histtype="step",
            linewidth=1.8,
            density=True,
            color=cmap(idx),
            label=f"d={distance:g}, CR={top_cjxl_pfm_cr[distance]:.2f}",
        )
        ax_bot_j2k_ebcc_jxl.hist(
            bot_j2k_for_ebcc_jxl_errors[distance],
            bins=bins,
            histtype="step",
            linewidth=1.8,
            density=True,
            color=cmap(idx),
            label=f"d={distance:g}, CR={bot_j2k_for_ebcc_jxl_cr[distance]:.2f}",
        )
        ax_bot_j2k_cjxl_pgm.hist(
            bot_j2k_for_cjxl_pgm_errors[distance],
            bins=bins,
            histtype="step",
            linewidth=1.8,
            density=True,
            color=cmap(idx),
            label=f"d={distance:g}, CR={bot_j2k_for_cjxl_pgm_cr[distance]:.2f}",
        )
        ax_bot_j2k_cjxl_pfm.hist(
            bot_j2k_for_cjxl_pfm_errors[distance],
            bins=bins,
            histtype="step",
            linewidth=1.8,
            density=True,
            color=cmap(idx),
            label=f"d={distance:g}, CR={bot_j2k_for_cjxl_pfm_cr[distance]:.2f}",
        )

    ax_top_ebcc_jxl.set_title("EBCC Pure JXL")
    ax_top_cjxl_pgm.set_title("cjxl (PGM uint16 input)")
    ax_top_cjxl_pfm.set_title("cjxl (PFM float input, normalized [0,1])")
    ax_bot_j2k_ebcc_jxl.set_title("EBCC JP2 Counterpart (match top-left CR)")
    ax_bot_j2k_cjxl_pgm.set_title("EBCC JP2 Counterpart (match top-middle CR)")
    ax_bot_j2k_cjxl_pfm.set_title("EBCC JP2 Counterpart (match top-right CR)")
    all_axes = (
        ax_top_ebcc_jxl,
        ax_top_cjxl_pgm,
        ax_top_cjxl_pfm,
        ax_bot_j2k_ebcc_jxl,
        ax_bot_j2k_cjxl_pgm,
        ax_bot_j2k_cjxl_pfm,
    )
    for ax in all_axes:
        ax.set_xlabel("Decoded - Original")
        ax.set_yscale("log")
        ax.grid(alpha=0.25, linestyle="--")

    # Range-relative absolute error thresholds: 0.1% and 0.5% of data range.
    rel_01 = 0.001 * data_range
    rel_05 = 0.005 * data_range
    for ax in all_axes:
        ax.axvline(rel_01, color="black", linestyle="--", linewidth=1.3)
        ax.axvline(-rel_01, color="black", linestyle="--", linewidth=1.3)
        ax.axvline(rel_05, color="gray", linestyle="--", linewidth=1.3)
        ax.axvline(-rel_05, color="gray", linestyle="--", linewidth=1.3)

    ax_top_ebcc_jxl.set_ylabel("Density (log scale)")
    ax_bot_j2k_ebcc_jxl.set_ylabel("Density (log scale)")
    for ax in all_axes:
        ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(f"Saved histogram: {out_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise
