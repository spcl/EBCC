#!/usr/bin/env python3
import ctypes
from pathlib import Path

import netCDF4 as nc
import numpy as np


RELATIVE_ERROR_BOUND = 0.01
BASE_CR = 20.0
RESIDUAL_RELATIVE_ERROR = 2  # enum residual_t::RELATIVE_ERROR


class CodecConfigT(ctypes.Structure):
    _fields_ = [
        ("dims", ctypes.c_size_t * 3),
        ("base_cr", ctypes.c_float),
        ("residual_compression_type", ctypes.c_int),
        ("residual_cr", ctypes.c_float),
        ("error", ctypes.c_float),
    ]


def float_to_u32(value: float) -> int:
    return int(np.array([value], dtype=np.float32).view(np.uint32)[0])


def load_temperature_array(root: Path) -> np.ndarray:
    path = root / "data" / "temperature.nc"
    ds = nc.Dataset(path, "r")
    try:
        arr = ds.variables["temperature"][:]
    finally:
        ds.close()
    return np.ascontiguousarray(arr.astype(np.float32).ravel())


def build_cd_values(height: int, width: int, base_cr: float, relative_error: float) -> np.ndarray:
    return np.array(
        [
            height,
            width,
            float_to_u32(base_cr),
            RESIDUAL_RELATIVE_ERROR,
            float_to_u32(relative_error),
        ],
        dtype=np.uint32,
    )


def load_codec_lib(path: Path) -> ctypes.CDLL:
    lib = ctypes.CDLL(str(path))
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
    if hasattr(lib, "ebcc_decode_legacy"):
        lib.ebcc_decode_legacy.restype = ctypes.c_size_t
        lib.ebcc_decode_legacy.argtypes = [
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.POINTER(ctypes.c_float)),
        ]
    lib.free_buffer.argtypes = [ctypes.c_void_p]
    return lib


def encode(lib: ctypes.CDLL, data: np.ndarray, cd_values: np.ndarray) -> bytes:
    config = CodecConfigT()
    out_buffer = ctypes.POINTER(ctypes.c_ubyte)()
    lib.populate_config(ctypes.byref(config), len(cd_values), cd_values, data.nbytes)
    encoded_size = lib.ebcc_encode(data, ctypes.byref(config), ctypes.byref(out_buffer))

    out = ctypes.cast(out_buffer, ctypes.POINTER(ctypes.c_ubyte * encoded_size)).contents
    encoded = bytes(out)
    lib.free_buffer(out_buffer)
    return encoded


def decode(lib: ctypes.CDLL, encoded: bytes) -> np.ndarray:
    in_array_type = ctypes.c_ubyte * len(encoded)
    in_array = in_array_type.from_buffer_copy(encoded)
    out_buffer = ctypes.POINTER(ctypes.c_float)()
    decoded_size = lib.ebcc_decode(ctypes.cast(in_array, ctypes.c_void_p), len(encoded), ctypes.byref(out_buffer))
    decoded = np.ctypeslib.as_array(out_buffer, shape=(decoded_size,)).copy()
    lib.free_buffer(out_buffer)
    return decoded


def decode_legacy(lib: ctypes.CDLL, encoded: bytes) -> np.ndarray:
    if not hasattr(lib, "ebcc_decode_legacy"):
        raise RuntimeError("Library does not export ebcc_decode_legacy")
    in_array_type = ctypes.c_ubyte * len(encoded)
    in_array = in_array_type.from_buffer_copy(encoded)
    out_buffer = ctypes.POINTER(ctypes.c_float)()
    decoded_size = lib.ebcc_decode_legacy(ctypes.cast(in_array, ctypes.c_void_p), len(encoded), ctypes.byref(out_buffer))
    decoded = np.ctypeslib.as_array(out_buffer, shape=(decoded_size,)).copy()
    lib.free_buffer(out_buffer)
    return decoded


def max_relative_error(original: np.ndarray, decoded: np.ndarray) -> float:
    data_range = float(np.max(original) - np.min(original))
    if data_range == 0.0:
        return 0.0
    return float(np.max(np.abs(original - decoded)) / data_range)


def assert_error_within_bound(name: str, original: np.ndarray, decoded: np.ndarray, bound: float) -> None:
    err = max_relative_error(original, decoded)
    print(f"{name}: max_relative_error={err:.8f} (bound={bound:.8f})")
    assert err <= bound + 1e-7, f"{name} exceeded bound: {err} > {bound}"


def report_same(name_a: str, arr_a: np.ndarray, name_b: str, arr_b: np.ndarray) -> None:
    equal = np.array_equal(arr_a, arr_b)
    all_close = np.allclose(arr_a, arr_b, rtol=0.0, atol=0.0)
    max_abs_diff = float(np.max(np.abs(arr_a - arr_b)))
    print(
        f"compare {name_a} vs {name_b}: "
        f"array_equal={equal}, all_close={all_close}, max_abs_diff={max_abs_diff:.12e}"
    )


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    data = load_temperature_array(root)
    height, width = 721, 1440
    cd_values = build_cd_values(height, width, BASE_CR, RELATIVE_ERROR_BOUND)

    legacy_lib = load_codec_lib(root / "data" / "test_legacy_decode" / "libh5z_ebcc_legacy.dylib")
    new_lib = load_codec_lib(root / "data" / "test_legacy_decode" / "libh5z_ebcc_newheader.dylib")

    legacy_encoded = encode(legacy_lib, data, cd_values)
    legacy_decoded_old = decode(legacy_lib, legacy_encoded)
    legacy_decoded_new = decode(new_lib, legacy_encoded)
    legacy_decoded_new_legacy = decode_legacy(new_lib, legacy_encoded)

    assert_error_within_bound("legacy_encode->legacy_decode", data, legacy_decoded_old, RELATIVE_ERROR_BOUND)
    assert_error_within_bound("legacy_encode->new_decode", data, legacy_decoded_new, RELATIVE_ERROR_BOUND)
    assert_error_within_bound("legacy_encode->new_decode_legacy", data, legacy_decoded_new_legacy, RELATIVE_ERROR_BOUND)
    report_same("legacy_decode", legacy_decoded_old, "new_decode_on_legacy_data", legacy_decoded_new)
    report_same("legacy_decode", legacy_decoded_old, "new_decode_legacy_on_legacy_data", legacy_decoded_new_legacy)

    new_encoded = encode(new_lib, data, cd_values)
    new_decoded_new = decode(new_lib, new_encoded)
    assert_error_within_bound("new_encode->new_decode", data, new_decoded_new, RELATIVE_ERROR_BOUND)

    report_same("legacy_encode->legacy_decode", legacy_decoded_old, "new_encode->new_decode", new_decoded_new)

    print("PASS")


if __name__ == "__main__":
    main()
