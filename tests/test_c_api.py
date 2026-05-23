import ctypes

import numpy as np
import pytest
from numpy.ctypeslib import ndpointer

try:
    from ebcc import EBCC_FILTER_PATH
    from ebcc.zarr_filter import CodecConfigT
except (FileNotFoundError, ImportError) as exc:
    pytest.skip(f"EBCC C API library is not available: {exc}", allow_module_level=True)


MAX_ERROR = 1
NDIMS = 3


class ChunkingHeader(ctypes.Structure):
    _fields_ = [
        ("magic", ctypes.c_ubyte * 4),
        ("version", ctypes.c_uint32),
        ("ndims", ctypes.c_uint32),
        ("dims", ctypes.c_uint64 * NDIMS),
        ("chunk_dims", ctypes.c_uint64 * NDIMS),
        ("num_chunks", ctypes.c_uint64),
        ("chunk_size", ctypes.c_uint64),
    ]


@pytest.fixture(scope="module")
def ebcc_lib():
    lib = ctypes.CDLL(EBCC_FILTER_PATH)
    lib.ebcc_encode.restype = ctypes.c_size_t
    lib.ebcc_encode.argtypes = [
        ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
        ctypes.POINTER(CodecConfigT),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
    ]
    lib.ebcc_decode_chunking.restype = ctypes.c_size_t
    lib.ebcc_decode_chunking.argtypes = [
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.POINTER(ctypes.c_float)),
    ]
    lib.ebcc_encode_chunking.restype = ctypes.c_size_t
    lib.ebcc_encode_chunking.argtypes = [
        ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
        ctypes.POINTER(CodecConfigT),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
    ]
    lib.free_buffer.argtypes = [ctypes.c_void_p]
    lib.free_buffer.restype = None
    return lib


def make_config(shape, chunk_shape=None, *, base_cr=2.0, error=0.01):
    config = CodecConfigT()
    config.dims[:] = shape
    config.base_cr = base_cr
    config.residual_compression_type = MAX_ERROR
    config.residual_cr = 0.0
    config.error = error
    config.chunk_dims[:] = chunk_shape or (0, 0, 0)
    return config


def make_data(shape):
    indices = np.indices(shape, dtype=np.float32)
    data = indices[0] * 100.0 + indices[1] * 1.5 + indices[2] * 0.25
    return np.ascontiguousarray(data, dtype=np.float32)


def encode_chunking(lib, data, config):
    out_buffer = ctypes.POINTER(ctypes.c_ubyte)()
    flattened = np.ascontiguousarray(data, dtype=np.float32).ravel()
    encoded_size = lib.ebcc_encode_chunking(
        flattened,
        ctypes.byref(config),
        ctypes.byref(out_buffer),
    )
    assert encoded_size > 0
    assert out_buffer

    try:
        out_array = ctypes.cast(
            out_buffer,
            ctypes.POINTER(ctypes.c_ubyte * encoded_size),
        ).contents
        return bytes(out_array)
    finally:
        lib.free_buffer(out_buffer)


def encode_plain(lib, data, config):
    out_buffer = ctypes.POINTER(ctypes.c_ubyte)()
    flattened = np.ascontiguousarray(data, dtype=np.float32).ravel()
    encoded_size = lib.ebcc_encode(
        flattened,
        ctypes.byref(config),
        ctypes.byref(out_buffer),
    )
    assert encoded_size > 0
    assert out_buffer

    try:
        out_array = ctypes.cast(
            out_buffer,
            ctypes.POINTER(ctypes.c_ubyte * encoded_size),
        ).contents
        return bytes(out_array)
    finally:
        lib.free_buffer(out_buffer)


def decode_chunking(lib, encoded, shape):
    out_buffer = ctypes.POINTER(ctypes.c_float)()
    encoded_array_type = ctypes.c_ubyte * len(encoded)
    encoded_buffer = encoded_array_type.from_buffer_copy(encoded)
    decoded_size = lib.ebcc_decode_chunking(
        ctypes.byref(encoded_buffer),
        len(encoded),
        ctypes.byref(out_buffer),
    )
    assert decoded_size == int(np.prod(shape))
    assert out_buffer

    try:
        decoded = np.ctypeslib.as_array(out_buffer, shape=(decoded_size,))
        return decoded.copy().reshape(shape)
    finally:
        lib.free_buffer(out_buffer)


def read_chunking_header(encoded):
    assert len(encoded) >= ctypes.sizeof(ChunkingHeader)
    return ChunkingHeader.from_buffer_copy(encoded)


def assert_roundtrip_close(decoded, data):
    assert decoded.shape == data.shape
    assert decoded.dtype == data.dtype
    assert np.allclose(decoded, data, atol=0.02)


def test_chunking_header_and_slab_roundtrip(ebcc_lib):
    shape = (2, 32, 32)
    chunk_shape = (1, 32, 32)
    data = make_data(shape)

    encoded = encode_chunking(ebcc_lib, data, make_config(shape, chunk_shape))
    header = read_chunking_header(encoded)

    assert bytes(header.magic) == b"EBCK"
    assert header.version == 1
    assert header.ndims == NDIMS
    assert tuple(header.dims) == shape
    assert tuple(header.chunk_dims) == chunk_shape
    assert header.num_chunks == 2
    assert header.chunk_size == int(np.prod(chunk_shape))

    decoded = decode_chunking(ebcc_lib, encoded, shape)
    assert_roundtrip_close(decoded, data)


def test_chunking_roundtrip_with_padded_edge_chunks(ebcc_lib):
    shape = (3, 33, 35)
    chunk_shape = (2, 32, 32)
    data = make_data(shape)

    encoded = encode_chunking(ebcc_lib, data, make_config(shape, chunk_shape))
    header = read_chunking_header(encoded)

    assert tuple(header.dims) == shape
    assert tuple(header.chunk_dims) == chunk_shape
    assert header.num_chunks == 8
    assert header.chunk_size == int(np.prod(chunk_shape))

    decoded = decode_chunking(ebcc_lib, encoded, shape)
    assert_roundtrip_close(decoded, data)


def test_contiguous_chunk_dim_larger_than_data_dim_roundtrip(ebcc_lib):
    shape = (2, 32, 32)
    chunk_shape = (4, 32, 32)
    data = make_data(shape)

    encoded = encode_chunking(ebcc_lib, data, make_config(shape, chunk_shape))
    header = read_chunking_header(encoded)

    assert tuple(header.dims) == shape
    assert tuple(header.chunk_dims) == chunk_shape
    assert header.num_chunks == 1
    assert header.chunk_size == int(np.prod(chunk_shape))

    decoded = decode_chunking(ebcc_lib, encoded, shape)
    assert_roundtrip_close(decoded, data)


def test_non_contiguous_chunk_dim_larger_than_data_dim_roundtrip(ebcc_lib):
    shape = (2, 33, 35)
    chunk_shape = (1, 64, 64)
    data = make_data(shape)

    encoded = encode_chunking(ebcc_lib, data, make_config(shape, chunk_shape))
    header = read_chunking_header(encoded)

    assert tuple(header.dims) == shape
    assert tuple(header.chunk_dims) == chunk_shape
    assert header.num_chunks == 2
    assert header.chunk_size == int(np.prod(chunk_shape))

    decoded = decode_chunking(ebcc_lib, encoded, shape)
    assert_roundtrip_close(decoded, data)


def test_zero_chunk_dims_default_to_full_array_chunk(ebcc_lib):
    shape = (2, 32, 32)
    data = make_data(shape)

    encoded = encode_chunking(ebcc_lib, data, make_config(shape))
    header = read_chunking_header(encoded)

    assert tuple(header.chunk_dims) == shape
    assert header.num_chunks == 1
    assert header.chunk_size == int(np.prod(shape))

    decoded = decode_chunking(ebcc_lib, encoded, shape)
    assert_roundtrip_close(decoded, data)


def test_decode_chunking_accepts_plain_ebcc_payload(ebcc_lib):
    shape = (2, 32, 32)
    data = make_data(shape)
    config = make_config(shape)

    encoded = encode_plain(ebcc_lib, data, config)
    assert not encoded.startswith(b"EBCK")

    decoded = decode_chunking(ebcc_lib, encoded, shape)
    assert_roundtrip_close(decoded, data)
