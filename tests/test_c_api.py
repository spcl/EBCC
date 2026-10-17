import ctypes
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import numpy as np
import pytest
from numpy.ctypeslib import ndpointer

try:
    from ebcc import EBCC_FILTER_PATH
    from ebcc.zarr_filter import CodecConfigT
except (FileNotFoundError, ImportError) as exc:
    pytest.skip(f"EBCC C API library is not available: {exc}", allow_module_level=True)


FIXED_CR = 0
MAX_ERROR = 1
RELATIVE_ERROR = 2
NDIMS = 3


class EbccHeader(ctypes.Structure):
    _fields_ = [
        ("magic", ctypes.c_ubyte * 4),
        ("version", ctypes.c_ubyte),
        ("flags", ctypes.c_ubyte),
        ("reserved", ctypes.c_uint16),
        ("minval_bits", ctypes.c_uint32),
        ("maxval_bits", ctypes.c_uint32),
        ("coeffs_size", ctypes.c_uint64),
        ("residual_minval_bits", ctypes.c_uint32),
        ("residual_maxval_bits", ctypes.c_uint32),
        ("compressed_size", ctypes.c_uint64),
        ("tail_size", ctypes.c_uint64),
    ]


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
    lib = ctypes.CDLL(os.environ.get("EBCC_TEST_LIBRARY", EBCC_FILTER_PATH))
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
    lib.ebcc_encode_chunking_compat.restype = ctypes.c_size_t
    lib.ebcc_encode_chunking_compat.argtypes = [
        ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
        ctypes.POINTER(CodecConfigT),
        ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
    ]
    lib.free_buffer.argtypes = [ctypes.c_void_p]
    lib.free_buffer.restype = None
    return lib


def make_config(shape, chunk_shape=None, *, base_cr=2.0, error=0.01, residual_type=MAX_ERROR):
    config = CodecConfigT()
    config.dims[:] = shape
    config.base_cr = base_cr
    config.residual_compression_type = residual_type
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


def encode_chunking_compat(lib, data, config):
    out_buffer = ctypes.POINTER(ctypes.c_ubyte)()
    flattened = np.ascontiguousarray(data, dtype=np.float32).ravel()
    encoded_size = lib.ebcc_encode_chunking_compat(
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


def codestream_transform(codestream):
    marker_offset = codestream.index(b"\xff\x52")
    segment_length = int.from_bytes(codestream[marker_offset + 2:marker_offset + 4], "big")
    cod_payload = codestream[marker_offset + 4:marker_offset + 2 + segment_length]
    return cod_payload[9]


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


@pytest.mark.parametrize('residual_type', [FIXED_CR, MAX_ERROR])
def test_chunking_roundtrip_with_padded_edge_chunks(ebcc_lib, residual_type):
    shape = (3, 33, 35)
    chunk_shape = (2, 32, 32)
    data = make_data(shape)

    encoded = encode_chunking(ebcc_lib, data, make_config(shape, chunk_shape, residual_type=residual_type))
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


def test_chunking_compat_zero_chunk_dims_default_to_internal_tiles(ebcc_lib):
    shape = (1, 2048, 32)
    data = make_data(shape)

    encoded = encode_chunking_compat(
        ebcc_lib,
        data,
        make_config(shape, residual_type=FIXED_CR),
    )
    header = read_chunking_header(encoded)

    assert tuple(header.chunk_dims) == (1, 1024, 32)
    assert header.num_chunks == 2
    assert header.chunk_size == 1024 * 32


def test_chunking_compat_accepts_relative_error_bound(ebcc_lib):
    shape = (2, 32, 32)
    data = make_data(shape)

    encoded = encode_chunking_compat(
        ebcc_lib,
        data,
        make_config(shape, error=0.01, residual_type=RELATIVE_ERROR),
    )
    header = read_chunking_header(encoded)

    assert tuple(header.chunk_dims) == (2, 32, 32)
    assert header.num_chunks == 1


def test_decode_chunking_accepts_plain_ebcc_payload(ebcc_lib):
    shape = (2, 32, 32)
    data = make_data(shape)
    config = make_config(shape)

    encoded = encode_plain(ebcc_lib, data, config)
    assert not encoded.startswith(b"EBCK")

    decoded = decode_chunking(ebcc_lib, encoded, shape)
    assert_roundtrip_close(decoded, data)


@pytest.mark.parametrize('truncate', [False, True])
def test_chunked_decode_preserves_caller_buffer(ebcc_lib, truncate):
    shape = (3, 67, 99)
    data = make_data(shape)
    config = make_config(shape, (2, 40, 64), error=.01)
    payload = encode_chunking(ebcc_lib, data, config)
    if truncate:
        payload = payload[:-1]
    storage = ctypes.create_string_buffer(payload)
    output = np.empty(shape, dtype=np.float32)
    pointer = output.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    size = ebcc_lib.ebcc_decode_chunking(storage, len(payload), ctypes.byref(pointer))
    assert ctypes.addressof(pointer.contents) == output.ctypes.data
    assert size == (0 if truncate else data.size)
    if not truncate:
        assert np.max(np.abs(output.astype(np.float64) - data)) <= config.error


def test_chunked_decode_rejects_null_arguments(ebcc_lib):
    output = ctypes.POINTER(ctypes.c_float)()
    assert ebcc_lib.ebcc_decode_chunking(None, 80, ctypes.byref(output)) == 0
    assert not output
    assert ebcc_lib.ebcc_decode_chunking(ctypes.create_string_buffer(b'EBCK'), 4, None) == 0


@pytest.mark.parametrize('components', [1, 3])
@pytest.mark.parametrize('disable_mct', [False, True])
def test_fixed_cr_mode_uses_base_compression_ratio(ebcc_lib, components, disable_mct):
    shape = (components, 128, 128)
    data = np.random.default_rng(71).normal(size=shape).astype(np.float32)
    sizes, errors = [], []
    # Error-bound allocation and ablation switches must not override FIXED_CR mode.
    # patch.dict restores the original environment when this block exits.
    with patch.dict(os.environ, {
        'EBCC_DISABLE_MCT': str(int(disable_mct)),
        'EBCC_KAPPA': 'invalid',
        'EBCC_INTERNAL_FIXED_QUANTIZATION': '1',
        'EBCC_INTERNAL_SAMPLE_RESIDUAL': '1',
    }):
        for ratio in (4, 16):
            config = make_config(shape, base_cr=ratio, residual_type=FIXED_CR, error=float('nan'))
            encoded = encode_plain(ebcc_lib, data, config)
            header = EbccHeader.from_buffer_copy(encoded)
            assert bytes(header.magic) == b'EBCC' and header.version == 2
            assert header.flags == 0 and header.reserved == components
            assert header.coeffs_size == header.residual_minval_bits == header.residual_maxval_bits == 0
            # Volumes retain an empty residual directory in the existing format.
            assert header.compressed_size == (24 * components if components > 1 else 0)
            start = ctypes.sizeof(EbccHeader)
            assert encoded[start:start + header.compressed_size] == bytes(header.compressed_size)
            assert len(encoded) == start + header.compressed_size + header.tail_size
            assert data.nbytes / header.tail_size == pytest.approx(ratio, rel=.15)
            decoded = decode_chunking(ebcc_lib, encoded, shape)
            assert np.isfinite(decoded).all()
            sizes.append(len(encoded))
            errors.append(np.mean((decoded.astype(np.float64) - data) ** 2))
    assert sizes[1] < sizes[0]
    assert errors[1] > errors[0]


@pytest.mark.parametrize('base_cr', [0, -1, float('nan'), float('inf'), -float('inf')])
@pytest.mark.parametrize('constant', [False, True])
def test_fixed_cr_mode_rejects_invalid_base_ratio(ebcc_lib, base_cr, constant):
    shape = (1, 32, 32)
    data = np.ones(shape, dtype=np.float32) if constant else make_data(shape)
    config = make_config(shape, base_cr=base_cr, residual_type=FIXED_CR)
    output = ctypes.POINTER(ctypes.c_ubyte)()
    assert ebcc_lib.ebcc_encode(data.ravel(), ctypes.byref(config), ctypes.byref(output)) == 0
    assert not output


@pytest.mark.parametrize('components', [1, 3])
def test_fixed_cr_mode_constant_field(ebcc_lib, components):
    data = np.full((components, 32, 32), 3.25, dtype=np.float32)
    encoded = encode_plain(ebcc_lib, data, make_config(data.shape, base_cr=30, residual_type=FIXED_CR, error=0))
    assert len(encoded) == ctypes.sizeof(EbccHeader) + 8
    np.testing.assert_array_equal(decode_chunking(ebcc_lib, encoded, data.shape), data)


@pytest.mark.parametrize('residual_type', [MAX_ERROR, RELATIVE_ERROR])
@pytest.mark.parametrize('base_cr', [0, -1, float('nan')])
def test_error_bound_modes_ignore_base_ratio(ebcc_lib, residual_type, base_cr):
    data = np.random.default_rng(19).normal(size=(3, 33, 65)).astype(np.float32)
    config = make_config(data.shape, base_cr=base_cr, residual_type=residual_type, error=.001)
    encoded = encode_plain(ebcc_lib, data, config)
    decoded = decode_chunking(ebcc_lib, encoded, data.shape)
    bound = config.error
    if residual_type == RELATIVE_ERROR:
        bound *= float(data.max()) - float(data.min())
    assert np.max(np.abs(decoded.astype(np.float64) - data)) <= bound


@pytest.mark.parametrize('args,base_cr,residual_opt', [
    (['--base_cr', '30'], 30, ('fixed_cr', 0)),
    (['-b', '10'], 10, ('fixed_cr', 0)),
    (['-m', '.5'], 200, ('max_error_target', .5)),
    (['-b', '0', '-m', '.5'], 0, ('max_error_target', .5)),
    (['--base_cr', '30', '--max_error_target', '.5'], 30, ('max_error_target', .5)),
    (['-r', '.02'], 200, ('relative_error_target', .02)),
    (['-b', '30', '-r', '.02'], 30, ('relative_error_target', .02)),
    (['--relative_error_target', '.02', '--base_cr', '30'], 30, ('relative_error_target', .02)),
    ([], 200, ('relative_error_target', .01)),
])
def test_cli_compression_mode(args, base_cr, residual_opt):
    from ebcc.filter_wrapper import EBCC_Filter
    script = Path(__file__).resolve().parents[1] / 'ebcc' / 'filter_wrapper.py'
    result = subprocess.run([sys.executable, str(script), '--height', '32', '--width', '32', *args],
                            capture_output=True, text=True, check=True)
    expected = (308,) + EBCC_Filter(base_cr, 32, 32, residual_opt).hdf_filter_opts
    assert tuple(map(int, result.stdout.strip().split(','))) == expected
    ignored_ratio = residual_opt[0] != 'fixed_cr' and any(arg in ('-b', '--base_cr') for arg in args)
    warning = 'Warning: --base_cr is ignored when an error bound is supplied.'
    assert (warning in result.stderr) == ignored_ratio
    assert warning not in result.stdout
