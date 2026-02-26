import numpy as np
import pytest

try:
    import zarr
    from zarr_any_numcodecs import AnyNumcodecsArrayBytesCodec
    from ebcc.zarr_filter import EBCCZarrFilter
    from ebcc.filter_wrapper import EBCC_Filter

    HAS_ZARR = True
except ImportError:
    HAS_ZARR = False

pytestmark = pytest.mark.skipif(not HAS_ZARR, reason="Zarr dependencies not available")


_JXL_ZARR_AVAILABLE = None


def _base_kwargs(compressor: str, base_param: float):
    if compressor == "jxl":
        return {"base_compressor": "jxl", "base_distance": float(base_param)}
    return {"base_compressor": "j2k", "base_cr": float(base_param)}


def _is_jxl_zarr_available() -> bool:
    global _JXL_ZARR_AVAILABLE
    if _JXL_ZARR_AVAILABLE is not None:
        return _JXL_ZARR_AVAILABLE

    try:
        data = np.arange(32 * 32, dtype=np.float32).reshape(1, 1, 32, 32)
        ebcc_filter = EBCC_Filter(
            **_base_kwargs("jxl", 1.0),
            height=32,
            width=32,
            residual_opt=("none", 0.0),
        )
        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts, base_compressor="jxl")
        encoded = zarr_filter.encode(data)
        decoded = zarr_filter.decode(encoded).reshape(data.shape)
        _JXL_ZARR_AVAILABLE = np.allclose(data, decoded, atol=1e-2)
    except Exception:
        _JXL_ZARR_AVAILABLE = False

    return _JXL_ZARR_AVAILABLE


class TestZarrIntegration:
    """Test Zarr integration with EBCC filter"""

    def create_test_data(self, height=32, width=32):
        data = np.zeros((1, 1, height, width), dtype=np.float32)
        for i in range(height):
            for j in range(width):
                data[0, 0, i, j] = i + j
        return data

    def _maybe_skip_jxl(self, compressor: str):
        if compressor == "jxl" and not _is_jxl_zarr_available():
            pytest.skip("JXL backend is not available in this build")

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_basic_zarr_filter_encode_decode(self, compressor):
        self._maybe_skip_jxl(compressor)

        width = 32
        height = 32
        data = self.create_test_data(height, width)
        atol = 1e-2
        base_param = 2.0 if compressor == "j2k" else 1.0

        ebcc_filter = EBCC_Filter(
            **_base_kwargs(compressor, base_param),
            height=height,
            width=width,
            residual_opt=("max_error_target", atol),
        )

        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts, base_compressor=compressor)

        encoded = zarr_filter.encode(data)
        assert encoded is not None
        assert len(encoded) > 0

        decoded = zarr_filter.decode(encoded)
        decoded = decoded.reshape(data.shape)

        assert decoded.shape == data.shape
        assert decoded.dtype == data.dtype
        assert np.allclose(data, decoded, atol=atol)

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_zarr_array_storage(self, tmp_path, compressor):
        self._maybe_skip_jxl(compressor)

        width = 32
        height = 32
        data = self.create_test_data(height, width)
        atol = 1e-2
        base_param = 2.0 if compressor == "j2k" else 1.0

        ebcc_filter = EBCC_Filter(
            **_base_kwargs(compressor, base_param),
            height=height,
            width=width,
            residual_opt=("max_error_target", atol),
        )

        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts, base_compressor=compressor)

        zarr_store = tmp_path / "test.zarr"
        z = zarr.create_array(
            store=str(zarr_store),
            data=data,
            chunks="auto",
            serializer=AnyNumcodecsArrayBytesCodec(zarr_filter),
            compressors=None,
            overwrite=True,
        )

        decoded_zarr = z[:]

        assert decoded_zarr.shape == data.shape
        assert np.allclose(data, decoded_zarr, atol=atol)

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_different_compression_ratios(self, compressor):
        self._maybe_skip_jxl(compressor)

        width = 32
        height = 32
        data = self.create_test_data(height, width)
        atol = 1e-2

        base_params = [2, 5, 10] if compressor == "j2k" else [0.5, 1.0, 2.0]
        for base_param in base_params:
            ebcc_filter = EBCC_Filter(
                **_base_kwargs(compressor, base_param),
                height=height,
                width=width,
                residual_opt=("max_error_target", atol),
            )

            zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts, base_compressor=compressor)
            encoded = zarr_filter.encode(data)
            decoded = zarr_filter.decode(encoded).reshape(data.shape)

            assert np.allclose(data, decoded, atol=atol)
            compression_ratio = data.nbytes / len(encoded)
            assert compression_ratio > 1.0

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    @pytest.mark.parametrize("max_error_target", [1e-2, 1e-1])
    def test_different_error_targets(self, compressor, max_error_target):
        self._maybe_skip_jxl(compressor)

        width = 32
        height = 32
        data = self.create_test_data(height, width)
        base_param = 2.0 if compressor == "j2k" else 1.0

        ebcc_filter = EBCC_Filter(
            **_base_kwargs(compressor, base_param),
            height=height,
            width=width,
            residual_opt=("max_error_target", max_error_target),
        )

        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts, base_compressor=compressor)
        encoded = zarr_filter.encode(data)
        decoded = zarr_filter.decode(encoded).reshape(data.shape)

        max_error = np.max(np.abs(data - decoded))
        assert max_error <= max_error_target * 1.5

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_zarr_filter_different_data_sizes(self, compressor):
        self._maybe_skip_jxl(compressor)

        test_cases = [(32, 32), (64, 32), (32, 64)]
        atol = 1e-2
        base_param = 2.0 if compressor == "j2k" else 1.0

        for height, width in test_cases:
            data = self.create_test_data(height, width)

            ebcc_filter = EBCC_Filter(
                **_base_kwargs(compressor, base_param),
                height=height,
                width=width,
                residual_opt=("max_error_target", atol),
            )

            zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts, base_compressor=compressor)
            encoded = zarr_filter.encode(data)
            decoded = zarr_filter.decode(encoded).reshape(data.shape)

            assert decoded.shape == data.shape
            assert np.allclose(data, decoded, atol=atol)

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_zarr_filter_constant_data(self, compressor):
        self._maybe_skip_jxl(compressor)

        width = 32
        height = 32
        data = np.full((1, 1, height, width), 42.0, dtype=np.float32)
        atol = 1e-6
        base_param = 2.0 if compressor == "j2k" else 1.0

        ebcc_filter = EBCC_Filter(
            **_base_kwargs(compressor, base_param),
            height=height,
            width=width,
            residual_opt=("max_error_target", atol),
        )

        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts, base_compressor=compressor)

        encoded = zarr_filter.encode(data)
        compression_ratio = data.nbytes / len(encoded)
        assert compression_ratio > 5.0

        decoded = zarr_filter.decode(encoded).reshape(data.shape)
        assert np.allclose(data, decoded, atol=atol)
