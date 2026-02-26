import tempfile

import h5py
import numpy as np
import pytest

from ebcc.filter_wrapper import EBCC_Filter


_JXL_HDF5_AVAILABLE = None


def _base_kwargs(compressor: str, base_param: float):
    if compressor == "jxl":
        return {"base_compressor": "jxl", "base_distance": float(base_param)}
    return {"base_compressor": "j2k", "base_cr": float(base_param)}


def _is_jxl_hdf5_available() -> bool:
    global _JXL_HDF5_AVAILABLE
    if _JXL_HDF5_AVAILABLE is not None:
        return _JXL_HDF5_AVAILABLE

    try:
        data = np.arange(32 * 32, dtype=np.float32).reshape(32, 32)
        with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
            ebcc_filter = EBCC_Filter(
                **_base_kwargs("jxl", 1.0),
                height=32,
                width=32,
                data_dim=2,
                residual_opt=("none", 0.0),
            )
            with h5py.File(tmp.name, "w") as f:
                f.create_dataset("d", shape=data.shape, **ebcc_filter)
                f["d"][:] = data
                _ = f["d"][:]
        _JXL_HDF5_AVAILABLE = True
    except Exception:
        _JXL_HDF5_AVAILABLE = False

    return _JXL_HDF5_AVAILABLE


class TestNetCDFIntegration:
    """Test NetCDF/HDF5 integration with EBCC filter"""

    def _maybe_skip_jxl(self, compressor: str):
        if compressor == "jxl" and not _is_jxl_hdf5_available():
            pytest.skip("JXL HDF5 filter is not available in this build")

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_basic_compression_decompression(self, climate_data, temp_hdf5_file, compressor):
        self._maybe_skip_jxl(compressor)

        base_param = 100.0 if compressor == "j2k" else 1.0
        ebcc_filter = EBCC_Filter(
            **_base_kwargs(compressor, base_param),
            height=climate_data.shape[0],
            width=climate_data.shape[1],
            data_dim=len(climate_data.shape),
            residual_opt=("relative_error_target", 0.009),
        )

        with h5py.File(temp_hdf5_file, 'w') as f:
            f.create_dataset('compressed', shape=climate_data.shape, **ebcc_filter)
            f['compressed'][:] = climate_data
            uncompressed = f['compressed'][:]

        assert uncompressed.shape == climate_data.shape
        assert uncompressed.dtype == climate_data.dtype

        data_range = np.max(climate_data) - np.min(climate_data)
        max_error = np.max(np.abs(climate_data - uncompressed))
        if data_range > 0:
            rel_error = max_error / data_range
            assert rel_error <= 0.02

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_compression_ratio(self, climate_data, temp_hdf5_file, temp_dir, compressor):
        self._maybe_skip_jxl(compressor)

        base_param = 100.0 if compressor == "j2k" else 1.0
        ebcc_filter = EBCC_Filter(
            **_base_kwargs(compressor, base_param),
            height=climate_data.shape[0],
            width=climate_data.shape[1],
            data_dim=len(climate_data.shape),
            residual_opt=("relative_error_target", 0.009),
        )

        original_file = temp_dir / "original.npy"
        np.save(original_file, climate_data)

        with h5py.File(temp_hdf5_file, 'w') as f:
            f.create_dataset('compressed', shape=climate_data.shape, **ebcc_filter)
            f['compressed'][:] = climate_data

        original_size = original_file.stat().st_size
        compressed_size = temp_hdf5_file.stat().st_size
        compression_ratio = original_size / compressed_size

        assert compression_ratio > 2

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_different_compression_ratios(self, test_data_medium, temp_hdf5_file, compressor):
        self._maybe_skip_jxl(compressor)

        base_params = [10, 50, 100, 200] if compressor == "j2k" else [0.5, 1.0, 2.0, 4.0]
        for base_param in base_params:
            ebcc_filter = EBCC_Filter(
                **_base_kwargs(compressor, base_param),
                height=test_data_medium.shape[0],
                width=test_data_medium.shape[1],
                data_dim=len(test_data_medium.shape),
                residual_opt=("relative_error_target", 0.01),
            )

            with h5py.File(temp_hdf5_file, 'w') as f:
                f.create_dataset('compressed', shape=test_data_medium.shape, **ebcc_filter)
                f['compressed'][:] = test_data_medium
                uncompressed = f['compressed'][:]

            assert uncompressed.shape == test_data_medium.shape

    @pytest.mark.parametrize("compressor", ["j2k", "jxl"])
    def test_multidimensional_data(self, test_data_3d, temp_hdf5_file, compressor):
        self._maybe_skip_jxl(compressor)

        base_param = 50.0 if compressor == "j2k" else 1.0
        ebcc_filter = EBCC_Filter(
            **_base_kwargs(compressor, base_param),
            height=test_data_3d.shape[1],
            width=test_data_3d.shape[2],
            data_dim=len(test_data_3d.shape),
            residual_opt=("max_error_target", 0.1),
        )

        with h5py.File(temp_hdf5_file, 'w') as f:
            f.create_dataset('compressed', shape=test_data_3d.shape, **ebcc_filter)
            f['compressed'][:] = test_data_3d
            uncompressed = f['compressed'][:]

        assert uncompressed.shape == test_data_3d.shape
        max_error = np.max(np.abs(test_data_3d - uncompressed))
        assert max_error <= 0.15
