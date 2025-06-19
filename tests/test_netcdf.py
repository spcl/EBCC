import pytest
import h5py
import numpy as np
from pathlib import Path
from ebcc.filter_wrapper import EBCC_Filter


class TestNetCDFIntegration:
    """Test NetCDF/HDF5 integration with EBCC filter"""
    
    def test_basic_compression_decompression(self, climate_data, temp_hdf5_file):
        """Test basic compression and decompression round-trip"""
        ebcc_filter = EBCC_Filter(
            base_cr=100,
            height=climate_data.shape[0],
            width=climate_data.shape[1],
            data_dim=len(climate_data.shape),
            residual_opt=("relative_error_target", 0.009)
        )
        
        with h5py.File(temp_hdf5_file, 'w') as f:
            f.create_dataset('compressed', shape=climate_data.shape, **ebcc_filter)
            f['compressed'][:] = climate_data
            uncompressed = f['compressed'][:]
        
        # Verify data integrity
        assert uncompressed.shape == climate_data.shape
        assert uncompressed.dtype == climate_data.dtype
        
        # Check error bounds
        data_range = np.max(climate_data) - np.min(climate_data)
        max_error = np.max(np.abs(climate_data - uncompressed))
        if data_range > 0:
            rel_error = max_error / data_range
            assert rel_error <= 0.02  # Allow some tolerance above target
        
    def test_compression_ratio(self, climate_data, temp_hdf5_file, temp_dir):
        """Test that compression achieves reasonable ratios"""
        ebcc_filter = EBCC_Filter(
            base_cr=100,
            height=climate_data.shape[0],
            width=climate_data.shape[1],
            data_dim=len(climate_data.shape),
            residual_opt=("relative_error_target", 0.009)
        )
        
        # Save original data
        original_file = temp_dir / "original.npy"
        np.save(original_file, climate_data)
        
        # Compress data
        with h5py.File(temp_hdf5_file, 'w') as f:
            f.create_dataset('compressed', shape=climate_data.shape, **ebcc_filter)
            f['compressed'][:] = climate_data
        
        # Check compression ratio
        original_size = original_file.stat().st_size
        compressed_size = temp_hdf5_file.stat().st_size
        compression_ratio = original_size / compressed_size
        
        assert compression_ratio > 5  # Should achieve at least 5:1 compression
        
    @pytest.mark.parametrize("base_cr", [10, 50, 100, 200])
    def test_different_compression_ratios(self, test_data_medium, temp_hdf5_file, base_cr):
        """Test different base compression ratios"""
        ebcc_filter = EBCC_Filter(
            base_cr=base_cr,
            height=test_data_medium.shape[0],
            width=test_data_medium.shape[1],
            data_dim=len(test_data_medium.shape),
            residual_opt=("relative_error_target", 0.01)
        )
        
        with h5py.File(temp_hdf5_file, 'w') as f:
            f.create_dataset('compressed', shape=test_data_medium.shape, **ebcc_filter)
            f['compressed'][:] = test_data_medium
            uncompressed = f['compressed'][:]
        
        # Higher compression ratios should still work
        assert uncompressed.shape == test_data_medium.shape
        
    def test_multidimensional_data(self, test_data_3d, temp_hdf5_file):
        """Test compression of 3D data"""
        ebcc_filter = EBCC_Filter(
            base_cr=50,
            height=test_data_3d.shape[1],
            width=test_data_3d.shape[2],
            data_dim=len(test_data_3d.shape),
            residual_opt=("max_error_target", 0.1)
        )
        
        with h5py.File(temp_hdf5_file, 'w') as f:
            f.create_dataset('compressed', shape=test_data_3d.shape, **ebcc_filter)
            f['compressed'][:] = test_data_3d
            uncompressed = f['compressed'][:]
        
        assert uncompressed.shape == test_data_3d.shape
        max_error = np.max(np.abs(test_data_3d - uncompressed))
        assert max_error <= 0.15  # Allow some tolerance
