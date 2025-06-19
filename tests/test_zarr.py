import pytest
import numpy as np
import tempfile
import shutil
from pathlib import Path

try:
    import zarr
    from zarr_any_numcodecs import AnyNumcodecsArrayArrayCodec, AnyNumcodecsArrayBytesCodec, AnyNumcodecsBytesBytesCodec
    from ebcc.zarr_filter import EBCCZarrFilter
    from ebcc.filter_wrapper import EBCC_Filter
    HAS_ZARR = True
except ImportError:
    HAS_ZARR = False

pytestmark = pytest.mark.skipif(not HAS_ZARR, reason="Zarr dependencies not available")


class TestZarrIntegration:
    """Test Zarr integration with EBCC filter"""
    
    def create_test_data(self, height=32, width=32):
        """Create test data with gradient pattern"""
        data = np.zeros((1, 1, height, width), dtype=np.float32)
        for i in range(height):
            for j in range(width):
                data[0, 0, i, j] = i + j
        return data
    
    def test_basic_zarr_filter_encode_decode(self):
        """Test basic EBCC Zarr filter encode/decode functionality"""
        width = 32
        height = 32
        data = self.create_test_data(height, width)
        atol = 1e-2
        
        # Create EBCC filter configuration
        ebcc_filter = EBCC_Filter(
            base_cr=2, 
            height=height, 
            width=width, 
            residual_opt=("max_error_target", atol)
        )
        
        # Create Zarr filter
        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)
        
        # Test encode
        encoded = zarr_filter.encode(data)
        assert encoded is not None
        assert len(encoded) > 0
        
        # Test decode
        decoded = zarr_filter.decode(encoded)
        decoded = decoded.reshape(data.shape)
        
        # Verify data integrity within tolerance
        assert decoded.shape == data.shape
        assert decoded.dtype == data.dtype
        assert np.allclose(data, decoded, atol=atol), "Decoded data does not match original within tolerance"
    
    def test_zarr_array_storage(self, tmp_path):
        """Test storing and retrieving data from Zarr array with EBCC filter"""
        width = 32
        height = 32
        data = self.create_test_data(height, width)
        atol = 1e-2
        
        # Create EBCC filter configuration
        ebcc_filter = EBCC_Filter(
            base_cr=2, 
            height=height, 
            width=width, 
            residual_opt=("max_error_target", atol)
        )
        
        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)
        
        # Create Zarr array with EBCC compression
        zarr_store = tmp_path / "test.zarr"
        z = zarr.create_array(
            store=str(zarr_store),
            data=data,
            chunks="auto",
            serializer=AnyNumcodecsArrayBytesCodec(zarr_filter),
            compressors=None,
            overwrite=True,
        )
        
        # Read back data
        decoded_zarr = z[:]
        
        # Verify data integrity
        assert decoded_zarr.shape == data.shape
        assert np.allclose(data, decoded_zarr, atol=atol), "Decoded Zarr data does not match original within tolerance"
    
    @pytest.mark.parametrize("base_cr", [2, 5, 10])
    def test_different_compression_ratios(self, base_cr):
        """Test Zarr filter with different compression ratios"""
        width = 32
        height = 32
        data = self.create_test_data(height, width)
        atol = 1e-2
        
        ebcc_filter = EBCC_Filter(
            base_cr=base_cr, 
            height=height, 
            width=width, 
            residual_opt=("max_error_target", atol)
        )
        
        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)
        
        # Test compression and decompression
        encoded = zarr_filter.encode(data)
        decoded = zarr_filter.decode(encoded).reshape(data.shape)
        
        # Verify data integrity
        assert np.allclose(data, decoded, atol=atol)
        
        # Higher compression ratios should generally produce smaller output
        original_size = data.nbytes
        compressed_size = len(encoded)
        compression_ratio = original_size / compressed_size
        
        # Should achieve some compression
        assert compression_ratio > 1.0, f"Poor compression ratio: {compression_ratio}"
    
    @pytest.mark.parametrize("max_error_target", [1e-2, 1e-1])
    def test_different_error_targets(self, max_error_target):
        """Test Zarr filter with different error targets"""
        width = 32
        height = 32
        data = self.create_test_data(height, width)
        
        ebcc_filter = EBCC_Filter(
            base_cr=2, 
            height=height, 
            width=width, 
            residual_opt=("max_error_target", max_error_target)
        )
        
        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)
        
        # Test compression and decompression
        encoded = zarr_filter.encode(data)
        decoded = zarr_filter.decode(encoded).reshape(data.shape)
        
        # Verify data integrity within specified tolerance
        max_error = np.max(np.abs(data - decoded))
        assert max_error <= max_error_target * 1.5, f"Max error {max_error} exceeds target {max_error_target}"

    def test_zarr_filter_different_data_sizes(self):
        """Test Zarr filter with different data dimensions"""
        test_cases = [
            (32, 32), 
            (64, 32),
            (32, 64)
        ]
        
        atol = 1e-2
        
        for height, width in test_cases:
            data = self.create_test_data(height, width)
            
            ebcc_filter = EBCC_Filter(
                base_cr=2, 
                height=height, 
                width=width, 
                residual_opt=("max_error_target", atol)
            )
            
            zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)
            
            # Test compression and decompression
            encoded = zarr_filter.encode(data)
            decoded = zarr_filter.decode(encoded).reshape(data.shape)
            
            # Verify data integrity
            assert decoded.shape == data.shape
            assert np.allclose(data, decoded, atol=atol)
    
    def test_zarr_filter_constant_data(self):
        """Test Zarr filter with constant data (should compress well)"""
        width = 32
        height = 32
        data = np.full((1, 1, height, width), 42.0, dtype=np.float32)
        atol = 1e-6  # Constant data should be very accurate
        
        ebcc_filter = EBCC_Filter(
            base_cr=2, 
            height=height, 
            width=width, 
            residual_opt=("max_error_target", atol)
        )
        
        zarr_filter = EBCCZarrFilter(ebcc_filter.hdf_filter_opts)
        
        # Test compression
        encoded = zarr_filter.encode(data)
        original_size = data.nbytes
        compressed_size = len(encoded)
        compression_ratio = original_size / compressed_size
        
        # Constant data should compress very well
        assert compression_ratio > 5.0, f"Poor compression for constant data: {compression_ratio}"
        
        # Test decompression
        decoded = zarr_filter.decode(encoded).reshape(data.shape)
        assert np.allclose(data, decoded, atol=atol)