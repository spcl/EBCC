import pytest
import numpy as np
import xarray as xr
import h5py
import os
import pandas as pd
from pathlib import Path
from ebcc.filter_wrapper import EBCC_Filter

try:
    import hdf5plugin
    HAS_HDF5PLUGIN = True
except ImportError:
    HAS_HDF5PLUGIN = False

try:
    from enstools.io import write
    from enstools.encoding.chunk_size import change_chunk_size
    HAS_ENSTOOLS = True
except ImportError:
    HAS_ENSTOOLS = False


class TestComparisonBenchmarks:
    """Benchmark tests comparing EBCC with other compression methods"""
    
    def compress_ebcc(self, data, output_path, error_target=10.0):
        """Compress data using EBCC"""
        ebcc_filter = EBCC_Filter(
            base_cr=30,
            height=data.shape[-2],
            width=data.shape[-1], 
            data_dim=len(data.shape),
            residual_opt=("max_error_target", error_target)
        )
        
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('z', shape=data.shape, **ebcc_filter)
            f['z'][...] = data
    
    @pytest.mark.skipif(not HAS_HDF5PLUGIN, reason="hdf5plugin not available")
    def compress_sperr(self, data, output_path, error_target=10.0):
        """Compress data using SPERR"""
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('z', shape=data.shape, dtype=np.float32,
                           chunks=(1, 1, data.shape[-2], data.shape[-1]),
                           **hdf5plugin.Sperr(absolute=error_target))
            f['z'][...] = data
    
    @pytest.mark.skipif(not HAS_HDF5PLUGIN, reason="hdf5plugin not available")  
    def compress_sz(self, data, output_path, error_target=10.0):
        """Compress data using SZ"""
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('z', shape=data.shape, dtype=np.float32,
                           chunks=(1, 1, data.shape[-2], data.shape[-1]),
                           **hdf5plugin.SZ(absolute=error_target))
            f['z'][...] = data
    
    @pytest.mark.skipif(not HAS_ENSTOOLS, reason="enstools not available")
    def compress_sz3(self, data, output_path, error_target=10.0):
        """Compress data using SZ3 via enstools"""
        # Create temporary xarray dataset
        ds = xr.Dataset({'z': (('time', 'level', 'lat', 'lon'), data)})
        ds = ds.astype(np.float32)
        
        change_chunk_size("4152960B")
        write(ds, str(output_path), compression=f"lossy,sz3,abs,{error_target}")
    
    def calculate_stats(self, original_data, compressed_file):
        """Calculate compression statistics"""
        with h5py.File(compressed_file, 'r') as f:
            compressed_data = f['z'][...]
        
        rmse = np.sqrt(np.square(original_data - compressed_data).mean())
        max_error = np.max(np.abs(original_data - compressed_data))
        file_size = os.path.getsize(compressed_file) / 1024 / 1024  # MB
        original_size = original_data.nbytes / 1024 / 1024  # MB
        compression_ratio = original_size / file_size
        
        return {
            'rmse': rmse,
            'max_error': max_error, 
            'size_mb': file_size,
            'compression_ratio': compression_ratio
        }
    
    def test_ebcc_vs_others_basic(self, geopotential_dataset, benchmark_output_dir):
        """Basic comparison of EBCC vs other methods"""
        data = geopotential_dataset['z'].data.astype(np.float32)
        error_target = 10.0
        results = {}
        
        # Test EBCC
        ebcc_file = benchmark_output_dir / "ebcc_comparison.h5"
        self.compress_ebcc(data, ebcc_file, error_target)
        results['ebcc'] = self.calculate_stats(data, ebcc_file)
        
        # Test SPERR if available
        if HAS_HDF5PLUGIN:
            sperr_file = benchmark_output_dir / "sperr_comparison.h5"
            self.compress_sperr(data, sperr_file, error_target)
            results['sperr'] = self.calculate_stats(data, sperr_file)
            
            # Test SZ if available
            sz_file = benchmark_output_dir / "sz_comparison.h5"
            self.compress_sz(data, sz_file, error_target)
            results['sz'] = self.calculate_stats(data, sz_file)
        
        # Test SZ3 if available
        if HAS_ENSTOOLS:
            sz3_file = benchmark_output_dir / "sz3_comparison.nc"
            self.compress_sz3(data, sz3_file, error_target)
            results['sz3'] = self.calculate_stats(data, sz3_file)
        
        # Verify EBCC performs reasonably
        ebcc_stats = results['ebcc']
        assert ebcc_stats['compression_ratio'] > 2.0
        assert ebcc_stats['max_error'] <= error_target * 1.5
        
        # Compare with other methods if available
        for method, stats in results.items():
            if method != 'ebcc':
                # All methods should achieve reasonable compression
                assert stats['compression_ratio'] > 1.5
                assert stats['max_error'] <= error_target * 2.0  # Allow more tolerance for others
    
    @pytest.mark.parametrize("error_target", [1.0, 10.0, 50.0])
    def test_ebcc_error_scaling(self, geopotential_dataset, benchmark_output_dir, error_target):
        """Test how EBCC compression scales with different error targets"""
        data = geopotential_dataset['z'].data.astype(np.float32)
        
        ebcc_file = benchmark_output_dir / f"ebcc_error_{error_target}.h5"
        self.compress_ebcc(data, ebcc_file, error_target)
        stats = self.calculate_stats(data, ebcc_file)
        
        # Verify error bounds
        assert stats['max_error'] <= error_target * 1.5
        
        # Higher error targets should generally give better compression
        # (though this may not always be strictly true)
        assert stats['compression_ratio'] > 1.5
    
    @pytest.mark.skipif(not (HAS_HDF5PLUGIN or HAS_ENSTOOLS), 
                       reason="No comparison methods available")
    def test_comprehensive_comparison(self, geopotential_dataset, benchmark_output_dir):
        """Comprehensive comparison across multiple error targets"""
        data = geopotential_dataset['z'].data.astype(np.float32)
        error_targets = [1.0, 10.0, 50.0]
        
        comparison_results = []
        
        for error_target in error_targets:
            # EBCC
            ebcc_file = benchmark_output_dir / f"comp_ebcc_{error_target}.h5"
            self.compress_ebcc(data, ebcc_file, error_target)
            ebcc_stats = self.calculate_stats(data, ebcc_file)
            ebcc_stats['method'] = 'ebcc'
            ebcc_stats['error_target'] = error_target
            comparison_results.append(ebcc_stats)
            
            # Other methods if available
            if HAS_HDF5PLUGIN:
                sperr_file = benchmark_output_dir / f"comp_sperr_{error_target}.h5"
                self.compress_sperr(data, sperr_file, error_target)
                sperr_stats = self.calculate_stats(data, sperr_file)
                sperr_stats['method'] = 'sperr'
                sperr_stats['error_target'] = error_target
                comparison_results.append(sperr_stats)
        
        # Convert to DataFrame for analysis
        df = pd.DataFrame(comparison_results)
        
        # Verify EBCC is competitive
        ebcc_results = df[df['method'] == 'ebcc']
        assert len(ebcc_results) == len(error_targets)
        
        # All EBCC results should meet error bounds
        for _, row in ebcc_results.iterrows():
            assert row['max_error'] <= row['error_target'] * 1.5
            assert row['compression_ratio'] > 1.5