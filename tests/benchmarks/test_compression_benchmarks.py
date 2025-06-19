import pytest
import numpy as np
import xarray as xr
import h5py
import os
import time
from pathlib import Path
from ebcc.filter_wrapper import EBCC_Filter


class TestCompressionBenchmarks:
    """Benchmark tests for EBCC compression performance"""
    
    @pytest.mark.parametrize("error_target", [0.01, 0.1])
    def test_ebcc_compression_performance(self, geopotential_dataset, benchmark_output_dir, 
                                        error_target, benchmark):
        """Benchmark EBCC compression with different error targets"""
        data = geopotential_dataset['z'].data.astype(np.float32)
        
        def compress_ebcc():
            ebcc_filter = EBCC_Filter(
                base_cr=30,
                height=data.shape[-2],
                width=data.shape[-1],
                data_dim=len(data.shape),
                residual_opt=("max_error_target", error_target)
            )
            
            output_file = benchmark_output_dir / f"ebcc_err{error_target}.h5"
            with h5py.File(output_file, 'w') as f:
                f.create_dataset('z', shape=data.shape, **ebcc_filter)
                f['z'][...] = data
                compressed_data = f['z'][...]
            
            return compressed_data, output_file.stat().st_size
        
        result = benchmark(compress_ebcc)
        compressed_data, file_size = result
        
        # Verify compression worked
        assert compressed_data.shape == data.shape
        max_error = np.max(np.abs(data - compressed_data))
        assert max_error <= error_target * 1.5  # Allow some tolerance
        
        # Calculate compression ratio
        original_size = data.nbytes
        compression_ratio = original_size / file_size
        assert compression_ratio > 2  # Should achieve at least 2:1
    
    def test_ebcc_memory_usage(self, test_data_medium, benchmark_output_dir):
        """Test memory usage during compression"""
        import psutil
        import gc
        
        process = psutil.Process()
        gc.collect()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        ebcc_filter = EBCC_Filter(
            base_cr=30,
            height=test_data_medium.shape[0],
            width=test_data_medium.shape[1],
            data_dim=len(test_data_medium.shape),
            residual_opt=("relative_error_target", 0.01)
        )
        
        output_file = benchmark_output_dir / "memory_test.h5"
        with h5py.File(output_file, 'w') as f:
            f.create_dataset('data', shape=test_data_medium.shape, **ebcc_filter)
            f['data'][...] = test_data_medium
            
            peak_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        gc.collect()
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        memory_used = peak_memory - initial_memory
        memory_per_mb = memory_used / (test_data_medium.nbytes / 1024 / 1024)
        
        # Memory should be properly freed (allow some memory retention)
        assert abs(final_memory - initial_memory) < memory_used * 1.5
    
    @pytest.mark.parametrize("data_shape", [
        (64, 64),      # Small
        (256, 256),    # Medium  
        (512, 512),    # Large
    ])
    def test_ebcc_scalability(self, base_test_data, benchmark_output_dir, data_shape, benchmark):
        """Test EBCC performance scaling with data size"""
        # Crop or repeat data to target shape
        if data_shape[0] <= base_test_data.shape[0] and data_shape[1] <= base_test_data.shape[1]:
            test_data = base_test_data[:data_shape[0], :data_shape[1]].astype(np.float32)
        else:
            # Tile the data if needed
            repeats_y = (data_shape[0] + base_test_data.shape[0] - 1) // base_test_data.shape[0]
            repeats_x = (data_shape[1] + base_test_data.shape[1] - 1) // base_test_data.shape[1]
            tiled = np.tile(base_test_data, (repeats_y, repeats_x))
            test_data = tiled[:data_shape[0], :data_shape[1]].astype(np.float32)
        
        def compress_data():
            ebcc_filter = EBCC_Filter(
                base_cr=30,
                height=test_data.shape[0],
                width=test_data.shape[1],
                data_dim=len(test_data.shape),
                residual_opt=("relative_error_target", 0.01)
            )
            
            output_file = benchmark_output_dir / f"scale_{data_shape[0]}x{data_shape[1]}.h5"
            with h5py.File(output_file, 'w') as f:
                f.create_dataset('data', shape=test_data.shape, **ebcc_filter)
                f['data'][...] = test_data
            
            return output_file.stat().st_size
        
        file_size = benchmark(compress_data)
        
        # Performance should scale reasonably with data size
        data_mb = test_data.nbytes / 1024 / 1024
        throughput = data_mb / benchmark.stats['mean']  # MB/s
        
        # Should achieve at least 1 MB/s throughput
        assert throughput > 1.0
    
    def test_error_bound_accuracy(self, geopotential_dataset, benchmark_output_dir):
        """Test accuracy of error bound enforcement"""
        data = geopotential_dataset['z'].data.astype(np.float32)
        error_targets = [0.001, 0.01, 0.1, 1.0]
        
        results = []
        for target in error_targets:
            ebcc_filter = EBCC_Filter(
                base_cr=30,
                height=data.shape[-2],
                width=data.shape[-1],
                data_dim=len(data.shape),
                residual_opt=("max_error_target", target)
            )
            
            output_file = benchmark_output_dir / f"error_test_{target}.h5"
            with h5py.File(output_file, 'w') as f:
                f.create_dataset('z', shape=data.shape, **ebcc_filter)
                f['z'][...] = data
                compressed_data = f['z'][...]
            
            max_error = np.max(np.abs(data - compressed_data))
            results.append((target, max_error))
            
            # Error should not exceed target by more than 50%
            assert max_error <= target * 1.5
        
        # Verify error bounds are meaningful (smaller targets = smaller errors)
        for i in range(len(results) - 1):
            assert results[i][1] <= results[i+1][1] * 2  # Allow some variation