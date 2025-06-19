import pytest
import numpy as np
import xarray as xr
import os
from pathlib import Path


@pytest.fixture(scope="session")
def benchmark_data_dir():
    """Path to benchmark data directory"""
    return Path(__file__).parent.parent.parent / "data"


@pytest.fixture(scope="session")
def geopotential_dataset(benchmark_data_dir):
    """Load geopotential test dataset"""
    file_path = benchmark_data_dir / "geopotential_pl_small.nc"
    return xr.open_dataset(file_path)


@pytest.fixture(scope="session")
def era5_pl_dataset(benchmark_data_dir):
    """Load ERA5 pressure level dataset if available"""
    file_path = benchmark_data_dir / "era5_pl_sample.nc"
    if file_path.exists():
        return xr.open_dataset(file_path)
    return None


@pytest.fixture(scope="session")
def era5_sfc_dataset(benchmark_data_dir):
    """Load ERA5 surface dataset if available"""
    file_path = benchmark_data_dir / "era5_sfc_sample.nc"
    if file_path.exists():
        return xr.open_dataset(file_path)
    return None


@pytest.fixture
def compression_methods():
    """Available compression methods for comparison"""
    return ["ebcc", "sz", "sz3", "sperr"]


@pytest.fixture
def error_targets():
    """Standard error targets for benchmarking"""
    return [0.001, 0.01, 0.1]


@pytest.fixture
def base_compression_ratios():
    """Base compression ratios for EBCC testing"""
    return [10, 30, 50, 100, 200]


@pytest.fixture
def era5_pressure_levels():
    """ERA5 pressure levels for comprehensive testing"""
    return [1000, 975, 950, 925, 900, 875, 
            850, 825, 800, 775, 750, 700, 
            650, 600, 550, 500, 450, 400, 
            350, 300, 250, 225, 200, 175, 
            150, 125, 100, 70, 50, 30, 
            20, 10, 7, 5, 3, 2, 1]


@pytest.fixture
def benchmark_output_dir(tmp_path):
    """Temporary directory for benchmark outputs"""
    output_dir = tmp_path / "benchmark_outputs"
    output_dir.mkdir()
    return output_dir