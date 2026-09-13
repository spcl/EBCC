import os
from pathlib import Path

from ebcc import EBCC_FILTER_DIR

# Set HDF5 plugin path BEFORE importing h5py
os.environ["HDF5_PLUGIN_PATH"] = EBCC_FILTER_DIR

import pytest
import numpy as np
import h5py
import tempfile
import shutil

# Set up HDF5 plugin path
@pytest.fixture(scope="session", autouse=True)
def setup_hdf5_plugin():
    """Set up HDF5 plugin path for all tests"""
    # Store original value if it exists
    original_path = os.environ.get("HDF5_PLUGIN_PATH")
    
    # Set the plugin path
    os.environ["HDF5_PLUGIN_PATH"] = EBCC_FILTER_DIR
    
    yield
    
    # Restore original value
    if original_path is not None:
        os.environ["HDF5_PLUGIN_PATH"] = original_path
    else:
        os.environ.pop("HDF5_PLUGIN_PATH", None)

@pytest.fixture
def temp_dir():
    """Create temporary directory for test files"""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)

@pytest.fixture(scope="session")
def base_test_data():
    """Load base test data from file"""
    data_path = Path(__file__).parent.parent / "data" / "test_data.npy"
    return np.load(data_path)

@pytest.fixture
def test_data_small(base_test_data):
    """Small 2D test data - crop from base data"""
    return base_test_data[:64, :64].astype(np.float32)

@pytest.fixture
def test_data_medium(base_test_data):
    """Medium 2D test data - crop from base data"""
    return base_test_data[:256, :256].astype(np.float32)

@pytest.fixture
def test_data_3d(base_test_data):
    """3D test data - stack 2D slices"""
    slice_2d = base_test_data[:64, :64].astype(np.float32)
    return np.stack([slice_2d] * 10, axis=0)

@pytest.fixture
def test_data_4d(base_test_data):
    """4D test data - stack 3D volumes"""
    slice_2d = base_test_data[:64, :64].astype(np.float32)
    volume_3d = np.stack([slice_2d] * 5, axis=0)
    return np.stack([volume_3d] * 3, axis=0)

@pytest.fixture
def climate_data(base_test_data):
    """Full climate test data"""
    return base_test_data.astype(np.float32)

@pytest.fixture
def constant_data():
    """Constant data for edge case testing"""
    return np.full((64, 64), 42.0, dtype=np.float32)

@pytest.fixture
def zero_data():
    """Zero data for edge case testing"""
    return np.zeros((64, 64), dtype=np.float32)

@pytest.fixture
def temp_hdf5_file(temp_dir):
    """Create temporary HDF5 file"""
    file_path = temp_dir / "test.hdf5"
    yield file_path
    # Cleanup handled by temp_dir fixture

@pytest.fixture
def compression_ratios():
    """Standard compression ratios for testing"""
    return [10, 30, 100, 500]

@pytest.fixture
def error_thresholds():
    """Standard error thresholds for testing"""
    return {
        "max_error": [0.001, 0.01, 0.1],
        "relative_error": [0.001, 0.01, 0.1]
    }

@pytest.fixture
def residual_modes():
    """Supported residual compression modes"""
    return [
        ("max_error_target", 0.1),
        ("relative_error_target", 0.01)
    ]
