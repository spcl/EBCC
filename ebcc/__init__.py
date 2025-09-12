import os
import sys

if sys.platform.startswith('linux'):
    lib_name = 'libh5z_ebcc.so'
elif sys.platform == 'darwin':
    lib_name = 'libh5z_ebcc.dylib'
#elif sys.platform == 'win32':
#    lib_name = 'h5z_ebcc.dll'
else:
    raise RuntimeError(f"Unsupported platform: {sys.platform}")


EBCC_FILTER_DIR = os.path.join(os.path.dirname(__file__), "../src/build/lib")
EBCC_FILTER_PATH = os.path.join(EBCC_FILTER_DIR, lib_name)

if not os.path.exists(EBCC_FILTER_PATH):
    raise FileNotFoundError(
        f"EBCC filter library not found at {EBCC_FILTER_PATH}. "
        "Please ensure the library is built and available in the expected directory."
    )