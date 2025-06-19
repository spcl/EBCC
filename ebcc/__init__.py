import os

EBCC_FILTER_DIR = os.path.join(os.path.dirname(__file__), "../src/build/lib")
EBCC_FILTER_PATH = os.path.join(EBCC_FILTER_DIR, "libh5z_j2k.so")

if not os.path.exists(EBCC_FILTER_PATH):
    raise FileNotFoundError(
        f"EBCC filter library not found at {EBCC_FILTER_PATH}. "
        "Please ensure the library is built and available in the expected directory."
    )