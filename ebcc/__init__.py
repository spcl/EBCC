import sys
from pathlib import Path


EBCC_FILTER_DIR = Path(__file__).resolve().parent
EBCC_FILTER_PATH = None
lib_name = None

if sys.platform == "darwin":
    _PATTERNS = ("libh5z_ebcc*.so", "libh5z_ebcc*.dylib")
elif sys.platform == "win32":
    _PATTERNS = ("h5z_ebcc*.dll",)
else:
    _PATTERNS = ("libh5z_ebcc*.so",)

for _pattern in _PATTERNS:
    matches = sorted(EBCC_FILTER_DIR.glob(_pattern))
    if matches:
        lib_name = matches[0].name
        EBCC_FILTER_PATH = str(matches[0])
        break

if EBCC_FILTER_PATH is None:
    raise FileNotFoundError(
        f"EBCC filter library not found in {EBCC_FILTER_DIR}. "
        f"Searched for: {', '.join(_PATTERNS)}. "
        "Please ensure the library is built and available in the expected directory."
    )

EBCC_FILTER_DIR = str(EBCC_FILTER_DIR)