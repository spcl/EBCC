from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path


def _candidate_lib_names():
    """
    Generate possible library filenames using Python's extension suffixes first,
    then fall back to legacy names for pre-existing builds.
    """
    seen = set()
    for suffix in EXTENSION_SUFFIXES:
        name = f"libh5z_ebcc{suffix}"
        if name not in seen:
            seen.add(name)
            yield name
    for name in ("libh5z_ebcc.so", "libh5z_ebcc.dylib", "h5z_ebcc.dll"):
        if name not in seen:
            yield name


EBCC_FILTER_DIR = Path(__file__).resolve().parent
EBCC_FILTER_PATH = None
lib_name = None

for _name in _candidate_lib_names():
    candidate = EBCC_FILTER_DIR / _name
    if candidate.exists():
        lib_name = _name
        EBCC_FILTER_PATH = str(candidate)
        break

if EBCC_FILTER_PATH is None:
    candidates = ", ".join(_candidate_lib_names())
    raise FileNotFoundError(
        f"EBCC filter library not found in {EBCC_FILTER_DIR}. "
        f"Tried: {candidates}. "
        "Please ensure the library is built and available in the expected directory."
    )
