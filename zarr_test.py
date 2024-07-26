import numpy as np
from src.zarr_filter import J2KFilter

data = np.ones((10000,), dtype=np.float32)

filter = J2KFilter((100, 100, 1128792064, 0))

encoded = filter.encode(data.tobytes())

decoded = np.frombuffer(filter.decode(encoded), dtype=np.float32)
