import numpy as np
import zarr_filter

data = np.ones((10000,), dtype=np.float32)

filter = zarr_filter.J2KFilter((100, 100, 1128792064, 0))

encoded = filter.encode(data.tobytes())

decoded = np.frombuffer(filter.decode(encoded), dtype=np.float32)
