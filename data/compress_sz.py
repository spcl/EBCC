from enstools.io import write
from enstools.encoding.chunk_size import change_chunk_size
from enstools.compression.analyzer.analyzer import analyze_dataset
import xarray
import numpy as np

#change_chunk_size("1024KB")
change_chunk_size("4152960B")
dataset = xarray.open_dataset(f'geopotential_pl_small.nc').astype(np.float32)
#encoding, metrics = analyze_dataset(dataset,
#                                    constrains="absolute:10.0",)
#write(dataset, "geopotential_pl_small_sz.nc", compression="lossy,sz3,abs,10.0")
write(dataset, "geopotential_pl_small_sz.nc", compression="lossy,sz,abs,10.0")
