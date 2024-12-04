import os
base_folder = "/home/huanglangwen/Documents/compression-filter"
os.environ["HDF5_PLUGIN_PATH"] = os.path.join(base_folder, 'src/build/lib')
import sys
sys.path.append(base_folder)
import xarray as xr
import h5py
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from filter_wrapper import JP2SPWV_Filter
from tqdm import tqdm
import hdf5plugin

METHODS = ["sz", "jp2spwv", "sz3", "sperr"]
LEVELS = [1000, 975, 950, 925, 900, 875, 
          850, 825, 800, 775, 750, 700, 
          650, 600, 550, 500, 450, 400, 
          350, 300, 250, 225, 200, 175, 
          150, 125, 100, 70, 50, 30, 
          20, 10, 7, 5, 3, 2, 1]


def prepare_encoding(method, abs_target, height, width, data_dim):
    #os.environ["HDF5_PLUGIN_PATH"] = hdf5plugin.PLUGIN_PATH
    if method == "jp2spwv":
        filter = JP2SPWV_Filter(
            base_cr=30, height=height, width=width,
            data_dim=data_dim, residual_opt=("max_error_target", abs_target),
            filter_path=os.path.join(base_folder, 'src/build/lib'))
    elif method == "sz":
        filter = hdf5plugin.SZ(absolute=abs_target)
    elif method == "sz3":
        filter = hdf5plugin.SZ3(absolute=abs_target)
    elif method == "sperr":
        filter = hdf5plugin.Sperr(absolute=abs_target)
    else:
        raise ValueError(f"Unknown method {method}")
    return dict(filter)

def compress(ds, name, output_name, encoding):
    nlon = ds.sizes["longitude"]
    nlat = ds.sizes["latitude"]
    if ds.ndim == 3:
        chunksizes = (1, nlat, nlon)
    elif ds.ndim == 4:
        chunksizes = (1, 1, nlat, nlon)
    else:
        raise ValueError("Only 3D and 4D data supported")
    chunksizes = encoding.pop("chunks") if "chunks" in encoding else chunksizes
    encoding.pop("dtype", None)
    encoding_nc = {name: {**encoding, "chunksizes": chunksizes}}
    ds = ds.astype(np.float32)
    ds.to_netcdf(output_name, encoding=encoding_nc, engine="h5netcdf")

def compare(ds_orig, filename):
    ds_comp = xr.open_dataarray(filename, engine="h5netcdf")
    delta = ds_comp - ds_orig
    rmse = np.sqrt(np.square(delta).mean()).item()
    max_diff = np.max(np.abs(delta)).item()
    size = os.path.getsize(filename) / 1024 / 1024
    original_size = ds_orig.nbytes / 1024 / 1024
    cr = original_size / size
    return rmse, max_diff, size, cr

def compress_and_compare(ds, var, level, output_folder, method, rel_target):
    height = ds['latitude'].size
    width = ds['longitude'].size
    da = ds[var]
    output_name = f"era5_{var}_{method}.nc"
    var_name = var
    if da.ndim == 4:
        da = da.sel(pressure_level=level)
        output_name = f"era5_{var}_{level}_{method}.nc"
        var_name = f"{var}_{level}"
    data_dim = len(da.shape)
    stats = []
    data_range = (da.max(dim=["latitude", "longitude"]) - da.min(dim=["latitude", "longitude"])).mean().item()
    abs_target = data_range * rel_target
    print(f"Compressing {var} with {method}, relative error {rel_target}, absolute error {abs_target}", flush=True)
    encoding = prepare_encoding(method, abs_target, height, width, data_dim)
    output_path = os.path.join(output_folder, output_name)
    compress(da, var, output_path, encoding)
    rmse, max_diff, size, cr = compare(ds[var], output_path)
    stats.append([method, var_name, rmse, max_diff, size, cr])
    stats_df = pd.DataFrame(stats, columns=["method", "variable", "rmse", "max_err", "size", "cr"])
    return stats_df




if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input_file_pl", type=str, default="era5_pl_sample.nc")
    parser.add_argument("--input_file_sfc", type=str, default="era5_sfc_sample.nc")
    parser.add_argument("--output_folder", type=str, default="./compressed")
    parser.add_argument("--rel_target", type=float, default=0.001)
    parser.add_argument("--variable", type=str, default="z")
    parser.add_argument("--method", type=str, default="jp2spwv")
    args = parser.parse_args()
    ds_pl = xr.open_dataset(args.input_file_pl)
    ds_sfc = xr.open_dataset(args.input_file_sfc)
    output_folder = args.output_folder + f"/{args.rel_target}"
    os.makedirs(output_folder, exist_ok=True)
    stats_df_list = []
    for level in LEVELS:
        stats_df = compress_and_compare(ds_pl, args.variable, level, output_folder, args.method, args.rel_target)
        stats_df_list.append(stats_df)
    stats_df = pd.concat(stats_df_list)
    print(stats_df)

#python compress_and_compare.py --rel_target=0.01 --variable=w --method=jp2spwv