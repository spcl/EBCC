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
    chunksizes = (1, 1, nlat, nlon)
    chunksizes = encoding.pop("chunks") if "chunks" in encoding else chunksizes
    encoding_nc = {name: {**encoding, "chunksizes": chunksizes}}
    ds.to_netcdf(output_name, encoding=encoding_nc, engine="h5netcdf")

def compare(ds_orig, filename):
    ds_comp = xr.open_dataarray(filename, engine="h5netcdf")
    delta = ds_comp - ds_orig
    rmse = np.sqrt(np.square(delta).mean()).item()
    max_diff = np.max(np.abs(delta)).item()
    size = os.path.getsize(filename) / 1024 / 1024
    return rmse, max_diff, size

def compress_pl(ds, output_folder, method, rel_target):
    height = ds['latitude'].size
    width = ds['longitude'].size
    data_dim = len(ds['z'].shape)
    stats = []
    for var in ds.data_vars:
        data_range = (ds[var].max(dim=["latitude", "longitude"]) - ds[var].min(dim=["latitude", "longitude"])).mean().item()
        abs_target = data_range * rel_target
        print(f"Compressing {var} with {method}, relative error {rel_target}, absolute error {abs_target}")
        encoding = prepare_encoding(method, abs_target, height, width, data_dim)
        output_path = os.path.join(output_folder, f"era5_pl_{var}_{method}.nc")
        compress(ds[var], var, output_path, encoding)
        rmse, max_diff, size = compare(ds[var], output_path)
        stats.append([method, var, rmse, max_diff, size])
    stats_df = pd.DataFrame(stats, columns=["method", "variable", "rmse", "max_err", "size"])
    return stats_df



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--input_file_pl", type=str, default="era5_pl_sample.nc")
    parser.add_argument("--input_file_sfc", type=str, default="era5_sfc_sample.nc")
    parser.add_argument("--output_folder", type=str, default="./compressed")
    parser.add_argument("--rel_target", type=float, default=0.001)
    args = parser.parse_args()
    ds_pl = xr.open_dataset(args.input_file_pl)
    ds_sfc = xr.open_dataset(args.input_file_sfc)
    output_folder = args.output_folder + f"/{args.rel_target}"
    os.makedirs(output_folder, exist_ok=True)
    stats_df_list = []
    for method in METHODS:
        stats_df = compress_pl(ds_pl, output_folder, method, args.rel_target)
        stats_df_list.append(stats_df)
        print(stats_df)
    stats_df = pd.concat(stats_df_list)
    stats_df["rel_target"] = args.rel_target
    stats_df.to_csv(os.path.join(args.output_folder, f"stats_{args.rel_target}.csv"), index=False)