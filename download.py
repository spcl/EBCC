import xarray as xr

ds = xr.open_zarr('gs://weatherbench2/datasets/era5/1959-2023_01_10-full_37-1h-0p25deg-chunk-1.zarr')

ds.sel(level=[850]).isel(time=slice(-1, None))['temperature'].to_netcdf('temperature.nc', mode='w')
