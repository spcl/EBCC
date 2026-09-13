#!/bin/bash

base_folder="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${PYTHON:-$base_folder/.venv/bin/python}"
HDF5_PLUGIN_PATH="$("$python_bin" -c 'from ebcc import EBCC_FILTER_DIR; print(EBCC_FILTER_DIR)')" || exit 1
export HDF5_PLUGIN_PATH

cdo -b F32 --filter "$("$python_bin" -m ebcc.filter_wrapper --base_cr 30 --height 721 --width 1440 -m 10)" copy geopotential_pl_small.nc geopotential_pl_small_cdo_ebcc.nc
ncdump -v z geopotential_pl_small_cdo_ebcc.nc | head -n 100
