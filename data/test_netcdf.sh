#!/bin/bash

export HDF5_PLUGIN_PATH=$(pwd)/../src/build/lib

cdo -b F32 --filter $(python ../filter_wrapper.py --base_cr 30 --height 721 --width 1440 -m 10) copy geopotential_pl_small.nc geopotential_pl_small_cdo_ebcc.nc
ncdump -v z geopotential_pl_small_cdo_ebcc.nc | head -n 100
