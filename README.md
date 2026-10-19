EBCC: Error Bounded Climate Compressor
======================================

[![arXiv](http://img.shields.io/badge/cs.CE-arXiv%3A2510.22265-B31B1B.svg)](https://arxiv.org/abs/2510.22265)

# Prebuilt binaries

Prebuilt wheels and HDF5 plugin archives are available on the [GitHub Releases](https://github.com/spcl/EBCC/releases) page for Linux (amd64, arm64) and macOS (arm64).

**Python wheel** (recommended — includes everything):
```bash
pip install ebcc-*.whl
```

**HDF5 plugin only** (for use with CDO or other HDF5 applications):
```bash
unzip ebcc-h5plugin-<platform>.zip
export HDF5_PLUGIN_PATH=$(pwd)
```

# Build with pip

Source builds require a C compiler, CMake, and Git (including builds from source archives). See [OpenJPEG integration](src/openjpeg_support/README.md) for maintenance details.

```bash
pip install -e .
# or with optional dependencies
pip install -e ".[dev]"           # Development dependencies
pip install -e ".[zarr]"          # Zarr support (requires python 3.11+ for zarr 3.0+)
```

# Manual build

```bash
git clone --recurse-submodules https://github.com/spcl/EBCC.git
mkdir -p EBCC/src/build
cd EBCC/src/build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make && make install
# the compiled filter is stored in `src/build/lib/libh5z_ebcc.so`
```

# User guide

EBCC compresses float32 arrays with at least two dimensions using a JPEG2000 base layer and an optional sparse wavelet residual. Three modes are supported:

1. `FIXED_CR` (`"fixed_cr"` or `residual_opt=None`): compress only the base layer using `base_cr` as the target compression ratio relative to the original float32 data. No error bound is enforced; the achieved ratio depends on the data and chunk size.
2. `MAX_ERROR` (`"max_error_target"`): maximum absolute error in the units of the data.
3. `RELATIVE_ERROR` (`"relative_error_target"`): maximum absolute error divided by the data range of the entire chunk.

`FIXED_CR` requires a finite positive `base_cr`. The other modes require a finite positive error bound, verify the reconstructed data against that bound, and ignore `base_cr`. Files written by earlier released versions remain readable.

The last two chunk axes are the spatial dimensions `(height, width)`, each of which must be between 32 and 2047. The product of the remaining chunk dimensions must not exceed 127. The Python wrapper defaults to one spatial frame per chunk. Specify `chunks` to compress several levels or time steps together:

```python
from ebcc.filter_wrapper import EBCC_Filter

# For an array with axes (time, level, latitude, longitude):
options = EBCC_Filter(
    height=721, width=1440,
    residual_opt=("max_error_target", 0.01),
    data_dim=4, chunks=(2, 37, 721, 1440),
)
```

`EBCC_Filter` converts these settings into HDF5 filter parameters. See the [HDF5 example](examples/example_hdf5.py) and [Zarr example](examples/example_zarr.py) for complete usage.

For `FIXED_CR` mode, use `base_cr=30` and `residual_opt=("fixed_cr", 0)`. The command-line equivalent is `python ebcc/filter_wrapper.py --base_cr 30`. Supplying `-m` or `-r` selects an error-bound mode and ignores an explicit `--base_cr` with a warning. With no options, the CLI defaults to a relative error bound of 0.01.

# CDO integration

The latest version of CDO can be downloaded here: `https://code.mpimet.mpg.de/projects/cdo/files`. It can be installed with:
```
tar -xvf <CDO-VERSION>.tar.gz
cd <CDO-VERSION>
./configure --enable-netcdf4 --with-netcdf=yes --with-hdf5=yes
make
sudo make install
```

Make sure the necessary libraries are loaded with `export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH`. Installation paths and the variables `CFLAGS` and `CPPFLAGS` may need to be set in order for the configuration and the compilation to succeed.

Set `HDF5_PLUGIN_PATH` to the plugin directory and use `ebcc/filter_wrapper.py` to generate CDO filter options. For example, with a maximum absolute error of 0.5 and tiles of 721x1440 (used by ERA5):

```bash
HDF5_PLUGIN_PATH=<path/to/filter> cdo -b F32 --filter 308,721,1440,1128792064,1,1056964608 copy temperature.nc compressed.nc
# or
HDF5_PLUGIN_PATH=<path/to/filter> cdo -b F32 --filter $(python ebcc/filter_wrapper.py --height 721 --width 1440 -m 0.5) copy temperature.nc compressed.nc
```

As an alternative, the `setfilter` function is also supported. The function allows the user to specify a filter for every variable in a netcdf file. Prepare a file `myfilter` containing the filter specification:

```
temperature="308,721,1440,1128792064,1,1056964608"
```

Then run:

```
HDF5_PLUGIN_PATH=<path/to/filter> cdo -b F32 setfilter,filename=myfilter temperature.nc compressed.nc
```

**CAUTION** Make sure to set output precision to float32 in cdo using `-b F32`! Otherwise, undefined behavior will occur (Segmentation Fault or incorrect result). Also make sure the chunksize of input netcdf file is a multiple of the tile size (height, width). If not, please either change height & width or rechunk the file using `nccopy -c ...` .

# Extra configurations through environment variables

- `EBCC_KAPPA`: positive base-layer allocation weight for error-bound modes, default `2.25`; ignored in `FIXED_CR` mode. Larger values retain more JPEG2000 coding passes. The final error bound is still verified; invalid values fall back to the default.
- `EBCC_DISABLE_MCT=1`: disable the component-axis transform during encoding. Other values keep MCT enabled for multicomponent chunks. Decoding requires no matching encoder settings.
- `EBCC_LOG_LEVEL`: integer from 0 to 5, default `3`: TRACE, DEBUG, INFO, WARN, ERROR, FATAL.
