EBCC: Error Bounded Climate Compressor
======================================

# How to install
Install required packages: `hdf5-devel`. Then:
```
git clone --recurse-submodules https://github.com/spcl/EBCC.git
mkdir compression-filter/src/build
cd compression-filter/src/build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make && make install
# the compiled filter is stored in `src/build/lib/libh5z_j2k.so`
```

# User guide
Our compression algorithm is implemented as an HDF5 filter. The filter expects a chunk size equal to a single 2D “frame”, and can scale to any number of dimensions greater or equan than 2.

Its base functionality compresses the data using JPEG2000. The user provides a "base" compression ratio for this functionality.
The user can also enable a compression of the residual in order to improve accuracy. The residual is the difference between the original frame and the uncompressed frame. The residual is compressed independently of the base frame, and is summed to the uncompressed base frame when decompressing. We support three modes of operation:
1. `NONE`: there is no residual
2. `MAX_ERROR`: the residual is wavelet encoded and sparsified. The sparsification factor is found through an iterative process that tries out several factors and selects the largest one that keeps the max error below the selected threshold
3. `RELATIVE_ERROR`: same as MAX_ERROR, but using the (data-range) relative error instead, rel_error = (x - ref) / (ref.max - ref.min)

The input parameters depend on the chosen mode of operation. As HDF5 filters support integer parameters only, a translation from float and double types to integer representation is required. We provide a python wrapper `EBCC_Filter` in `filter_wrapper.py` to simply the process. An example of how to use the filter using the python `h5py` library can be found in `test.py`

# CDO integration
The latest version of CDO can be downloaded here: `https://code.mpimet.mpg.de/projects/cdo/files`. It can be installed with:
```
tar -xvf <CDO-VERSION>.tar.gz
cd <CDO-VERSION>
./configure --enable-netcdf4 --with-netcdf=yes --with-hdf5=yes
make
sudo make install
```

Make sure the necessary libraries are loaded with `export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH`. Installation paths and the variables `CFLAGS` and `CPPFLAGS` may need to be set in order for the configuration and the compliation to succeed.

CDO can be now used with the `--filter` options as output by `filter_wrapper.py`. For HDF5 to see the filter, the `HDF5_PLUGIN_PATH` needs to be set.
As an example, with the default filter configuration and tiles of 721x1440 (used by ERA5):
```bash
HDF5_PLUGIN_PATH=<path/to/filter> cdo -b F32 --filter 308,721,1440,1128792064,3,1008981770 copy temperature.nc compressed.nc
# or
HDF5_PLUGIN_PATH=<path/to/filter> cdo -b F32 --filter $(python filter_wrapper.py --base_cr 30 --height 721 --width 1440 -m 0.5) copy temperature.nc compressed.nc
```

As an alternative, the `setfilter` function is also supported. The function allows the user to specify a filter for every variable in a netcdf file. Prepare a file `myfilter` containing the filter specification>
```
temperature="308,721,1440,1128792064,3,1008981770"
```

Then run:
```
HDF5_PLUGIN_PATH=<path/to/filter> cdo -b F32 setfilter,filename=myfilter temperature.nc compressed.nc
```
**CAUTION** Make sure to set output precision to float32 in cdo using `-b F32`! Otherwise, undefined behavior will occur (Segmentation Fault or incorrect result). Also make sure the chunksize of input netcdf file is a multiple of the tile size (height, width). If not, please either change height & width or rechunk the file using `nccopy -c ...` .

# Extra configurations through environment variables
- `EBCC_LOG_LEVEL`: valid value int [0, 5], default to 3, 0 - TRACE, 1 - DEBUG, 2 - INFO, 3 - WARN, 4 - ERROR, 5 - FATAL
- `EBCC_INIT_BASE_ERROR_QUANTILE`: valid value float [0, 1), default to 1e-6, set to 0 to turn off residual compression layer
- `EBCC_DISABLE_PURE_BASE_COMPRESSION_FALLBACK`: when set, turn off pure JP2 fallback (not recommended)
