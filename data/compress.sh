#!/bin/bash
#enstools-compression compress "geopotential_pl_small.nc" -o "geopotential_pl_small_sz.nc" --compression "lossy,sz,abs,10.0"
python compress_sz.py
python compress_jp2.py
python compress_sperr.py
python compare.py
