import numpy as np
import xarray as xr
from numcodecs.abc import Codec
import numcodecs
import ctypes
from numpy.ctypeslib import ndpointer

class CodecConfigT(ctypes.Structure):
    _fields_ = [
        ('dims', ctypes.c_size_t * 3),
        ('base_cr', ctypes.c_float),
        ('residual_compression_type', ctypes.c_int),
        ('residual_cr', ctypes.c_float),
        ('error', ctypes.c_float),
        ('quantile', ctypes.c_double),
    ]

class EBCCZarrFilter(Codec):
    codec_id = 'ebcc_filter'
    
    def __init__(self, arglist):
        self.arglist = np.array(arglist, dtype=np.uint32)
        self.c_stdlib = ctypes.CDLL('libc.so.6')
        self.lib = ctypes.CDLL('libh5z_j2k.so')
        self.lib.populate_config.argtypes = [
            ctypes.POINTER(CodecConfigT),
            ctypes.c_size_t,
            ndpointer(ctypes.c_uint, flags="C_CONTIGUOUS"),
            ctypes.c_size_t
        ]
        self.lib.encode_climate_variable.restype = ctypes.c_size_t
        self.lib.encode_climate_variable.argytpes = [
            ctypes.c_void_p,
            ctypes.POINTER(CodecConfigT),
            ctypes.POINTER(ctypes.c_void_p)
        ]
        self.lib.decode_climate_variable.restype = ctypes.c_size_t
        self.lib.decode_climate_variable.argytpes = [
            ctypes.c_void_p,
            ctypes.POINTER(CodecConfigT),
            ctypes.POINTER(ctypes.c_void_p)
        ]

    def encode(self, buf):
        out_buffer = ctypes.c_void_p()
        codec_config = CodecConfigT()
        self.lib.populate_config(ctypes.byref(codec_config), len(self.arglist), self.arglist, len(buf))

        encoded_size = self.lib.encode_climate_variable(buf, ctypes.byref(codec_config), ctypes.byref(out_buffer))

        out_array = ctypes.cast(out_buffer, ctypes.POINTER(ctypes.c_ubyte * encoded_size)).contents

        out_bytes = bytes(out_array)

        self.c_stdlib.free(out_buffer)

        return out_bytes

    def decode(self, buf, out=None):
        if out is not None:
            out_buffer = out.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        else:
            out_buffer = ctypes.POINTER(ctypes.c_float)()
        
        decoded_size = self.lib.decode_climate_variable(buf, len(buf), ctypes.byref(out_buffer))
        
        if out is not None:
            return out
        
        out_array = np.ctypeslib.as_array(out_buffer, shape=(decoded_size,))
    
        self.c_stdlib.free(out_buffer)
        
        return out_array.tobytes()

    def get_config(self):
        return {'id': self.codec_id, 'arglist': self.arglist}

    @classmethod
    def from_config(cls, config):
        return cls(config['arglist'])

# Register the custom filter
numcodecs.register_codec(EBCCZarrFilter)