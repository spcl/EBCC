import numpy as np
import sys
from numcodecs.abc import Codec
import numcodecs
import ctypes
from numpy.ctypeslib import ndpointer
from ebcc import EBCC_FILTER_PATH

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
        if sys.platform.startswith('linux'):
            self.c_stdlib = ctypes.CDLL('libc.so.6')
        elif sys.platform == 'darwin':
            self.c_stdlib = ctypes.CDLL('libc.dylib')
        else:
            raise RuntimeError("Unsupported platform: " + sys.platform)
        self.lib = ctypes.CDLL(EBCC_FILTER_PATH)
        self.lib.populate_config.argtypes = [
            ctypes.POINTER(CodecConfigT),
            ctypes.c_size_t,
            ndpointer(ctypes.c_uint, flags="C_CONTIGUOUS"),
            ctypes.c_size_t
        ]
        self.lib.ebcc_encode.restype = ctypes.c_size_t
        self.lib.ebcc_encode.argtypes = [
            ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
            ctypes.POINTER(CodecConfigT),
            ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte))
        ]
        self.lib.ebcc_decode.restype = ctypes.c_size_t
        self.lib.ebcc_decode.argtypes = [
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.POINTER(ctypes.c_float))
        ]
        self.lib.free_buffer.argtypes = [ctypes.c_void_p]

    def encode(self, buf):
        assert isinstance(buf, np.ndarray), "Input buffer must be a numpy array"
        assert buf.dtype == np.float32, "Input buffer must be of dtype float32"
        out_buffer = ctypes.POINTER(ctypes.c_ubyte)()
        codec_config = CodecConfigT()
        buf = np.ascontiguousarray(buf, dtype=np.float32).ravel()
        self.lib.populate_config(ctypes.byref(codec_config), len(self.arglist), self.arglist, buf.nbytes)
        encoded_size = self.lib.ebcc_encode(buf, ctypes.byref(codec_config), ctypes.byref(out_buffer))

        out_array = ctypes.cast(out_buffer, ctypes.POINTER(ctypes.c_ubyte * encoded_size)).contents

        out_bytes = bytes(out_array)

        self.lib.free_buffer(out_buffer)

        return out_bytes

    def decode(self, buf, out=None):
        if out is not None:
            out_buffer = out.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        else:
            out_buffer = ctypes.POINTER(ctypes.c_float)()
        array_type = ctypes.c_ubyte * len(buf)
        buf_c = array_type.from_buffer_copy(buf)
        decoded_size = self.lib.ebcc_decode(ctypes.byref(buf_c), len(buf), ctypes.byref(out_buffer))

        if out is not None:
            return out
        
        out_array = np.ctypeslib.as_array(out_buffer, shape=(decoded_size,))
        out_array = out_array.copy()
        self.lib.free_buffer(out_buffer)
        
        return out_array

    def get_config(self):
        return {'id': self.codec_id, 'arglist': self.arglist.astype(int).tolist()}

    @classmethod
    def from_config(cls, config):
        return cls(config['arglist'])

# Register the custom filter
numcodecs.register_codec(EBCCZarrFilter)