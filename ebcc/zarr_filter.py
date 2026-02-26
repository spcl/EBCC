import ctypes
import numpy as np
import numcodecs
from numcodecs.abc import Codec
from numpy.ctypeslib import ndpointer
from ebcc import EBCC_FILTER_PATH


class CodecConfigT(ctypes.Structure):
    _fields_ = [
        ('dims', ctypes.c_size_t * 3),
        ('base_param', ctypes.c_float),
        ('base_compressor', ctypes.c_int),
        ('residual_compression_type', ctypes.c_int),
        ('residual_cr', ctypes.c_float),
        ('error', ctypes.c_float),
    ]


class EBCCZarrFilter(Codec):
    codec_id = 'ebcc_j2k'
    _COMPRESSOR_TO_ID = {'j2k': 0, 'jxl': 1}

    def __init__(self, arglist, base_compressor='j2k'):
        self.arglist = np.array(arglist, dtype=np.uint32)
        self.base_compressor = base_compressor.lower()
        if self.base_compressor not in self._COMPRESSOR_TO_ID:
            raise ValueError("base_compressor must be 'j2k' or 'jxl'")

        self.lib = ctypes.CDLL(EBCC_FILTER_PATH)
        self.lib.populate_config.argtypes = [
            ctypes.POINTER(CodecConfigT),
            ctypes.c_size_t,
            ndpointer(ctypes.c_uint, flags='C_CONTIGUOUS'),
            ctypes.c_size_t,
        ]
        self.lib.ebcc_encode.restype = ctypes.c_size_t
        self.lib.ebcc_encode.argtypes = [
            ndpointer(ctypes.c_float, flags='C_CONTIGUOUS'),
            ctypes.POINTER(CodecConfigT),
            ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
        ]
        self.lib.ebcc_decode.restype = ctypes.c_size_t
        self.lib.ebcc_decode.argtypes = [
            ctypes.c_void_p,
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.POINTER(ctypes.c_float)),
        ]
        self.lib.free_buffer.argtypes = [ctypes.c_void_p]

    @property
    def codec_name(self):
        return 'ebcc_jxl' if self.base_compressor == 'jxl' else 'ebcc_j2k'

    def encode(self, buf):
        assert isinstance(buf, np.ndarray), 'Input buffer must be a numpy array'
        assert buf.dtype == np.float32, 'Input buffer must be of dtype float32'

        out_buffer = ctypes.POINTER(ctypes.c_ubyte)()
        codec_config = CodecConfigT()
        flat = np.ascontiguousarray(buf, dtype=np.float32).ravel()

        self.lib.populate_config(ctypes.byref(codec_config), len(self.arglist), self.arglist, flat.nbytes)
        codec_config.base_compressor = self._COMPRESSOR_TO_ID[self.base_compressor]

        encoded_size = self.lib.ebcc_encode(flat, ctypes.byref(codec_config), ctypes.byref(out_buffer))
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
        return {
            'id': self.codec_name,
            'arglist': self.arglist.astype(int).tolist(),
            'base_compressor': self.base_compressor,
        }

    @classmethod
    def from_config(cls, config):
        compressor = config.get('base_compressor')
        if compressor is None:
            compressor = 'jxl' if config.get('id') == 'ebcc_jxl' else 'j2k'
        return cls(config['arglist'], base_compressor=compressor)


class EBCCZarrFilterJXL(EBCCZarrFilter):
    codec_id = 'ebcc_jxl'

    def __init__(self, arglist, base_compressor='jxl'):
        super().__init__(arglist=arglist, base_compressor=base_compressor)


class EBCCZarrFilterLegacy(EBCCZarrFilter):
    codec_id = 'ebcc_filter'

    def __init__(self, arglist, base_compressor='j2k'):
        super().__init__(arglist=arglist, base_compressor=base_compressor)


numcodecs.register_codec(EBCCZarrFilter)
numcodecs.register_codec(EBCCZarrFilterJXL)
numcodecs.register_codec(EBCCZarrFilterLegacy)
