from ctypes import *
from ctypes.util import find_library
import numpy as np
import os

# flags copied from fftw3.h
FFTW_MEASURE = 0
FFTW_ESTIMATE = 1 << 6

V1_V2_INPUT_PARAMETERS = [
  {"size":     8192, "sparsity":   50},
  {"size":    16384, "sparsity":   50},
  {"size":    32768, "sparsity":   50},
  {"size":    65536, "sparsity":   50},
  {"size":   131072, "sparsity":   50},
  {"size":   262144, "sparsity":   50},
  {"size":   524288, "sparsity":   50},
  {"size":  1048576, "sparsity":   50},
  {"size":  2097152, "sparsity":   50},
  {"size":  4194304, "sparsity":   50},
  {"size":  8388608, "sparsity":   50},
  {"size": 16777216, "sparsity":   50},
  {"size":  4194304, "sparsity":   50},
  {"size":  4194304, "sparsity":  100},
  {"size":  4194304, "sparsity":  200},
  {"size":  4194304, "sparsity":  500},
  {"size":  4194304, "sparsity": 1000},
  {"size":  4194304, "sparsity": 2000},
  {"size":  4194304, "sparsity": 2500},
  {"size":  4194304, "sparsity": 4000},
]

class sfft:
  libsfft = cdll.LoadLibrary('libsfft.so')
  sfft_make_plan = libsfft["sfft_make_plan"]
  sfft_make_plan.restype = c_void_p
  sfft_make_plan.argtypes = [c_int, c_int, c_int, c_int]
  sfft_malloc = libsfft["sfft_malloc"]
  sfft_malloc.restype = c_void_p
  sfft_malloc.argtypes = [c_size_t]
  sfft_exec = libsfft["sfft_exec"]
  sfft_exec.restype = None
  sfft_exec.argtypes = [c_void_p, np.ctypeslib.ndpointer(np.complex_, ndim=1, flags='C'), \
			np.ctypeslib.ndpointer(np.complex_, ndim=1, flags='C')]
  sfft_free_plan = libsfft["sfft_free_plan"]
  sfft_free_plan.restype = None
  sfft_free_plan.argtypes = [c_void_p]
  sfft_free = libsfft["sfft_free"]
  sfft_free.restype = None
  sfft_free.argtypes = [c_void_p]

  def __init__(self, length=16384, sparsity=50, version=1, optimization=FFTW_ESTIMATE):
    if not isinstance(length, (int, long)):
      raise TypeError("length is not an integer")
    if not isinstance(sparsity, (int, long)):
      raise TypeError("sparsity is not an integer")
    if not isinstance(version, (int, long)):
      raise TypeError("version is not an integer")
    if not isinstance(optimization, (int, long)):
      raise TypeError("optimization is not an integer")

    if version not in [1, 2, 3]:
      raise ValueError("sFFT version %d is not valid.  Try 1, 2, or 3." % (version))
    if version in [1, 2] and {"size": length, "sparsity": sparsity} not in V1_V2_INPUT_PARAMETERS:
      raise ValueError("n = %d and k = %d is not a valid input parameter combination for sFFT version %d." % (length, sparsity, version))
    if optimization not in [FFTW_MEASURE, FFTW_ESTIMATE]:
      raise ValueError("FFTW optimization %d is not valid." % (optimization))

    self.length = length
    self.sparsity = sparsity
    self.version = version
    self.optimization = optimization

    self.sfft_plan = sfft.sfft_make_plan(c_int(length), c_int(sparsity), c_int(version-1), c_int(FFTW_ESTIMATE))

  def doit(self, a):
    requires = ['C', 'ALIGNED']
    a = np.asanyarray(a)
    a = np.require(a, np.complex_, requires)
    b = np.empty_like(a)
    sfft.sfft_exec(self.sfft_plan, a, b)
    return b

#  def __delete__(self):
#    sfft.sfft_free_plan(self.sfft_plan)
#    print "destroyed"
#
#    sfft.sfft_free(self.input_vector)
#    print "input vector freed"
#
#    sfft.sfft_free(self.output_vector)
#    print "output vector freed"
#
