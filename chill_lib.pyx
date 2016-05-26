# cython function that takes a numpy array of particle positions and returns their chill+ status.

import numpy as np
cimport numpy as np
from libcpp.vector cimport vector


cdef extern from "chiller.h":
  cdef cppclass Chiller:
    Chiller() except +
    Chiller(double*, int, double*) except +
    vector[int] getStatus()

cdef class PyChiller:
  cdef Chiller c_chill
  def __cinit__(self, \
                np.ndarray[double, ndim=2, mode="c"] positions not None, int numParticles, \
                np.ndarray[double, mode="c"] simulationCell not None):
    self.c_chill = Chiller(&positions[0,0], numParticles, &simulationCell[0])
  def get_status(self):
    return self.c_chill.getStatus()


print("Hello world")
