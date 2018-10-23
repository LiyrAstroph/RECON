
from libc.string cimport *  
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from cython.operator cimport dereference

cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

import numpy as np
cimport numpy as np
np.import_array()

cdef class cyrecon:
  cdef int argc               # argc and argv for dnest
  cdef char **argv

  def __cinit__(self):

    # setup argc and argv
    cdef int i
    self.argv = <char **>PyMem_Malloc(9*sizeof(char *))
    for i in range(9):
      self.argv[i] = <char *>PyMem_Malloc(200*sizeof(char))
    
    self.argc = 0
    self.argv[self.argc] = 'recon'
    self.argc += 1
    self.argv[self.argc] = '-s'
    self.argc += 1
    self.argv[self.argc] = 'resetart_dnest.txt'
    self.argc += 1
    self.argv[self.argc] = '-l'
    self.argc += 1
    self.argv[self.argc] = 'param'
    self.argc += 1

    return

  def __cdealloc__(self):
    cdef int i
    for i in range(9):
      PyMem_Free(self.argv[i])
    
    PyMem_Free(self.argv)
  
    return

  def run(self):
    logz = recon_run(self.argc, self.argv)
    return logz