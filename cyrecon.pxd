
#!/usr/bin/python
#cython: initializedcheck=False, boundscheck=False, wraparound=False, cdivision=True, profile=False

cdef extern from "PyFuncs.h":
  pass

cdef extern from "proto.h":
  double recon_run(int argc, char **argv)

cdef extern from "allvars.h":
  pass