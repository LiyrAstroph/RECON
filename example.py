
from mpi4py import MPI
import numpy as np
import cyrecon

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# create a dnest sampler
sample = cyrecon.cyrecon()

# run sampler
logz = sample.run()

if rank == 0:
  print(logz)

