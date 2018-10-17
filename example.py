#
# an example from the DNest4 package by Brendon J. Brewer, with minor modifications
# 
#

from mpi4py import MPI
import numpy as np
import cyrecon

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# create a dnest sampler
sample = cyrecon.cyrecon()

# run sampler
sample.run()

comm.Barrier()
