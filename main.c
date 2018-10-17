/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

/*! \file main.c
 *  \brief start of the program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"


int main(int argc, char **argv)
{

  MPI_Init(&argc, &argv);

  recon_run(argc, argv);

  MPI_Finalize();
  return 0;
}