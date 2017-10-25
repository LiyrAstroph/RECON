#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

#include "dnestvars.h"

int main(int argc, char **argv)
{
  int opt;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);

  if(thistask == roottask)
  {
    opterr = 0;
    optind = 0; // reset getopt. 
    flag_postprc = 0; /* default value, 0 means postprocessing after runing MCMC sampling. */
    temperature = 1.0; /* default value */
    flag_restart = 0;
    flag_sample_info = 0;
    flag_temp = 0;
    flag_sim = 0;
    
    while( (opt = getopt(argc, argv, "pgt:rch")) != -1)
    {
      switch(opt)
      {
        case 'p':  /* only do postprocessing */
          flag_postprc = 1;
          temperature = 1.0;
          printf("# MCMC samples available, only do post-processing.\n");
          break;
        case 'g':
          flag_sim = 1;
          printf("# generate mock time series.\n");
          break;
        case 't': /* temperature for postprocessing */
          flag_temp = 1;
          temperature = atof(optarg);
          printf("# Set a temperature %f.\n", temperature);
          if(temperature == 0.0)
          {
            printf("# Incorrect option -t %s.\n", optarg);
            exit(0);
          }
          if(temperature < 1.0)
          {
            printf("# Temperature should >= 1.0\n");
            exit(0);
          }
          break;
        case 'r':
          flag_restart = 1;
          printf("# Restart run.\n");
          break;

        case 'c':
          printf("# Recalculate the sample info.\n");
          flag_sample_info = 1;
          break;

        case 'h':
          flag_help = 1;
          //print_help();
          break;

        case '?':
          printf("# Incorrect option -%c %s.\n", optopt, optarg);
          exit(0);
          break;

        default:
          break;
      }
    }
  }

  MPI_Bcast(&flag_sim, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_postprc, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_restart, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_temp, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_sample_info, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  recon();

  MPI_Finalize();
  return 0;
}