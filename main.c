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
  double t0, t1, dt;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);

  if(thistask == roottask)
  {
    t0 = second();
    printf("===============RECON==================\n");
    printf("Starts to run...\n");
    printf("%d cores used.\n", totaltask);
  }

  /* cope with command argument */
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
    flag_end = 0;

    while( (opt = getopt(argc, argv, "pgt:rch")) != -1)
    {
      switch(opt)
      {
        case 'p':  /* only do postprocessing */
          flag_postprc = 1;
          temperature = 1.0;
          printf("# MCMC samples available, only do post-processing.\n");
          break;
        case 'g':  /* generate mock light curve */
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
        case 'r':   /* restart */
          flag_restart = 1;
          printf("# Restart run.\n");
          break;

        case 'c': 
          printf("# Recalculate the sample info.\n");
          flag_sample_info = 1;
          break;

        case 'h':
          flag_help = 1;
          print_help();
          break;

        case '?':
          printf("# Incorrect option -%c %s.\n", optopt, optarg);
          exit(0);
          break;

        default:
          break;
      }
    }

    if(flag_help == 0) // not only print help.
    {
      if(argv[optind] != NULL) // parameter file is specified 
        strcpy(parset.param_file, argv[optind]); /* copy input parameter file */
      else
      {
        flag_end = 1;
        fprintf(stderr, "# Error: No parameter file specified!\n");
      }
    }
  }

  /* broadcast flags */
  MPI_Bcast(&flag_sim, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_postprc, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_restart, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_temp, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_sample_info, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_help, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_end, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  if(flag_end == 1 && flag_help == 0 )
  {
    if(thistask == roottask)
    {
      fprintf(stdout, "Ends incorrectly.\n");
    }

    MPI_Finalize();
    return 0;
  }

  /* run the code */
  if(flag_help == 0)
  {
    read_parset();
    recon();
  }
  
  MPI_Finalize();

  if(thistask == roottask)
  {
    int ht, mt;
    double st;
    t1 = second();
    dt = timediff(t0, t1);
    get_hms(dt, &ht, &mt, &st);
    printf("Time used: %dh %dm %fs.\n", ht, mt, st);
    printf("Ends successfully.\n");
    printf("===============RECON==================\n");
  }
  return 0;
}