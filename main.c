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

#include "dnestvars.h"


int main(int argc, char **argv)
{
  int opt;
  double t0=0.0, t1=0.0, dt;
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
    recon_flag_postprc = 0; /* default value, 0 means postprocessing after runing MCMC sampling. */
    recon_temperature = 1.0; /* default value */
    recon_flag_restart = 0;
    recon_flag_sample_info = 0;
    recon_flag_temp = 0;
    recon_flag_sim = 0;
    recon_flag_end = 0;
    recon_flag_limits=0;
    recon_flag_prior_exam=0;

    while( (opt = getopt(argc, argv, "pgt:rchle")) != -1)
    {
      switch(opt)
      {
        case 'p':  /* only do postprocessing */
          recon_flag_postprc = 1;
          recon_temperature = 1.0;
          printf("# MCMC samples available, only do post-processing.\n");
          break;
        case 'g':  /* generate mock light curve */
          recon_flag_sim = 1;
          printf("# generate mock time series.\n");
          break;
        case 't': /* temperature for postprocessing */
          recon_flag_temp = 1;
          recon_temperature = atof(optarg);
          printf("# Set a temperature %f.\n", recon_temperature);
          if(recon_temperature == 0.0)
          {
            printf("# Incorrect option -t %s.\n", optarg);
            exit(0);
          }
          if(recon_temperature < 1.0)
          {
            printf("# Temperature should >= 1.0\n");
            exit(0);
          }
          break;
        case 'r':   /* restart */
          recon_flag_restart = 1;
          printf("# Restart run.\n");
          break;

        case 'c': 
          printf("# Recalculate the sample info.\n");
          recon_flag_sample_info = 1;
          break;

        case 'l':
          recon_flag_limits = 1;
          //printf("# level-dependent sampling.\n");
          break;

        case 'e':
          recon_flag_prior_exam = 1;
          printf("# Examine priors.\n");
          break;

        case 'h':
          recon_flag_help = 1;
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

    if(recon_flag_help == 0) // not only print help.
    {
      if(argv[optind] != NULL) // parameter file is specified 
        strcpy(parset.param_file, argv[optind]); /* copy input parameter file */
      else
      {
        recon_flag_end = 1;
        fprintf(stderr, "# Error: No parameter file specified!\n");
      }
    }
  }

  /* broadcast flags */
  MPI_Bcast(&recon_flag_sim, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_postprc, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_restart, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_temp, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_temperature, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_sample_info, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_help, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_end, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_limits, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&recon_flag_prior_exam, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  if(recon_flag_end == 1 && recon_flag_help == 0 )
  {
    if(thistask == roottask)
    {
      fprintf(stdout, "Ends incorrectly.\n");
    }

    MPI_Finalize();
    return 0;
  }

  /* run the code */
  if(recon_flag_help == 0)
  {
    read_parset();
    recon();
  }
  

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

  MPI_Finalize();
  return 0;
}