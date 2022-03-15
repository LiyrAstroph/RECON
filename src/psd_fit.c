/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

/*!
 *  \file psd_fit.c
 *  \brief perform periodogram fitting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
#include <string.h>

#include "dnest.h"

#include "allvars.h"
#include "proto.h"

/* function set for DNest */
DNestFptrSet *fptrset_fit;

double psd_fit()
{
  int i, argc=0, narg=9;
  char **argv;
  double logz = 0.0;
  
  strcpy(options_file, "OPTIONSFIT");

  argv = malloc(narg*sizeof(char *));
  for(i=0; i<narg; i++)
  {
    argv[i] = malloc(200*sizeof(char));
  }
  
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc++], "data/restart0_dnest.txt");

  if(recon_flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restart0_dnest.txt");
  }

  if(recon_flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
  }
  if(recon_flag_sample_info == 1)
  {
    strcpy(argv[argc++], "-c");
  }
  if(recon_flag_limits == 1)
  {
    strcpy(argv[argc++], "-l");
  }

  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], "0");

  psd_fit_init();

  logz = dnest(argc, argv, fptrset_fit, num_params, NULL, NULL, NULL, "./data/", options_file, NULL, NULL);
  
  psd_fit_postproc();

  psd_fit_end();

  for(i=0; i<narg; i++)
    free(argv[i]);

  free(argv);
  
  return logz;
}

/*!
 *  posterior process 
 */
void psd_fit_postproc()
{
  if(thistask == roottask)
  {
    int num_ps, i, j;
    char fname[200];
    void *posterior_sample, *post_model;
    FILE *fp;

    char posterior_sample_file[200];
    int size_of_modeltype = num_params * sizeof(double);
    double *pm;

    dnest_get_posterior_sample_file(posterior_sample_file);

    //file for posterior sample
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }

    // read numbers of points in posterior sample
    if(fscanf(fp, "# %d\n", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.%d %d\n", posterior_sample_file, i, j);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);
    }
    fclose(fp);

    /* calculate best parameter values */
    for(j=0; j<num_params; j++)
    {
      par_fit_best[j] = 0.0;
      par_fit_best_std[j] = 0.0;
    }
    for(i=0; i<num_ps; i++)
    {
      pm = (double *)(posterior_sample + i*size_of_modeltype);
      for(j=0; j<num_params; j++)
      {
        par_fit_best[j] += pm[j];
      }
    }
    for(j=0; j<num_params; j++)
    {
      par_fit_best[j] /= num_ps;
    }

    for(i=0; i<num_ps; i++)
    {
      pm = (double *)(posterior_sample + i*size_of_modeltype);
      for(j=0; j<num_params; j++)
      {
        par_fit_best_std[j] += pow(pm[j] - par_fit_best[j], 2);
      }
    }
    for(j=0; j<num_params; j++)
    {
      par_fit_best_std[j] = sqrt(par_fit_best_std[j]/num_ps);
      printf("best param %d: %f %f\n", j, par_fit_best[j], par_fit_best_std[j]);
    }

    free(post_model);
    free(posterior_sample);
  }

  MPI_Bcast(par_fit_best, num_params, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(par_fit_best_std, num_params, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  return;
}

/*!
 *  initiate fitting.
 */
void psd_fit_init()
{
  int i, j;

  /* number of parameters */
  num_params = num_params_psd_tot;

  DF = (freq_data[1] - freq_data[0]);

  fptrset_fit = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset_fit->from_prior = from_prior_fit;
  fptrset_fit->perturb = perturb_fit;
  fptrset_fit->log_likelihoods_cal = log_likelihoods_cal_fit;
  fptrset_fit->log_likelihoods_cal_initial = log_likelihoods_cal_initial_fit;
  fptrset_fit->log_likelihoods_cal_restart = log_likelihoods_cal_initial_fit;
  fptrset_fit->accept_action = accept_action_fit;
  fptrset_fit->kill_action = kill_action_fit;
  fptrset_fit->restart_action = restart_action_fit;

  /* array for parameter ranges */
  par_range_model = malloc((num_params)*sizeof(double *));
  for(i=0; i<num_params; i++)
  {
    par_range_model[i] = malloc(2*sizeof(double));
  }

  /* arrays for parameter fixing */
  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));

  set_par_range_fit();
  set_par_fix_fit();
  
  psd_data_model = malloc(nf_data * sizeof(double));
  psd_data_period_model = malloc(nf_data * sizeof(double));

  workspace_psd = (double *)malloc( 100 * sizeof(complex) );
  workspace_complex = (complex *)malloc( 100 * sizeof(complex) );

  // factor (2*PI*f*I)**k for CARMA process 
  if(parset.psd_model_enum == carma)
  {
    freq_array_pow = malloc((parset.carma_p) * sizeof(complex *));
    for(i=0; i<parset.carma_p; i++)
    {
      freq_array_pow[i] = (complex *)malloc(nf_data * sizeof(complex));
    }

    for(i=0; i<nf_data; i++)
    {
      freq_array_pow[0][i] = 2.0*PI*freq_data[i] * I;
      for(j=1; j<parset.carma_p; j++)
      {
        freq_array_pow[j][i] = freq_array_pow[j-1][i] * (2.0*PI*freq_data[i] * I);
      }
    }
  }

  return;
}

void psd_fit_end()
{
  int i;

  dnest_free_fptrset(fptrset_fit);
  
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(par_fix);
  free(par_fix_val);
  
  free(psd_data_model);
  free(psd_data_period_model);

  free(workspace_psd);
  free(workspace_complex);

  if(parset.psd_model_enum == carma)
  {
    for(i=0; i<parset.carma_p; i++)
    {
      free(freq_array_pow[i]);
    }
    free(freq_array_pow);
  }
}

void from_prior_fit(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
  }
  
  for(i=0; i<num_params; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
  return;
}

double perturb_fit(void *model)
{
  double logH=0.0, width;
  int which;
  double *pm = (double *)model;
 
  do
  {
    which = dnest_rand_int(num_params);
  }while(par_fix[which] == 1);

  which_parameter_update = which;
  
  width = ( par_range_model[which][1] - par_range_model[which][0] );

  pm[which] += dnest_randh() * width;
  dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);

  return logH;
}

double log_likelihoods_cal_fit(const void *model)
{
  double logL;
  logL = prob_fit(model);
  return logL;
}

double log_likelihoods_cal_initial_fit(const void *model)
{
  double logL;
  logL = prob_fit(model);
  return logL;
}

double prob_fit(const void *model)
{
  int i, i0;
  double *arg;
  double *pm = (double *)model;
  double prob = 0.0, x1, x2, P;

  arg = pm;
  
  psdfunc_array(freq_data, arg, psd_data_model, nf_data);
  
  if(parset.psdperiod_enum > delta)
  {
    psdfunc_period_array(freq_data, arg+parset.num_params_psd, psd_data_period_model, nf_data);
    for(i=0; i<nf_data; i++)
    {
      x1 = psd_data[i]/psd_data_model[i];
      x2 = psd_data_period_model[i]/psd_data_model[i];

      if(x2 < 0.01)
      {
        prob += log(x1) - (x1+x2);
      }
      else if(x1*x2 > 10.0)
      {
        prob += log(gsl_sf_bessel_In_scaled(0, 2.0*sqrt(x1*x2))) + 2.0*sqrt(x1*x2)
              + log(x1) - (x1+x2);
      }
      else
      {
        prob += log(gsl_sf_bessel_In(0, 2.0*sqrt(x1*x2)))
              + log(x1) - (x1+x2);
      }
    }
  }
  else if(parset.psdperiod_enum > none)
  {
    i0 = (int)( exp(pm[num_params_psd_tot - 2]) / DF );
    P = exp(2.0*pm[num_params_psd_tot-3]);
    for(i=0; i<nf_data; i++)
    {
      x1 = psd_data[i]/psd_data_model[i];
      prob += log(x1) - x1;
      if(i==i0)
      {
        x2 = P/psd_data_model[i];
        if(x2 < 0.01)
        {
          prob += log(x1) - (x1+x2);
        }
        else if(x1*x2 > 10.0)
        {
          prob += log(gsl_sf_bessel_In_scaled(0, 2.0*sqrt(x1*x2))) + 2.0*sqrt(x1*x2)
                + log(x1) - (x1+x2);
        }
        else
        {
          prob += log(gsl_sf_bessel_In(0, 2.0*sqrt(x1*x2)))
                + log(x1) - (x1+x2);
        }
      }
    }
  }
  else 
  {
    for(i=0; i<nf_data; i++)
    {
      x1 = psd_data[i]/psd_data_model[i];
      prob += log(x1) - x1;
    }
  }

  return prob;
}

void set_par_range_fit()
{
  int i;
  // variability parameters
  for(i=0; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[i][0];
    par_range_model[i][1] = var_range_model[i][1];
  }

  return;
}

void set_par_fix_fit()
{
  int i;
  //setup fixed parameters
  for(i=0; i<num_params; i++)
    par_fix[i] = 0;

  if(parset.flag_whitenoise == 0)
  {
    if(parset.psdperiod_enum == none) // no periodic component
    {
      par_fix[num_params_psd_tot-1] = 1;
      par_fix_val[num_params_psd_tot-1] = -DBL_MAX;
    }
    else // periodic component included
    {
      par_fix[num_params_psd_tot-3 - 1] = 1;
      par_fix_val[num_params_psd_tot-3 - 1] = -DBL_MAX;
    }
    
    if(thistask == roottask)
    {
      printf("# Exclude white noise.\n");
    }
  }

  return;
}

void restart_action_fit(int iflag)
{
  return;
}

void kill_action_fit(int i, int i_copy)
{
  return;
}

void accept_action_fit()
{
  return;
}