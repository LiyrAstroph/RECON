/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>

#include "dnestvars.h"

#include "allvars.h"
#include "proto.h"

// function set for DNest.
DNestFptrSet *fptrset;

/*
 *========================================================================
 * do MCMC smapling
 *
 *========================================================================
 */

double recon()
{
  int i, argc=0, narg=9;
  char **argv;
  double logz = 0.0;
  
  strcpy(options_file, "OPTIONS");

  argv = malloc(narg*sizeof(char *));
  for(i=0; i<narg; i++)
  {
    argv[i] = malloc(200*sizeof(char));
  }
  
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc++], "data/restart_dnest.txt");

  if(recon_flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restart_dnest.txt");
  }

  if(recon_flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
  }
  if(recon_flag_temp == 1)
  {
    sprintf(argv[argc++], "-t%f", recon_temperature);
  }
  if(recon_flag_sample_info == 1)
  {
    strcpy(argv[argc++], "-c");
  }
  if(recon_flag_limits == 1)
  {
    strcpy(argv[argc++], "-l");
  }

  recon_init();

  if(recon_flag_sim == 1)
  {
    if(thistask == roottask)
    {
      sim();
    }
  }
  else
  {
    if(recon_flag_cal_psd != 1)
    {
      logz = dnest(argc, argv, fptrset, num_params, options_file);
      if(parset.flag_saveoutput == 0)
      {
        recon_postproc();
      }
    }
  }
  
  recon_end();

  for(i=0; i<narg; i++)
    free(argv[i]);

  free(argv);
  
  return logz;
}

/*
 *========================================================================
 * postprocess
 *
 *========================================================================
 */
int recon_postproc()
{
  if(thistask == roottask)
  {
    int num_ps, i, j, i1, i2;
    char fname[200];
    void *posterior_sample, *post_model;
    FILE *fp, *fcon, *fcon_all, *fcon_mean;

    char posterior_sample_file[200];
    int size_of_modeltype = num_params * sizeof(double);
    double *psd;

    get_posterior_sample_file(options_file, posterior_sample_file);

    //file for posterior sample
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    //file for continuum reconstruction
    sprintf(fname, "%s", "data/recon.txt");
    fcon = fopen(fname, "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    sprintf(fname, "%s", "data/recon_all.txt");
    fcon_all = fopen(fname, "w");
    if(fcon_all == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
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
    
    which_parameter_update = -1;
    which_particle_update = 0;

    psd = workspace_genlc[which_particle_update];

    for(j=0; j<nd_sim; j++)
    {
      flux_sim_mean[j] = 0.0;
      err_sim_mean[j] = 0.0;
    }

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
      
      genlc_array_initial(post_model);

      gsl_interp_init(gsl_linear_sim, time_sim, flux_sim, nd_sim);
  
      for(j=0; j<nd_sim; j++)
      {
        flux_sim_mean[j] += flux_sim[j];
        err_sim_mean[j] += flux_sim[j]*flux_sim[j];
      }

      i1 = ceil((time_data[0] - 0.2*(time_data[ndata-1] - time_data[0]) - time_sim[0])/DT);
      if(i1 < 0)
        i1 = 0;
      i2 = ceil((time_data[ndata-1]  + 0.2*(time_data[ndata-1] - time_data[0]) - time_sim[0])/DT);
      if(i2 > nd_sim)
        i2 = nd_sim;
      for(j=i1; j<i2; j++)
      {
        //flux_data_sim[j] = gsl_interp_eval(gsl_linear_sim, time_sim, flux_sim, time_data[j], gsl_acc_sim);

        fprintf(fcon, "%f  %f\n", time_sim[j], flux_sim[j]*flux_scale+flux_mean 
          + parset.slope_endmatch * (time_sim[j] - time_data[0]));
      }
      fprintf(fcon, "\n");

      for(j=0; j<nd_sim; j++)
      {
        fprintf(fcon_all, "%f  %f\n", time_sim[j], flux_sim[j]*flux_scale + flux_mean
          + parset.slope_endmatch * (time_sim[j] - time_data[0]) );
      }
      fprintf(fcon_all, "\n");

    }

    for(j=0; j<nd_sim; j++)
    {
      flux_sim_mean[j] /= num_ps;
      err_sim_mean[j] = sqrt( err_sim_mean[j]/num_ps - flux_sim_mean[j]*flux_sim_mean[j] );
    }

    sprintf(fname, "%s", "data/recon_mean.txt");
    fcon_mean = fopen(fname, "w");
    if(fcon_all == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    for(j=0; j<nd_sim; j++)
    {
      fprintf(fcon_mean, "%f %e %e\n", time_sim[j], flux_sim_mean[j]*flux_scale+flux_mean
        + parset.slope_endmatch * (time_sim[j] - time_data[0]), err_sim_mean[j]*flux_scale);
    }

    fclose(fp);
    fclose(fcon);
    fclose(fcon_all);
    fclose(fcon_mean);

    free(post_model);
    free(posterior_sample);
  }

  return 0;
}

/*
 *========================================================================
 * setup recon
 *
 *========================================================================
 */
int recon_init()
{
  char fname[200];
  int i, j;
  double Tall, Tmin, Tmax;
 
  /* set the root task */
  roottask = 0;

  fptrset = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset->from_prior = from_prior_recon;

  if(recon_flag_limits == 0)
    fptrset->perturb = perturb_recon;
  else
    fptrset->perturb = perturb_recon_limits;

  fptrset->restart_action = restart_action_recon;
  if(parset.psdperiod_enum > delta)
  {
    fptrset->accept_action = accept_action_recon_period;
    fptrset->kill_action = kill_action_recon_period;
  }
  else
  {
    fptrset->accept_action = accept_action_recon_other;
    fptrset->kill_action = kill_action_recon_other;
  }
  if(parset.flag_saveoutput == 0)
  {
    fptrset->print_particle = print_particle_recon;
    fptrset->read_particle = read_particle_recon;
  }
  else
  {
    fptrset->print_particle = print_particle_recon_saveoutput;
    fptrset->read_particle = read_particle_recon_saveoutput;
  }

  if(recon_flag_prior_exam == 0)
  {
    fptrset->log_likelihoods_cal = log_likelihoods_cal_recon;
    fptrset->log_likelihoods_cal_initial = log_likelihoods_cal_initial_recon;
    fptrset->log_likelihoods_cal_restart = log_likelihoods_cal_initial_recon;
  }
  else
  {
    fptrset->log_likelihoods_cal = log_likelihoods_cal_recon_exam;
    fptrset->log_likelihoods_cal_initial = log_likelihoods_cal_recon_exam;
    fptrset->log_likelihoods_cal_restart = log_likelihoods_cal_recon_exam;
  }

  set_psd_functions();

  // open file for output recon configurations 
  // read number of particles
  if(thistask == roottask)
  {
    sprintf(fname, "%s/data/recon_info.txt", parset.file_dir);
    finfo = fopen(fname, "w");
    if(finfo == NULL)
    {
      printf("Cannot open file %s.\n", fname);
      exit(0);
    }

    get_num_particles(options_file);
    get_max_num_levels(options_file);
  }
  MPI_Bcast(&num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&max_num_levels, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  num_params_psd_tot = parset.num_params_psd + parset.num_params_psdperiod;

  /* fft */
  if(recon_flag_sim==1)
  {
    V = 10;
    W = 10;
    time_media = 0.0;
    DT = parset.DT/W;
    nd_sim = parset.nd_sim * W * V;

    freq_limit_data_lower = 1.0/(nd_sim * DT/W);
    freq_limit_data_upper = 1.0/(DT/W)/2.0;

    read_sim_arg();
  }
  else
  {
  
    if(thistask == roottask)
    {
      sprintf(fname, "%s/%s", parset.file_dir, parset.file_name);
      ndata = get_line_number(fname);
    }
 
    MPI_Bcast(&ndata, 1, MPI_INT, roottask, MPI_COMM_WORLD);
 
    time_data = malloc(ndata * sizeof(double));
    flux_data = malloc(ndata * sizeof(double));
    err_data = malloc(ndata * sizeof(double));
    flux_data_sim = malloc(ndata * sizeof(double));
 
    if(thistask == roottask)
    {
      read_data(fname, ndata, time_data, flux_data, err_data);

      psddata_cal();
      
      time_cad_cal();
    }

    MPI_Bcast(time_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(flux_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(err_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    
    MPI_Bcast(&time_cad_media, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(&parset.slope_endmatch, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    
    time_media = (time_data[ndata-1] + time_data[0])/2.0;
    time_cad_min = (time_data[ndata-1] - time_data[0])/2.0;
    flux_data_min = flux_data_max = flux_data[0];
    flux_mean = flux_data[0];
    for(i=1; i<ndata; i++)
    {
      flux_mean += flux_data[i];
 
      if(time_data[i] - time_data[i-1] < time_cad_min)
        time_cad_min = time_data[i] - time_data[i-1];
 
      if(flux_data_min > flux_data[i])
        flux_data_min = flux_data[i];
 
      if(flux_data_max < flux_data[i])
        flux_data_max = flux_data[i];
    }
 
    flux_mean /= ndata;
    flux_scale = (flux_data_max - flux_data_min)/2.0;
    for(i=0; i<ndata; i++)
    {
      flux_data[i] = (flux_data[i] - flux_mean)/flux_scale;
      err_data[i] /= flux_scale;
    }

    flux_var = 0.0;
    for(i=0; i<ndata; i++)
    {
      flux_var += (flux_data[i] - flux_mean/flux_scale) * (flux_data[i] - flux_mean/flux_scale);
    }
    flux_var /= (ndata-1);

    err_mean = 0.0;
    for(i=0; i<ndata; i++)
    {
      err_mean += err_data[i];
    }
    err_mean /= ndata;

  
    V = parset.V;
    W = parset.W;
    Tall = time_data[ndata-1] - time_data[0];
    Tmin = time_data[0] - 0.5*(V-1.0)*Tall;
    Tmax = time_data[ndata-1] + 0.5*(V-1.0)*Tall;
    Tall = Tmax - Tmin;
  
    //DT = (time_data[ndata-1] - time_data[0])/(ndata -1)/W;
    //DT = time_cad_min/W;
    DT = time_cad_media/W;
    nd_sim = ceil(Tall/DT) + 1;
    nd_sim = (nd_sim/2) * 2; //make sure it is an even number.

    noise_power = 2.0*DT * err_mean*err_mean;

    if(thistask == roottask)
    {
      printf("N=%d, DT=%f, Tall=%f\n", nd_sim, DT, Tall);

      fprintf(finfo, "N=%d, DT=%f, Tall=%f\n", nd_sim, DT, Tall);
      fprintf(finfo, "time media: %f\n", time_media);
      fprintf(finfo, "flux mean: %f\n", flux_mean);
      fprintf(finfo, "flux scale:%f\n", flux_scale);
      fprintf(finfo, "min. data cad.: %f\n", time_cad_min);
      fprintf(finfo, "med. data cad.: %f\n", time_cad_media);
      fprintf(finfo, "mean data cad.: %f\n", (time_data[ndata-1] - time_data[0])/(ndata -1));

      fprintf(finfo, "end mathcing flag: %d\n", parset.flag_endmatch);
      fprintf(finfo, "end mathcing slope: %f\n", parset.slope_endmatch);
      fprintf(finfo, "level-dependent sampling flag: %d\n", recon_flag_limits);
    }

    freq_limit_data_lower = 1.0/(time_data[ndata-1] - time_data[0]);
    freq_limit_data_upper = ndata/(time_data[ndata-1] - time_data[0])/2.0;

    if(parset.freq_limit >= freq_limit_data_lower)
    {
      if(thistask == roottask)
      {
        printf("freq_limit is not good, larger than the minimal frequency limit of data.\n");
      }
      exit(0);
    }
  }

  num_recon = nd_sim;
  if(parset.psdperiod_enum > delta)
  {
    num_recon += nd_sim/2;
  }

  num_params = num_recon + num_params_psd_tot;

  var_range_model = malloc((num_params_psd_tot+2)*sizeof(double *));
  for(i=0; i<num_params_psd_tot+2; i++)
  {
    var_range_model[i] = malloc(2*sizeof(double));
  }

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  
  set_par_range();
  set_par_fix();

  // normalization factor in line with the continuous transform
  norm_psd = 1.0/sqrt(nd_sim) * sqrt(nd_sim/(2.0 *nd_sim * DT));
  // normalization factor for likelihood probability
  norm_prob = -0.5*ndata*log(2.0*PI);
  for(i=0; i<ndata; i++)
    norm_prob += (-log(err_data[i]));

  DF = 1.0/(nd_sim * DT);

  time_sim = malloc(nd_sim * sizeof(double));
  flux_sim = malloc(nd_sim * sizeof(double));
  flux_sim_mean = malloc(nd_sim * sizeof(double));
  err_sim_mean = malloc(nd_sim * sizeof(double));
  for(i=0; i<nd_sim; i++)
  {
    time_sim[i] = (i - nd_sim/2.0) * DT  + time_media;
  }
  
  gsl_acc_sim = gsl_interp_accel_alloc();
  gsl_linear_sim = gsl_interp_alloc(gsl_interp_linear, nd_sim);

  fft_work = (fftw_complex *)fftw_malloc( (nd_sim/2+1) * sizeof(fftw_complex));
  pfft = fftw_plan_dft_c2r_1d(nd_sim, fft_work, flux_sim, FFTW_MEASURE);

  freq_array = (double *)malloc(nd_sim/2*sizeof(double));
  workspace_psd = (double *)malloc( 100 * sizeof(complex) );
  workspace_complex = (complex *)malloc( 100 * sizeof(complex) );

  // initialize frequency grid array
  for(i=0; i<nd_sim/2; i++)
  {
    freq_array[i] = (i+1)*DF;
  }
  // initialize fft_work
  for(i=0; i<nd_sim/2+1; i++)
  {
    fft_work[i][0] = fft_work[i][1] = 0.0;
  }

  // determine the index for freq_limit in the frequency array
  if(parset.freq_limit <= freq_array[0] )
  {
    idx_limit = 0;
  }
  else if(parset.freq_limit >= freq_array[nd_sim/2 - 1])
  {
    idx_limit = nd_sim/2;
  }
  else
  {
    idx_limit = (int)((parset.freq_limit - freq_array[0])/(freq_array[1] - freq_array[0])) + 1;
  }

  // factor (2*PI*f*I)**k for CARMA process 
  if(parset.psd_model_enum == carma)
  {
    freq_array_pow = malloc((parset.carma_p) * sizeof(complex *));
    for(i=0; i<parset.carma_p; i++)
    {
      freq_array_pow[i] = (complex *)malloc(nd_sim/2 * sizeof(complex));
    }

    for(i=0; i<nd_sim/2; i++)
    {
      freq_array_pow[0][i] = 2.0*PI*freq_array[i] * I;
      for(j=1; j<parset.carma_p; j++)
      {
        freq_array_pow[j][i] = freq_array_pow[j-1][i] * (2.0*PI*freq_array[i] * I);
      }
    }
  }

  workspace_genlc = malloc(num_particles * sizeof(double *));
  workspace_genlc_perturb = malloc(num_particles * sizeof(double *));
  for(i=0; i<num_particles; i++)
  {
    workspace_genlc[i] = malloc(nd_sim/2 * sizeof(double));
    workspace_genlc_perturb[i] = malloc(nd_sim/2 * sizeof(double));
  }

  if(parset.psdperiod_enum > delta)
  {
    workspace_genlc_period = malloc(num_particles * sizeof(double *));
    workspace_genlc_period_perturb = malloc(num_particles * sizeof(double *));
    for(i=0; i<num_particles; i++)
    {
      workspace_genlc_period[i] = malloc(nd_sim/2 * sizeof(double));
      workspace_genlc_period_perturb[i] = malloc(nd_sim/2 * sizeof(double));
    }
  }

  if(thistask == roottask)
  {
    fclose(finfo);
    
    if(recon_flag_restart == 0)
    {
      remove_restart_file(); /* remove restart file. */
    }
  }

  return 0;
}

/*
 *========================================================================
 * finalize recon
 *
 *========================================================================
 */
int recon_end()
{
  int i;
 
  dnest_free_fptrset(fptrset);

  fftw_destroy_plan(pfft);

  gsl_interp_free(gsl_linear_sim);
  gsl_interp_accel_free(gsl_acc_sim);

  fftw_free(fft_work);
  
  for(i=0; i<num_params_psd_tot+2; i++)
  {
    free(var_range_model[i]);
  }
  free(var_range_model);

  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(par_fix);
  free(par_fix_val);
  
  free(time_sim);
  free(flux_sim);
  free(flux_sim_mean);
  free(err_sim_mean);

  if(recon_flag_sim != 1)
  {
    free(time_data);
    free(flux_data);
    free(err_data);
    free(flux_data_sim);
  }
  free(freq_array);
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

  for(i=0; i<num_particles; i++)
  {
    free(workspace_genlc[i]);
    free(workspace_genlc_perturb[i]);
  }
  free(workspace_genlc);
  free(workspace_genlc_perturb);

  if(parset.psdperiod_enum > delta)
  {
    for(i=0; i<num_particles; i++)
    {
      free(workspace_genlc_period[i]);
      free(workspace_genlc_period_perturb[i]);
    }
    free(workspace_genlc_period);
    free(workspace_genlc_period_perturb);
  }
  
  return 0;
}

/*
 *========================================================================
 * generate a light curve for a given model.
 *
 *========================================================================
 */
void genlc(const void *model)
{
  int i, i0;
  double *arg, freq, psd_sqrt;
  double *pm = (double *)model;

  arg = pm;

  fft_work[0][0] = pm[num_params_psd_tot+0]; //zero-frequency power.
  fft_work[0][1] = 0.0;

  for(i=0; i<nd_sim/2-1; i++)
  {
    freq = (i+1)*DF;
    psd_sqrt = psdfunc_sqrt(freq, arg)/sqrt(2.0);
    fft_work[i+1][0] = pm[num_params_psd_tot+1+2*i] * psd_sqrt;
    fft_work[i+1][1] = pm[num_params_psd_tot+1+2*i+1] * psd_sqrt;
  }
  freq = nd_sim/2*DF;
  psd_sqrt = psdfunc_sqrt(freq, arg);
  fft_work[nd_sim/2][0] = pm[num_params_psd_tot + nd_sim-1] * psd_sqrt;
  fft_work[nd_sim/2][1] = 0.0;
  
  /*for(i=1; i<nd_sim/2; i++)
  {
    freq = i*1.0/(nd_sim * DT);
    psd_sqrt = psdfunc_sqrt(freq, arg)/sqrt(2.0);
    fft_work[i][0] *= psd_sqrt;
    fft_work[i][1] *= psd_sqrt;
  }
  freq = nd_sim/2*1.0/(nd_sim * DT);
  psd_sqrt = psdfunc_sqrt(freq, arg);
  fft_work[nd_sim/2][0] *= psd_sqrt; */
  
  // add periodic component
  if(parset.psdperiod_enum > delta)
  {
    for(i=1; i<nd_sim/2; i++)
    {
      freq = i*DF;
      psd_sqrt = psdfunc_period_sqrt(freq, arg+num_params_psd_tot-3); // the last 3 vars
      fft_work[i][0] += psd_sqrt * cos(pm[num_params_psd_tot + nd_sim-1+i] * 2.0*PI);
      fft_work[i][1] += psd_sqrt * sin(pm[num_params_psd_tot + nd_sim-1+i] * 2.0*PI);
    }
    i = nd_sim/2;
    freq = i*DF;
    psd_sqrt = psdfunc_period_sqrt(freq, arg+num_params_psd_tot-3); // the last 3 vars
    fft_work[i][0] += psd_sqrt;
  }
  else if(parset.psdperiod_enum > none)
  {
    i0 = (int)( exp(pm[num_params_psd_tot - 2]) / DF );
    fft_work[i0][0] += exp(pm[num_params_psd_tot-3]) * cos(pm[num_params_psd_tot -1] * 2.0*PI);
    fft_work[i0][1] += exp(pm[num_params_psd_tot-3]) * sin(pm[num_params_psd_tot -1] * 2.0*PI);
  }
  
  fftw_execute(pfft);

  //normalization
  for(i=0; i<nd_sim; i++)
  {
    flux_sim[i] = flux_sim[i] * norm_psd;
  }

  return;
}

/*
 *========================================================================
 * fast version of genlc.
 * calculate PSD for a frequecy series at the same calling.
 *
 *========================================================================
 */
void genlc_array(const void *model)
{
  int i;
  double *arg, *psd_sqrt;
  double *pm = (double *)model;

  arg = pm;

  if(which_parameter_update < parset.num_params_psd)
  {
    psd_sqrt = workspace_genlc_perturb[which_particle_update];
    psdfunc_sqrt_array(freq_array, arg, psd_sqrt, nd_sim/2);
  }
  else
  {
    psd_sqrt = workspace_genlc[which_particle_update];
  }

  fft_work[0][0] = pm[num_params_psd_tot+0]; //zero-frequency power.
  fft_work[0][1] = 0.0;

  for(i=0; i<nd_sim/2-1; i++)
  {
    fft_work[i+1][0] = pm[num_params_psd_tot+1+2*i] * psd_sqrt[i]/sqrt(2.0);
    fft_work[i+1][1] = pm[num_params_psd_tot+1+2*i+1] *  psd_sqrt[i]/sqrt(2.0);
  }
  fft_work[nd_sim/2][0] = pm[num_params_psd_tot + nd_sim-1] * psd_sqrt[nd_sim/2-1];
  fft_work[nd_sim/2][1] = 0.0;
  
  fftw_execute(pfft);

  // normalization
  for(i=0; i<nd_sim; i++)
  {
    flux_sim[i] = flux_sim[i] * norm_psd;
  }

  return;
}

/*
 *========================================================================
 * fast version of genlc for cases of simple periodic PSD
 * calculate PSD for a frequecy series at the same calling.
 *
 *========================================================================
 */
void genlc_array_simple_period(const void *model)
{
  int i;
  double *arg, *psd_sqrt, *psdperiod_sqrt;
  double *pm = (double *)model;

  arg = pm;

  if(which_parameter_update < parset.num_params_psd)
  {
    psd_sqrt = workspace_genlc_perturb[which_particle_update];
    psdfunc_sqrt_array(freq_array, arg, psd_sqrt, nd_sim/2);

    psdperiod_sqrt = workspace_genlc_period[which_particle_update];
  }
  else if(which_parameter_update < num_params_psd_tot)
  {
    psd_sqrt = workspace_genlc[which_particle_update];

    psdperiod_sqrt = workspace_genlc_period_perturb[which_particle_update];
    psdfunc_period_sqrt_array(freq_array, arg+parset.num_params_psd, psdperiod_sqrt, nd_sim/2);
  }
  else
  {
    psd_sqrt = workspace_genlc[which_particle_update];
    psdperiod_sqrt = workspace_genlc_period[which_particle_update];
  }

  fft_work[0][0] = pm[num_params_psd_tot+0]; //zero-frequency power.
  fft_work[0][1] = 0.0;

  for(i=0; i<nd_sim/2-1; i++)
  {
    fft_work[i+1][0] = pm[num_params_psd_tot+1+2*i] * psd_sqrt[i]/sqrt(2.0);
    fft_work[i+1][1] = pm[num_params_psd_tot+1+2*i+1] *  psd_sqrt[i]/sqrt(2.0);
  }
  fft_work[nd_sim/2][0] = pm[num_params_psd_tot + nd_sim-1] * psd_sqrt[nd_sim/2-1];
  fft_work[nd_sim/2][1] = 0.0;

  // add periodic component
  for(i=1; i<nd_sim/2; i++)
  {
    fft_work[i][0] += psdperiod_sqrt[i-1] * cos(pm[num_params_psd_tot + nd_sim-1+i] * 2.0*PI);
    fft_work[i][1] += psdperiod_sqrt[i-1] * sin(pm[num_params_psd_tot + nd_sim-1+i] * 2.0*PI);
  }
  fft_work[nd_sim/2][0] += psdperiod_sqrt[nd_sim/2-1] ;
  
  fftw_execute(pfft);

  // normalization
  for(i=0; i<nd_sim; i++)
  {
    flux_sim[i] = flux_sim[i] * norm_psd;
  }

  return;
}

/*
 *========================================================================
 * fast version of genlc for cases of Delta periodic PSD
 * calculate PSD for a frequecy series at the same calling.
 *
 *========================================================================
 */
void genlc_array_simple_delta(const void *model)
{
  int i, i0;
  double *arg, *psd_sqrt;
  double *pm = (double *)model;

  arg = pm;

  if(which_parameter_update < parset.num_params_psd)
  {
    psd_sqrt = workspace_genlc_perturb[which_particle_update];
    psdfunc_sqrt_array(freq_array, arg, psd_sqrt, nd_sim/2);
  }
  else
  {
    psd_sqrt = workspace_genlc[which_particle_update];
  }

  fft_work[0][0] = pm[num_params_psd_tot+0]; //zero-frequency power.
  fft_work[0][1] = 0.0;

  for(i=0; i<nd_sim/2-1; i++)
  {
    fft_work[i+1][0] = pm[num_params_psd_tot+1+2*i] * psd_sqrt[i]/sqrt(2.0);
    fft_work[i+1][1] = pm[num_params_psd_tot+1+2*i+1] *  psd_sqrt[i]/sqrt(2.0);
  }
  fft_work[nd_sim/2][0] = pm[num_params_psd_tot + nd_sim-1] * psd_sqrt[nd_sim/2-1];
  fft_work[nd_sim/2][1] = 0.0;

  // add periodic component
  i0 = (int)( exp(pm[num_params_psd_tot - 2])/DF );
  
  fft_work[i0][0] += exp(pm[num_params_psd_tot-3]) * cos(pm[num_params_psd_tot -1] * 2.0*PI);
  fft_work[i0][1] += exp(pm[num_params_psd_tot-3]) * sin(pm[num_params_psd_tot -1] * 2.0*PI);
  
  fftw_execute(pfft);

  // normalization
  for(i=0; i<nd_sim; i++)
  {
    flux_sim[i] = flux_sim[i] * norm_psd;
  }

  return;
}

/*
 * action when perturbed move is accepted for periodical PSD cases.
 *
 */
void accept_action_recon_period()
{
  double *ptemp = NULL;
  int param;
  
  param = which_parameter_update;

  if(param < parset.num_params_psd)
  {
    ptemp = workspace_genlc[which_particle_update];
    workspace_genlc[which_particle_update] = workspace_genlc_perturb[which_particle_update];
    workspace_genlc_perturb[which_particle_update] = ptemp;
  }
  else if( param < num_params_psd_tot ) // only update when periodic psd parameters changed
  {
    ptemp = workspace_genlc_period[which_particle_update];
    workspace_genlc_period[which_particle_update] = workspace_genlc_period_perturb[which_particle_update];
    workspace_genlc_period_perturb[which_particle_update] = ptemp;
  }
  return;
}
/*
 * action when perturbed move is accepted for other PSD cases.
 *
 */
void accept_action_recon_other()
{
  double *ptemp = NULL;
  int param;
  
  param = which_parameter_update;
  if(param < parset.num_params_psd)
  {
    ptemp = workspace_genlc[which_particle_update];
    workspace_genlc[which_particle_update] = workspace_genlc_perturb[which_particle_update];
    workspace_genlc_perturb[which_particle_update] = ptemp;
  }
  return;
}

/*
 * action when DNest particle is killed for periodic PSD cases.
 */
void kill_action_recon_period(int i, int i_copy)
{
  memcpy(workspace_genlc[i], workspace_genlc[i_copy], nd_sim/2*sizeof(double));
  memcpy(workspace_genlc_period[i], workspace_genlc_period[i_copy], nd_sim/2*sizeof(double));
}
void kill_action_recon_other(int i, int i_copy)
{
  memcpy(workspace_genlc[i], workspace_genlc[i_copy], nd_sim/2*sizeof(double));
}

/*
 *========================================================================
 * fast version of genlc.
 * generate lc at the initial step.
 *
 *========================================================================
 */
void genlc_array_initial(const void *model)
{
  int i, i0;
  double *arg, *psd_sqrt, *psdperiod_sqrt;
  double *pm = (double *)model;

  arg = pm;
  psd_sqrt = workspace_genlc[which_particle_update];
  psdfunc_sqrt_array(freq_array, arg, psd_sqrt, nd_sim/2);
  
  fft_work[0][0] = pm[num_params_psd_tot+0]; //zero-frequency power.
  fft_work[0][1] = 0.0;

  for(i=0; i<nd_sim/2-1; i++)
  {
    fft_work[i+1][0] = pm[num_params_psd_tot+1+2*i] * psd_sqrt[i]/sqrt(2.0);
    fft_work[i+1][1] = pm[num_params_psd_tot+1+2*i+1] *  psd_sqrt[i]/sqrt(2.0);
  }
  fft_work[nd_sim/2][0] = pm[num_params_psd_tot + nd_sim-1] * psd_sqrt[nd_sim/2-1];
  fft_work[nd_sim/2][1] = 0.0;

  // add periodic component
  if(parset.psdperiod_enum > delta)
  {
    psdperiod_sqrt = workspace_genlc_period[which_particle_update];
    psdfunc_period_sqrt_array(freq_array, arg+parset.num_params_psd, psdperiod_sqrt, nd_sim/2);
    
    for(i=1; i<nd_sim/2; i++)
    {
      fft_work[i][0] += psdperiod_sqrt[i-1] * cos(pm[num_params_psd_tot + nd_sim-1+i] * 2.0*PI);
      fft_work[i][1] += psdperiod_sqrt[i-1] * sin(pm[num_params_psd_tot + nd_sim-1+i] * 2.0*PI);
    }
    fft_work[nd_sim/2][0] += psdperiod_sqrt[nd_sim/2-1] ;
  }
  else if(parset.psdperiod_enum > none)
  {
    i0 = (int)( exp(pm[num_params_psd_tot - 2]) / DF );
    fft_work[i0][0] += exp(pm[num_params_psd_tot-3]) * cos(pm[num_params_psd_tot -1] * 2.0*PI);
    fft_work[i0][1] += exp(pm[num_params_psd_tot-3]) * sin(pm[num_params_psd_tot -1] * 2.0*PI);
  }
  
  fftw_execute(pfft);

  // normalization
  for(i=0; i<nd_sim; i++)
  {
    flux_sim[i] = flux_sim[i] * norm_psd;
  }

  return;
}

/*
 *========================================================================
 * likelihood probability
 *
 *========================================================================
 */
double prob_recon(const void *model)
{
  double prob=0.0;
  int i;

  which_particle_update = dnest_get_which_particle_update();

  func_genlc_array(model);
  
  gsl_interp_init(gsl_linear_sim, time_sim, flux_sim, nd_sim);
  
  for(i=0; i<ndata; i++)
  {
    flux_data_sim[i] = gsl_interp_eval(gsl_linear_sim, time_sim, flux_sim, time_data[i], gsl_acc_sim);
  }

  for(i=0; i<ndata; i++)
  {
    prob += -0.5*pow( flux_data_sim[i] - flux_data[i], 2.0)/(err_data[i] *err_data[i]);
  }
  prob += norm_prob;

  return prob;
}

/*
 *========================================================================
 * likelihood probability at initial step
 *
 *========================================================================
 */
double prob_initial_recon(const void *model)
{
  double prob=0.0;
  int i;

  which_particle_update = dnest_get_which_particle_update();

  genlc_array_initial(model);
  
  gsl_interp_init(gsl_linear_sim, time_sim, flux_sim, nd_sim);
  
  for(i=0; i<ndata; i++)
  {
    flux_data_sim[i] = gsl_interp_eval(gsl_linear_sim, time_sim, flux_sim, time_data[i], gsl_acc_sim);
  }

  for(i=0; i<ndata; i++)
  {
    prob += -0.5*pow( flux_data_sim[i] - flux_data[i], 2.0)/(err_data[i] *err_data[i]);
  }
  prob += norm_prob;

  return prob;
}

void from_prior_recon(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params_psd_tot+1; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
  }

  for(i=1; i<nd_sim; i++)
    pm[i+num_params_psd_tot] = dnest_randn();

  for(i=num_params_psd_tot+nd_sim; i<num_params; i++)
    pm[i] = dnest_rand();

  for(i=0; i<num_params; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
}

void print_particle_recon(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

void print_particle_recon_saveoutput(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params_psd_tot; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

void read_particle_recon(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "Error: Cannot read sample file.\n");
      exit(0);
    }
  }

  return;
}

void read_particle_recon_saveoutput(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < num_params_psd_tot; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "Error: Cannot read sample file.\n");
      exit(0);
    }
  }

  return;
}

double log_likelihoods_cal_recon(const void *model)
{
  double logL;
  logL = prob_recon(model);
  return logL;
}

double log_likelihoods_cal_initial_recon(const void *model)
{
  double logL;
  logL = prob_initial_recon(model);
  return logL;
}

double log_likelihoods_cal_recon_exam(const void *model)
{
  return 0.0;
}

double perturb_recon(void *model)
{
  double logH=0.0, width, rnd;
  int which;
  double *pm = (double *)model;

  /* sample variability parameters more frequently */
  do
  {
    rnd = dnest_rand();
    if(rnd < 0.1)
      which = dnest_rand_int(num_params_psd_tot);
    else
      which = dnest_rand_int(num_recon) + num_params_psd_tot;
  }while(par_fix[which] == 1);

  which_parameter_update = which;
  
  width = ( par_range_model[which][1] - par_range_model[which][0] );
      
  if(which < num_params_psd_tot + 1)
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_psd_tot + nd_sim)
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), 0.0, 1.0);
  }

  return logH;
  
}

double perturb_recon_limits(void *model)
{
  double logH=0.0, width, rnd, limit1, limit2;
  int which, which_level;
  double *pm = (double *)model;

  /* sample variability parameters more frequently */
  do
  {
    rnd = dnest_rand();
    if(rnd < 0.1)
      which = dnest_rand_int(num_params_psd_tot);
    else
      which = dnest_rand_int(num_recon) + num_params_psd_tot;
  }while(par_fix[which] == 1);

  which_parameter_update = which;
  
  which_level = dnest_get_which_level_update();

  if( which_level > 0 )
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    width = ( par_range_model[which][1] - par_range_model[which][0] );
  }
      
  if(which < num_params_psd_tot + 1)
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_psd_tot + nd_sim)
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), 0.0, 1.0);
  }

  return logH;
  
}

void restart_action_recon(int iflag)
{
  return;
}


int get_line_number(char *fname)
{
  char buf[200], buf1[200], buf2[200], buf3[200];
  int nc;
  FILE *fp;

  fp = fopen(fname, "r");
  if(fp==NULL)
  {
    printf("Cannot open file %s.\n", fname);
    exit(0);
  }
  nc = 0;
  while(1)
  {
    fgets(buf, 200, fp);
    //make sure at least three columns each line
    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3)<3)
      break;
    if(feof(fp)!=0)
      break;
    nc++;
  }
  fclose(fp);
  return nc;
}

/*!
 * get file name of posterior sample. 
 */
void get_posterior_sample_file(char *fname, char *samplefile)
{
  FILE *fp;
  char buf[200], buf1[200];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, 200, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  fgets(buf, 200, fp);
//  sscanf(buf, "%d", &options.new_level_interval);

  fgets(buf, 200, fp);
//  sscanf(buf, "%d", &options.save_interval);

  fgets(buf, 200, fp);
//  sscanf(buf, "%d", &options.thread_steps);

  fgets(buf, 200, fp);
//  sscanf(buf, "%d", &options.max_num_levels);

  fgets(buf, 200, fp);
//  sscanf(buf, "%lf", &options.lambda);

  fgets(buf, 200, fp);
//  sscanf(buf, "%lf", &options.beta);

  fgets(buf, 200, fp);
//  sscanf(buf, "%d", &options.max_num_saves);

  fgets(buf, 200, fp);
//  sscanf(buf, "%s", options.sample_file);
  
  fgets(buf, 200, fp);
//  sscanf(buf, "%s", options.sample_info_file);
  
  fgets(buf, 200, fp);
//  sscanf(buf, "%s", options.levels_file);
  
  fgets(buf, 200, fp);
//  sscanf(buf, "%s", options.sampler_state_file);
  
  fgets(buf, 200, fp);
  sscanf(buf, "%s", samplefile);
  fclose(fp);
}


void set_par_fix()
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

/*!
 * this function set the parameter range.
 */
void set_par_range()
{
  int i=0, j;
  
  if(parset.psd_model_enum == simple )
  {
    i=0;
    var_range_model[i][0] = log(1.0e-10); // A
    var_range_model[i++][1] = log(1.0e6);
  
    switch(parset.psd_type)
    {
      case 0:   // single power-law
        var_range_model[i][0] = 0.0;  //slope
        var_range_model[i++][1] = 5.0;
        break;

      case 1:   // damped random walk
        var_range_model[i][0] = log(freq_limit_data_lower/(2.0*PI)); //characteristic frequency
        var_range_model[i++][1] = log(freq_limit_data_upper/(2.0*PI));
        break;
  
      case 2:  // bending power-law
        var_range_model[i][0] = 1.0; //alpha_hi
        var_range_model[i++][1] = 5.0;
 
        var_range_model[i][0] = 0.0; //alpha_hi-alpha_lo
        var_range_model[i++][1] = 4.0;

        var_range_model[i][0] = log(freq_limit_data_lower); //the smallest freq as bending frequency
        var_range_model[i++][1] = log(freq_limit_data_upper); // the largest freq
        break;

      default:
        var_range_model[i][0] = 0.0;
        var_range_model[i++][1] = 5.0;
    }
  }
  else if(parset.psd_model_enum == harmonic)
  {
    i=0;
    var_range_model[i][0] = -10.0;
    var_range_model[i++][1] =  10.0;

    var_range_model[i][0] = fmin(-10.0, log(freq_limit_data_lower));
    var_range_model[i++][1] = fmax(10.0, log(freq_limit_data_upper));

    for(j=1; j<parset.harmonic_term_num; j++)
    {
      var_range_model[i][0] = -10.0;
      var_range_model[i++][1] =  10.0;

      var_range_model[i][0] = log(freq_limit_data_lower) 
           + (j-1) * (log(freq_limit_data_upper) - log(freq_limit_data_lower))/(parset.harmonic_term_num-1);
      var_range_model[i++][1] =  log(freq_limit_data_lower) 
           + (j) * (log(freq_limit_data_upper) - log(freq_limit_data_lower))/(parset.harmonic_term_num-1);

      var_range_model[i][0] = log(1.0001);
      var_range_model[i++][1] =  log(200.0);
    }

  }
  else if(parset.psd_model_enum == carma)
  {
    i=0;
    
    var_range_model[i][0] = -10.0; //sigma
    var_range_model[i++][1] =  10.0;

    for(j=0; j<parset.carma_p/2; j++)
    {
      var_range_model[i][0] = log(freq_limit_data_lower/(2.0*PI));  // Lorentzian width
      var_range_model[i++][1] =  log(freq_limit_data_upper/(2.0*PI));

      //Lorentz centriod
      var_range_model[i][0] = log(freq_limit_data_lower) 
           + (j) * (log(freq_limit_data_upper) - log(freq_limit_data_lower))/((int)(parset.carma_p/2));
      var_range_model[i++][1] =  log(freq_limit_data_lower) 
           + (j+1) * (log(freq_limit_data_upper) - log(freq_limit_data_lower))/((int)(parset.carma_p/2));
    }
    if(parset.carma_p%2 == 1)
    {
      var_range_model[i][0] = -10.0;
      var_range_model[i++][1] =  10.0;
    }

    // mean-averaging coefficients
    for(j=0; j<parset.carma_q; j++)
    {
      var_range_model[i][0] = -10.0;
      var_range_model[i++][1] =  log(1000.0);
    }
    
  }
  else
  {

  }
  
  var_range_model[i][0] = log(1.0e-10); //noise
  var_range_model[i++][1] = log(1.0e3);

  if(parset.psdperiod_enum > delta)
  {
    if(noise_power < flux_var)
    {
      var_range_model[i][0] = log(noise_power);  //Ap
      var_range_model[i++][1] = log(flux_var);
    }
    else
    {
      var_range_model[i][0] = log(noise_power*0.001);  //Ap
      var_range_model[i++][1] = log(noise_power*10.0);
    }

    var_range_model[i][0] = log(freq_limit_data_lower); //center
    var_range_model[i++][1] = log(freq_limit_data_upper);

    var_range_model[i][0] = log(freq_limit_data_lower*0.001);   //sigma
    var_range_model[i++][1] = log(freq_limit_data_upper);
  }
  else if(parset.psdperiod_enum > none)
  {
    var_range_model[i][0] = log(1.0e-10);  //Ap
    var_range_model[i++][1] = log(1.0e6);

    var_range_model[i][0] = log(freq_limit_data_lower); //center
    var_range_model[i++][1] = log(freq_limit_data_upper);

    var_range_model[i][0] = 0.0;   //sigma
    var_range_model[i++][1] = 1.0;
  }
  
  var_range_model[i][0] = -10.0;  //zero-frequency power
  var_range_model[i++][1] = 10.0;

  var_range_model[i][0] = -5.0;  //frequency series
  var_range_model[i++][1] = 5.0;

  // variability parameters
  for(i=0; i<num_params_psd_tot+1; i++)
  {
    par_range_model[i][0] = var_range_model[i][0];
    par_range_model[i][1] = var_range_model[i][1];
  }

  // continuum light curve parameters
  for(i=num_params_psd_tot+1; i<num_params_psd_tot+nd_sim; i++)
  {
    par_range_model[i][0] = var_range_model[num_params_psd_tot+1][0];
    par_range_model[i][1] = var_range_model[num_params_psd_tot+1][1];
  }

  for(i=num_params_psd_tot+nd_sim; i< num_params; i++)
  {
    par_range_model[i][0] = 0.0;
    par_range_model[i][1] = 1.0;
  }

  return;
}

/*
 *  set psd functions
 */
void set_psd_functions()
{
  //initialize periodic PSD
  switch(parset.psdperiod_enum)
  {
    case none:
      psdfunc_period_array = NULL;
      psdfunc_period_sqrt_array = NULL;
      func_genlc_array = genlc_array;
      parset.num_params_psdperiod = 0;
      break;

    case delta:
      psdfunc_period_array = NULL;
      psdfunc_period_sqrt_array = NULL;
      func_genlc_array = genlc_array_simple_delta;
      parset.num_params_psdperiod = 3;
      break;

    case gaussian:
      psdfunc_period_array = psd_gaussian_array;
      psdfunc_period_sqrt_array = psd_gaussian_sqrt_array;
      func_genlc_array = genlc_array_simple_period;
      parset.num_params_psdperiod = 3;
      break;

    case lorentzian:
      psdfunc_period_array = psd_lorentz_array;
      psdfunc_period_sqrt_array = psd_lorentz_sqrt_array;
      func_genlc_array = genlc_array_simple_period;
      parset.num_params_psdperiod = 3; 
      break;

    default:
      psdfunc_period_array = NULL;
      psdfunc_period_sqrt_array = NULL;
      func_genlc_array = genlc_array;
      parset.num_params_psdperiod = 0;
      break;
  }

  //initalize PSD functions and number of PSD parameters
  if(parset.psd_model_enum == simple)
  {
    switch(parset.psd_type)
    {
      case 0: // single power-law
        psdfunc_array = psd_power_law_array;
        psdfunc_sqrt_array = psd_power_law_sqrt_array;
        parset.num_params_psd = 3;
        break;
  
      case 1: // damped random walk
        psdfunc_array = psd_drw_array;
        psdfunc_sqrt_array = psd_drw_sqrt_array;
        parset.num_params_psd = 3;
        break;
      
      case 2: // bending power-law
        psdfunc_array = psd_bending_power_law_array;
        psdfunc_sqrt_array = psd_bending_power_law_sqrt_array;   
        parset.num_params_psd = 5;
        break;

      default:
        psdfunc_array = psd_power_law_array;
        psdfunc_sqrt_array = psd_power_law_sqrt_array;
        parset.num_params_psd = 3;
        break;
    }
  }
  else if(parset.psd_model_enum == harmonic)
  {
    psdfunc = psd_harmonic;
    psdfunc_sqrt = psd_harmonic_sqrt; 

    psdfunc_array = psd_harmonic_array;
    psdfunc_sqrt_array = psd_harmonic_sqrt_array;   
 
    parset.num_params_psd = 2+3*(parset.harmonic_term_num-1) + 1; //including noise
  }
  else if(parset.psd_model_enum == carma)
  {
    psdfunc_array = psd_carma_array;
    psdfunc_sqrt_array = psd_carma_sqrt_array;
    parset.num_params_psd = parset.carma_p + parset.carma_q + 1 + 1; 
  }
  else
  {

  }

  return;
}

void time_cad_cal()
{
  int i;
  double *cad;
  cad = malloc((ndata -1)*sizeof(double));

  for(i=0; i<ndata-1; i++)
    cad[i] = time_data[i+1] - time_data[i];
  
  qsort(cad, ndata-1, sizeof(double), recon_cmp);

  if(ndata<2)
    time_cad_media = cad[0];
  else
    time_cad_media = cad[ndata/2-1];

  free(cad);

  return;
}

int recon_cmp(const void *a, const void *b)
{
  return *(double *)a > *(double *)b?1:-1;
}

void sim()
{
  FILE *fp;
  int i, i1, i2;
  void *model;
  double *pm;
  char fname[200];

  const gsl_rng_type * gsl_T;
  gsl_rng * gsl_r;
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);

  if(recon_flag_seed == 1)
  {
    gsl_rng_set(gsl_r, recon_seed+thistask*10);
  }
  else
  {
    recon_seed = time(NULL);
    gsl_rng_set(gsl_r, recon_seed+thistask*10);
    printf("# random seed: %d\n", recon_seed);
  }
  

  model = malloc(num_params * sizeof(double));
  
  pm = (double *)model;

  memcpy(pm, parset.psd_arg, num_params_psd_tot * sizeof(double));

  pm[num_params_psd_tot] = 0.0;
  
  for(i=1; i<nd_sim; i++)
  {
    pm[num_params_psd_tot + i] = gsl_ran_ugaussian(gsl_r);
  }

  for(i=num_params_psd_tot+nd_sim; i<num_params; i++)
  {
    pm[i] = gsl_rng_uniform(gsl_r);
  }
  
  which_parameter_update = -1;
  which_particle_update = 0;
  genlc_array_initial(model);

  // add Gaussian noise
  for(i=0; i<nd_sim; i++)
  {
    flux_sim[i] += parset.ferr*gsl_ran_ugaussian(gsl_r);
  }
  
  sprintf(fname, "%s/%s", parset.file_dir, parset.file_sim);
  fp = fopen(fname, "w");
  if(fp==NULL)
  {
    printf("Cannot open file fname.\n");
    exit(0);
  }
  
  i1 = nd_sim/2 - nd_sim/V/2;
  i2 = nd_sim/2 + nd_sim/V/2;
  for(i=i1; i<i2; i=i+W)
  {
    if(gsl_rng_uniform(gsl_r) > parset.fbad)
    {
      fprintf(fp, "%f %f %f\n", time_sim[i], flux_sim[i], parset.ferr);
    }
  }
  fclose(fp);

  sprintf(fname, "%s/%s_all", parset.file_dir, parset.file_sim);
  fp = fopen(fname, "w");
  if(fp==NULL)
  {
    printf("Cannot open file fname.\n");
    exit(0);
  }
  
  i1 = nd_sim/2 - nd_sim/V/2;
  i2 = nd_sim/2 + nd_sim/V/2;
  for(i=i1; i<i2; i=i+W)
  {
    fprintf(fp, "%f %f %f\n", time_sim[i], flux_sim[i], parset.ferr);
  }
  fclose(fp);

  sprintf(fname, "%s/%s_full", parset.file_dir, parset.file_sim);
  fp = fopen(fname, "w");
  if(fp==NULL)
  {
    printf("Cannot open file fname.\n");
    exit(0);
  }
  
  for(i=0; i<nd_sim; i++)
  {
    fprintf(fp, "%f %f %f\n", time_sim[i], flux_sim[i], parset.ferr);
  }
  fclose(fp);
  
  free(model);
  gsl_rng_free(gsl_r);

  return;

}



void test()
{
  int i;
  void *model;
  double *pm;
  
  const gsl_rng_type * gsl_T;
  gsl_rng * gsl_r;
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);
  if(recon_flag_seed == 1)
  {
    gsl_rng_set(gsl_r, recon_seed+thistask*10);
  }
  else
  {
    gsl_rng_set(gsl_r, time(NULL)+thistask*10);
  }

  model = malloc(num_params * sizeof(double));
  
  pm = (double *)model;

  pm[0] = log(1.0e2);
  pm[1] = 2.0;
  pm[2] = log(1.0e-50);
  pm[3] = 1.0;
  
  for(i=1; i<num_recon; i++)
  {
    pm[num_params_psd_tot + i] = gsl_ran_ugaussian(gsl_r);
  }

  genlc_array_initial(model);
  
  for(i=0; i<nd_sim; i++)
  {
    printf("%f %f\n", time_sim[i], flux_sim[i]);
  }

  free(model);
  gsl_rng_free(gsl_r);
  return;
}
