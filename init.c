/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

/*! \file init.c
 *  \brief initiats the running.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "dnestvars.h"

#include "allvars.h"
#include "proto.h"

/*
 * this function initiates the running.
 * 
 */ 
void init()
{
  char fname[200];
  int i;
  double Tmax, Tmin, Tall;
  
  /* set the root task */
  roottask = 0;

  /* read parameter set */
  read_parset();
  /* set PSD functions */
  set_psd_functions();

  /* total number of PSD parameters */
  num_params_psd_tot = parset.num_params_psd + parset.num_params_psdperiod;
  
  /* open file for output recon configurations */
  if(thistask == roottask)
  {
    sprintf(fname, "%s/data/recon_info.txt", parset.file_dir);
    finfo = fopen(fname, "w");
    if(finfo == NULL)
    {
      printf("Cannot open file %s.\n", fname);
      exit(0);
    }
  }

  /* fft related configuration */
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
    /* get number of data points */
    if(thistask == roottask)
    {
      sprintf(fname, "%s/%s", parset.file_dir, parset.file_name);
      ndata = get_line_number(fname);
    }
    MPI_Bcast(&ndata, 1, MPI_INT, roottask, MPI_COMM_WORLD);
 
    /* allocate data array momery */
    time_data = malloc(ndata * sizeof(double));
    flux_data = malloc(ndata * sizeof(double));
    err_data = malloc(ndata * sizeof(double));
    flux_data_sim = malloc(ndata * sizeof(double));
    
    freq_data = malloc(ndata * sizeof(double));
    psd_data = malloc(ndata * sizeof(double));

    /* read data, calculate periodogram, and calculate cadence */
    if(thistask == roottask)
    {
      read_data(fname, ndata, time_data, flux_data, err_data);
      psddata_cal();
      time_cad_cal();
    }

    MPI_Bcast(time_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(flux_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(err_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

    MPI_Bcast(freq_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(psd_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(&nf_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);
    
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
      /* subtract the mean, and then scale */
      flux_data[i] = (flux_data[i] - flux_mean)/flux_scale;
      err_data[i] /= flux_scale;
    }

    /* scale periodogram data */
    for(i=0; i<nf_data; i++)
    {
      psd_data[i] /= (flux_scale * flux_scale);
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

    par_fit_best = malloc(num_params_psd_tot * sizeof(double));
    par_fit_best_std = malloc(num_params_psd_tot * sizeof(double));
  }
  if(thistask == roottask)
  {
    fclose(finfo);
  }

  var_range_model = malloc((num_params_psd_tot+2)*sizeof(double *));
  for(i=0; i<num_params_psd_tot+2; i++)
  {
    var_range_model[i] = malloc(2*sizeof(double));
  }

  set_var_range();

  return;
}

void end()
{
  int i;
  for(i=0; i<num_params_psd_tot+2; i++)
  {
    free(var_range_model[i]);
  }
  free(var_range_model);

  if(recon_flag_sim != 1)
  {
    free(time_data);
    free(flux_data);
    free(err_data);
    free(flux_data_sim);
    free(freq_data);
    free(psd_data);
    free(par_fit_best);
    free(par_fit_best_std);
  }

  return;
}

void set_var_range()
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
  
  var_range_model[i][0] = log(1.0e-6); //noise
  var_range_model[i++][1] = log(1.0e3);

  if(parset.psdperiod_enum > delta)
  {
    if(noise_power < flux_var)
    {
      var_range_model[i][0] = log(noise_power);  //Ap
      var_range_model[i++][1] = log(flux_var*1000.0);
    }
    else
    {
      var_range_model[i][0] = log(noise_power*0.001);  //Ap
      var_range_model[i++][1] = log(noise_power*100.0);
    }

    var_range_model[i][0] = log(freq_limit_data_lower); //center
    var_range_model[i++][1] = log(freq_limit_data_upper);

    var_range_model[i][0] = log(freq_limit_data_lower*0.01);   //sigma
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
