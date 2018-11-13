/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

#ifndef _ALLVARS_H
#define _ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

#define PI  (M_PI)

extern int roottask, thistask, totaltask;
extern int namelen;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];

extern int ndata;
extern double *flux_data, *time_data, *err_data;
extern double time_media, time_cad_min, time_cad_media, flux_data_min, flux_data_max, flux_scale, flux_mean;

extern double freq_limit_data_lower, freq_limit_data_upper;

extern int nd_sim;
extern double DT, V, W;
extern double *time_sim, *flux_sim, *flux_data_sim;
extern double *flux_sim_mean, *err_sim_mean;
extern double *freq_array;
extern complex **freq_array_pow;
extern int idx_limit;
extern double norm_psd, norm_prob;
extern double *workspace_psd;
extern double **workspace_genlc, **workspace_genlc_perturb;
extern double **workspace_genlc_period, **workspace_genlc_period_perturb;
extern int *which_parameter_update_prev;


extern gsl_interp_accel *gsl_acc_sim;
extern gsl_interp  *gsl_linear_sim;
extern fftw_plan pfft;
extern fftw_complex *fft_work;

extern int num_params, num_params_psd_tot, num_recon;
extern double **par_range_model, **var_range_model, *par_fix_val;
extern int *par_fix;

extern int which_level_update, which_particle_update, which_parameter_update, num_particles;

extern int recon_flag_postprc; 
extern double recon_temperature;
extern int recon_flag_restart;
extern int recon_flag_sample_info;
extern int recon_flag_temp;
extern int recon_flag_sim;
extern int recon_flag_help;
extern int recon_flag_end;
extern int recon_flag_limits;
extern int recon_flag_prior_exam;
extern int recon_flag_seed;
extern int recon_flag_cal_psd;
extern int recon_seed;

typedef struct
{
  char file_dir[256]; 
  char param_file[256]; 
  char file_name[256];
  char psd_model[32];
  char psdperiod_model[32];
  char str_psd_arg[256];
  char file_sim[256];
  char psd_type_str[32];

  double V, W;
  double freq_limit;

  int psd_type, harmonic_term_num, carma_p, carma_q;
  double slope_endmatch;

  int psd_model_enum, psdperiod_enum;
  int num_params_psd, num_params_psdperiod;
  double fbad, ferr;

  int nd_sim;
  double DT;
  double psd_arg[32];  // maximal number of arguments is 20.

  int flag_domain;
  int flag_endmatch;
  int flag_whitenoise;
  int flag_saveoutput;
}PARSET;
extern PARSET parset;

extern FILE *finfo;

enum PSDMODEL {simple, harmonic, carma, celerite};
enum PSDPERIODMODEL {none, gaussian, lorentz};

// for harmonic

// for carma
extern complex *workspace_complex;
#endif