
#ifndef _ALLVARS_H
#define _ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

#define PI  (M_PI)

extern int thistask, totaltask;
extern int roottask;
extern int namelen;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];

extern int ndata;
extern double *flux_data, *time_data, *err_data;
extern double time_media, time_cad_min, flux_data_min, flux_data_max, flux_scale, flux_mean;

extern int nd_sim;
extern double DT, V, W;
extern double *time_sim, *flux_sim, *flux_data_sim;
extern double *flux_sim_mean, *err_sim_mean;

extern gsl_interp_accel *gsl_acc_sim;
extern gsl_interp  *gsl_linear_sim;
extern fftw_plan pfft;
extern fftw_complex *fft_work;

extern int num_params, num_params_psd, num_recon;
extern double **par_range_model, **var_range_model, *par_fix_val;
extern int *par_fix;


/* size of model type, defined in dnest */
extern int size_of_modeltype;


extern int which_particle_update, which_parameter_update;
extern int which_level_update, num_particles;
extern unsigned long long int which_mcmc_steps;
extern int *perturb_accept;
extern double *limits;

extern int flag_postprc; 
extern double temperature;
extern int flag_restart;
extern int flag_sample_info;
extern int flag_temp;
extern int flag_sim;
extern int flag_help;

#endif