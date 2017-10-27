#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <fftw3.h>
#include <gsl/gsl_interp.h>

#include "allvars.h"

int thistask, totaltask;
int roottask;
int namelen;
char proc_name[MPI_MAX_PROCESSOR_NAME];

int ndata;
double *flux_data, *time_data, *err_data;
double time_media, time_cad_min, flux_data_min, flux_data_max, flux_scale, flux_mean;

int nd_sim;
double DT, V, W;
double *time_sim, *flux_sim, *flux_data_sim;
double *flux_sim_mean, *err_sim_mean;

gsl_interp_accel *gsl_acc_sim;
gsl_interp  *gsl_linear_sim;

fftw_plan pfft;
fftw_complex *fft_work;

int num_params, num_params_psd, num_recon;
double **par_range_model, **var_range_model, *par_fix_val;
int *par_fix;

int which_particle_update, which_parameter_update;
int which_level_update, num_particles;
unsigned long long int which_mcmc_steps;
int *perturb_accept;


int flag_postprc; 
double temperature;
int flag_restart;
int flag_sample_info;
int flag_temp;
int flag_sim;
int flag_help;
int flag_end;

PARSET parset;