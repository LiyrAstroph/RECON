
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <float.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <dnestvars.h>

#include "allvars.h"
#include "proto.h"

int recon()
{
  int i, argc=0;
  char **argv;
  
  strcpy(options_file, "OPTIONS");

  argv = malloc(8*sizeof(char *));
  for(i=0; i<8; i++)
  {
    argv[i] = malloc(200*sizeof(char));
  }
  
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc++], "data/restart_dnest.txt");

  if(flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
  }
  if(flag_temp == 1)
  {
    sprintf(argv[argc++], "-t%f", temperature);
  }
  if(flag_sample_info == 1)
  {
    strcpy(argv[argc++], "-c");
  }

  recon_init();

  if(flag_sim == 1)
  {
    if(thistask == roottask)
    {
      sim();
    }

  }
  else
  {
    dnest(argc, argv);
    recon_postproc();
  }
  
  recon_end();
  return 0;
}


int recon_postproc()
{
  if(thistask == roottask)
  {
    int num_ps, i, j;
    char fname[200];
    void *posterior_sample, *post_model;
    FILE *fp, *fcon, *fcon_all, *fcon_mean;

    char posterior_sample_file[200];

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
      
      genlc(post_model);

      gsl_interp_init(gsl_linear_sim, time_sim, flux_sim, nd_sim);
  
      for(j=0; j<nd_sim; j++)
      {
        flux_sim_mean[j] += flux_sim[j];
        err_sim_mean[j] += flux_sim[j]*flux_sim[j];
      }

      for(j=0; j<ndata; j++)
      {
        flux_data_sim[j] = gsl_interp_eval(gsl_linear_sim, time_sim, flux_sim, time_data[j], gsl_acc_sim);

        fprintf(fcon, "%f  %f\n", time_data[j], flux_data_sim[j]*flux_scale+flux_mean);
      }
      fprintf(fcon, "\n");

      for(j=0; j<nd_sim; j++)
      {
        fprintf(fcon_all, "%f  %f\n", time_sim[j], flux_sim[j]*flux_scale + flux_mean);
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
      fprintf(fcon_mean, "%f %f %f\n", time_sim[j], flux_sim_mean[j]*flux_scale+flux_mean, err_sim_mean[j]*flux_scale);
    }

    fclose(fp);
    fclose(fcon);
    fclose(fcon_all);
    fclose(fcon_mean);
  }

  return 0;
}

int recon_init()
{
  char fname[200];
  int i;
  double Tall, Tmin, Tmax;

  /* setup functions used for dnest*/
  from_prior = from_prior_recon;
  log_likelihoods_cal = log_likelihoods_cal_recon;
  log_likelihoods_cal_initial = log_likelihoods_cal_initial_recon;
  log_likelihoods_cal_restart = log_likelihoods_cal_restart_recon;
  perturb = perturb_recon;
  print_particle = print_particle_recon;
  copy_model = copy_model_recon;
  create_model = create_model_recon;
  get_num_params = get_num_params_recon;
  restart_clouds = restart_clouds_recon;

  roottask = 0;

  psdfunc = psd_power_law;

  strcpy(fname, "data/sim.txt");

  if(thistask == roottask)
  {
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
  }

  MPI_Bcast(time_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(flux_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(err_data, ndata, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

  time_media = (time_data[0] + time_data[ndata-1])/2.0;
  time_cad_min = time_media;
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

  /* fft */
  if(flag_sim==1)
  {
    time_media = 0.0;
    DT = 1.0;
    nd_sim = 1000;
  }
  else
  {
    V = 2.0;
    W = 2.0;
    Tall = time_data[ndata-1] - time_data[0];
    Tmin = time_data[0] - 0.5*(V-1.0)*Tall;
    Tmax = time_data[ndata-1] + 0.5*(V-1.0)*Tall;
    Tall = Tmax - Tmin;
  
    DT = (time_data[ndata-1] - time_data[0])/(ndata -1)/W;
    nd_sim = Tall/DT + 1;
    nd_sim = (nd_sim/2) * 2; //make sure it is an even number.
  }
  
  num_recon = nd_sim;
  num_params_psd = 3;

  num_params = num_recon + num_params_psd;
  size_of_modeltype = num_params * sizeof(double);

  var_range_model = malloc((num_params_psd+2)*sizeof(double *));
  for(i=0; i<num_params_psd+2; i++)
  {
    var_range_model[i] = malloc(2*sizeof(double));
  }
  i=0;
  var_range_model[i][0] = log(1.0e-10);
  var_range_model[i++][1] = log(1.0e6);

  var_range_model[i][0] = 0.0;
  var_range_model[i++][1] = 5.0;

  var_range_model[i][0] = log(1.0e-10);
  var_range_model[i++][1] = log(1.0e3);

  var_range_model[i][0] = -100.0;
  var_range_model[i++][1] = 100.0;

  var_range_model[i][0] = -10.0;
  var_range_model[i++][1] = 10.0;

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));

  set_par_range();
  //setup fixed parameters
  for(i=0; i<num_params; i++)
    par_fix[i] = 0;

  time_sim = malloc(nd_sim * sizeof(double));
  flux_sim = malloc(nd_sim * sizeof(double));
  flux_sim_mean = malloc(nd_sim * sizeof(double));
  err_sim_mean = malloc(nd_sim * sizeof(double));
  
  gsl_acc_sim = gsl_interp_accel_alloc();
  gsl_linear_sim = gsl_interp_alloc(gsl_interp_linear, nd_sim);

  fft_work = (fftw_complex *)fftw_malloc( nd_sim * sizeof(fftw_complex));

  pfft = fftw_plan_dft_c2r_1d(nd_sim, fft_work, flux_sim, FFTW_MEASURE);

  if(thistask == roottask)
  {
    get_num_particles(options_file);
  }
  MPI_Bcast(&num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  perturb_accept = malloc(num_particles * sizeof(int));

  return 0;
}

int recon_end()
{
  int i;

  fftw_destroy_plan(pfft);

  gsl_interp_free(gsl_linear_sim);
  gsl_interp_accel_free(gsl_acc_sim);

  fftw_free(fft_work);
  free(perturb_accept);
  
  for(i=0; i<num_params_psd+2; i++)
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
  free(time_data);
  free(flux_data);
  free(flux_data_sim);

  return 0;
}

int genlc(const void *model)
{
  int i;
  double *arg, freq;
  double *pm = (double *)model;

  arg = malloc(num_params_psd * sizeof(double));
  memcpy(arg, pm, num_params_psd * sizeof(double));

  fft_work[0][0] = pm[num_params_psd+0]; //zero-frequency power.
  fft_work[0][1] = 0.0;

  for(i=0; i<nd_sim/2-1; i++)
  {
    fft_work[i+1][0] = pm[num_params_psd+1+2*i];
    fft_work[i+1][1] = pm[num_params_psd+1+2*i+1];
  }
  fft_work[nd_sim/2][0] = pm[num_params_psd + nd_sim-1];
  fft_work[nd_sim/2][1] = 0.0;

  for(i=1; i<nd_sim/2+1; i++)
  {
    freq = i*1.0/(nd_sim * DT);
    fft_work[i][0] *= sqrt(psdfunc(freq, arg)/2.0);
    fft_work[i][1] *= sqrt(psdfunc(freq, arg)/2.0);
  }

  fftw_execute(pfft);

  for(i=0; i<nd_sim; i++)
  {
    //printf("%f\n", flux_workspace[i]);
    time_sim[i] = (i - nd_sim/2) * DT  + time_media;
    flux_sim[i] = flux_sim[i]/sqrt(nd_sim);
  }

  return 0;
}

double prob_recon(const void *model)
{
  double prob;
  int i;

  genlc(model);

  gsl_interp_init(gsl_linear_sim, time_sim, flux_sim, nd_sim);
  
  for(i=0; i<ndata; i++)
  {
    flux_data_sim[i] = gsl_interp_eval(gsl_linear_sim, time_sim, flux_sim, time_data[i], gsl_acc_sim);
  }

  prob = 0.0;
  for(i=0; i<ndata; i++)
  {
    prob += -0.5*pow( flux_data_sim[i] - flux_data[i], 2.0)/err_data[i]/err_data[i] 
            -0.5*log(2.0*PI) - log(err_data[i]);
  }
  return prob;
}

void from_prior_recon(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params_psd+1; i++)
  {
    pm[i] = var_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
  }

  for(i=1; i<num_recon; i++)
    pm[i+num_params_psd] = dnest_randn();

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
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
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
  logL = prob_recon(model);
  return logL;
}

double log_likelihoods_cal_restart_recon(const void *model)
{
  double logL;
  logL = prob_recon(model);
  return logL;
}

double perturb_recon(void *model)
{
  double logH=0.0, width, rnd, limit1, limit2;
  int which;
  double *pm = (double *)model;

  /* sample variability parameters more frequently */
  do
  {
    rnd = dnest_rand();
    if(rnd < 0.1)
      which = dnest_rand_int(num_params_psd);
    else
      which = dnest_rand_int(num_recon) + num_params_psd;
  }while(par_fix[which] == 1);

  which_parameter_update = which;

  /* level-dependent width */
  which_level_update = which_level_update > (size_levels - 10)?(size_levels-10):which_level_update;
  which_level_update = which_level_update <0?0:which_level_update;
  if( which_level_update != 0)
  {
    limit1 = limits[(which_level_update-1) * num_params *2 + which *2];
    limit2 = limits[(which_level_update-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    width = ( par_range_model[which][1] - par_range_model[which][0] );
  }

  //width = ( par_range_model[which][1] - par_range_model[which][0] );
  
  if(which < num_params_psd + 1)
  {
    pm[which] += dnest_randh() * width;
    wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }

  return logH;
  
}

void copy_model_recon(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void* create_model_recon()
{
  return (void *)malloc( size_of_modeltype );
}

int get_num_params_recon()
{
  return num_params;
}

void restart_clouds_recon(int iflag)
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
 * get number of particles from the option file.
 */
void get_num_particles(char *fname)
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
  sscanf(buf, "%d", &num_particles);
  fclose(fp);
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

int read_data(char *fname, int n, double *t, double *f, double *e)
{
  FILE *fp;
  int i;
  char buf[200];

  fp = fopen(fname, "r");
  if(fp==NULL)
  {
    printf("Cannot open file %s.\n", fname);
    exit(0);
  }

  for(i=0; i<n; i++)
  {
    fgets(buf, 200, fp);
    sscanf(buf, "%lf %lf %lf\n", &t[i], &f[i], &e[i]);
  }
  fclose(fp);

  return 0;
}

/*!
 * this function set the parameter range.
 */
void set_par_range()
{
  int i;

  // variability parameters

  for(i=0; i<num_params_psd; i++)
  {
    par_range_model[i][0] = var_range_model[i][0];
    par_range_model[i][1] = var_range_model[i][1];
  }

  // continuum light curve parameters
  for(i=num_params_psd; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[num_params_psd][0];
    par_range_model[i][1] = var_range_model[num_params_psd][1];
  }
  return;
}

double psd_drw(double fk, double *arg)
{
  double A=exp(arg[0]), fknee=exp(arg[1]), cnoise = exp(arg[2]);

  return A/(1.0 + pow(fk/fknee, 2.0)) + cnoise;
}

double psd_power_law(double fk, double *arg)
{
  double A=exp(arg[0]), alpha=arg[1], cnoise=exp(arg[2]);

  return A * pow(fk, -alpha) + cnoise;
}

void sim()
{
  FILE *fp;
  int i;
  void *model;
  double *pm;
  
  const gsl_rng_type * gsl_T;
  gsl_rng * gsl_r;
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);
  gsl_rng_set(gsl_r, time(NULL));

  model = malloc(size_of_modeltype);
  
  pm = (double *)model;

  pm[0] = log(1.0);
  pm[1] = 1.5;
  pm[2] = log(1.0e-50);
  pm[3] = 0.0;
  
  for(i=1; i<num_recon; i++)
  {
    pm[num_params_psd + i] = gsl_ran_gaussian(gsl_r, 1.0);
  }

  genlc(model);
  
  fp = fopen("./data/sim.txt", "w");
  if(fp==NULL)
  {
    printf("Cannot open file sim.txt.\n");
    exit(0);
  }

  for(i=0; i<nd_sim; i+=(int)(nd_sim/100))
  {
    fprintf(fp, "%f %f %f\n", time_sim[i], flux_sim[i], 0.1);
  }
  
  fclose(fp);
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
  gsl_rng_set(gsl_r, time(NULL));

  model = malloc(size_of_modeltype);
  
  pm = (double *)model;

  pm[0] = log(1.0e2);
  pm[1] = 2.0;
  pm[2] = log(1.0e-50);
  pm[3] = 1.0;
  
  for(i=1; i<num_recon; i++)
  {
    pm[num_params_psd + i] = gsl_ran_gaussian(gsl_r, 1.0);
  }

  genlc(model);
  
  for(i=0; i<nd_sim; i++)
  {
    printf("%f %f\n", time_sim[i], flux_sim[i]);
  }

}