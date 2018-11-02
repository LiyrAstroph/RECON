/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_fit.h>

#include "proto.h"
#include "allvars.h"


/*!
 *  This function calculates PSD for data.
 */
int psddata_cal()
{
  FILE *fp;
  int i, nf;
  double *freq, *psd;
  char fname[200], str1[200], str2[200], *pstr;

  freq = malloc(ndata * sizeof(double));
  psd = malloc(ndata * sizeof(double));

  psd_fft(time_data, flux_data, ndata, freq, psd, &nf);

  strcpy(str1, parset.file_name);

  pstr = strchr(str1, '/');
  strcpy(str2, pstr+1);

  *pstr = '\0';

  sprintf(fname, "%s/%s/psd_%s", parset.file_dir, str1, str2);

  fp = fopen(fname, "w");
  
  if(fp == NULL)
  {
    printf("Cannot open file %s.\n", fname);
    exit(0);
  }

  for(i=0; i<nf; i++)
  {
    fprintf(fp, "%e %e\n", freq[i], psd[i]);
  }
  fclose(fp);

  psd_fit_check(freq, psd, nf);

  free(freq);
  free(psd);
  return 0;
}

int psd_fft(double *t, double *f, int n, double *freq, double *psd, int *nf)
{
  int i;
  double mean;
  double *ftemp, *ttemp;
  fftw_plan ppsd;
  fftw_complex *psd_work;

  ttemp = malloc(n*sizeof(double));
  ftemp = malloc(n*sizeof(double));
  psd_work = fftw_malloc(n*sizeof(fftw_complex));
  ppsd = fftw_plan_dft_r2c_1d(n, ftemp, psd_work, FFTW_MEASURE);

  resample(t, f, n, ttemp, ftemp);

  // subtract mean
  mean = 0.0;
  for(i=0; i<n; i++)
  {
    mean += ftemp[i];
  }
  mean /= n;

  for(i=0; i<n; i++)
  {
    ftemp[i] -= mean;
  }

  fftw_execute(ppsd);
  
  //omit the zero frequency
  *nf = n/2;
  for(i=1; i< *nf + 1; i++)
  {
    freq[i-1] = i * 1.0/(ttemp[n-1] - ttemp[0]);
    psd[i-1] = psd_work[i][0]*psd_work[i][0] + psd_work[i][1]*psd_work[i][1];

    // normalize
    psd[i-1] /= n;
    psd[i-1] = psd[i-1] * 2.0 * (ttemp[n-1] - ttemp[0])/n;
  }

  fftw_destroy_plan(ppsd);
  fftw_free(psd_work);
  free(ftemp);
  free(ttemp);
  return 0;
}


int psd_fft_rebin(double *freq, double *psd, int nf, double *freq_rebin, double *psd_rebin, int *nf_rebin)
{
  int i, i0, i1, ic, nbins_min;
  double width, width_log, f0, f1, freq0, psum, fsum;
  double *freq_work, *psd_work;
  
  freq_work = malloc(nf*sizeof(double));
  psd_work = malloc(nf*sizeof(double));

  freq0 = freq[0];
  width = freq[1] - freq[0];
  width_log = log10(1.3);
  nbins_min = 2;

  for(i=0; i<nf; i++)
  {
    freq_work[i] = log10(freq[i]);
    psd_work[i] = log10(psd[i]);
  }

  ic = 0;
  f0 = freq_work[0];
  i0 = 0;
  while(f0 < freq_work[nf-1])
  {
    f1 = f0 + width_log;
    //asume that freq is uniformly spaced
    i1 = (int)((pow(10.0, f1) - freq0)/width);
 
    //two bins at least
    if(i1 - i0 < nbins_min)
    {
      i1 = i0 + nbins_min;
    }

    // take care of the end
    if( (nf-1 - i1 < nbins_min) || (freq_work[nf-1] - f1 < width_log) )
      i1 = nf - 1;

    // update f1
    f1 = freq_work[i1];

    freq_rebin[ic] = 0.5 * (freq_work[i0] + freq_work[i1]);
    fsum = freq_work[i1] - freq_work[i0];
    psum = 0.0;
    for(i=i0; i<i1; i++)
    {
      psum += 0.5*(psd_work[i] + psd_work[i+1]) 
                 * (freq_work[i+1] - freq_work[i]);
    }

    psd_rebin[ic] = psum/fsum;

    freq_rebin[ic] = pow(10.0, freq_rebin[ic]);
    psd_rebin[ic] = pow(10.0, psd_rebin[ic]);

    i0 = i1;
    f0 = f1;
    ic += 1;
  }
  
  *nf_rebin = ic;

  free(psd_work);
  free(freq_work);
  return 0;
}

int resample(double *t, double *f, int n, double *ts, double *fs)
{
  int i;
  
  gsl_interp_accel *gsl_acc;
  gsl_interp *gsl_intp;

  gsl_acc = gsl_interp_accel_alloc();
  gsl_intp = gsl_interp_alloc(gsl_interp_linear, n);

  gsl_interp_init(gsl_intp, t, f, n);

  for(i=1; i<n-1; i++)
  {
    ts[i] = i * (t[n-1] - t[0])/(n-1.0) + t[0];
  }
  ts[0] = t[0];
  ts[n-1] = t[n-1];

  for(i=0; i<n; i++)
  {
    fs[i] = gsl_interp_eval(gsl_intp, t, f, ts[i], gsl_acc);
  }

  gsl_interp_free(gsl_intp);
  gsl_interp_accel_free(gsl_acc);

  return 0;
}

void psd_fit_check(double *freq, double *psd, int nf)
{
  FILE *fp;
  int i, nf_rebin;
  double c0, c1, cov00, cov01, cov11, res;
  double *freq_rebin, *psd_rebin;
  char fname[200], str1[200], str2[200], *pstr;
  
  freq_rebin = malloc(nf*sizeof(double));
  psd_rebin = malloc(nf*sizeof(double));

  psd_fft_rebin(freq, psd, nf, freq_rebin, psd_rebin, &nf_rebin);

  for(i=0; i<nf_rebin; i++)
  {
    freq_rebin[i] = log10(freq_rebin[i]);
    psd_rebin[i] = log10(psd_rebin[i]);
  }
  gsl_fit_linear(freq_rebin, 1, psd_rebin, 1, nf_rebin, &c0, &c1, &cov00, &cov01, &cov11, &res);

  printf("# Single power-law PSD fit, slope: %f.\n", c1);
  if( -c1 > 1.8 && -c1 < 2.2)
  {
    printf("# **** Leakage WARNING! ****\n");
  }

  strcpy(str1, parset.file_name);

  pstr = strchr(str1, '/');
  strcpy(str2, pstr+1);

  *pstr = '\0';

  sprintf(fname, "%s/%s/psd_%s_rebin", parset.file_dir, str1, str2);

  fp = fopen(fname, "w");
  
  if(fp == NULL)
  {
    printf("Cannot open file %s.\n", fname);
    exit(0);
  }

  for(i=0; i<nf_rebin; i++)
  {
    fprintf(fp, "%e %e\n", pow(10.0, freq_rebin[i]), pow(10.0, psd_rebin[i]));
  }
  fclose(fp);

  free(freq_rebin);
  free(psd_rebin);

  return;
}

//======================================================================

double psd_drw(double fk, double *arg)
{
  double A=exp(arg[0]), fknee=exp(arg[1]), cnoise = exp(arg[2]);

  return A/(1.0 + pow(fk/fknee, 2.0)) + cnoise;
}

double psd_drw_sqrt(double fk, double *arg)
{
  double A=exp(arg[0]), fknee=exp(arg[1]), cnoise = exp(arg[2]);

  return sqrt(A/(1.0 + pow(fk/fknee, 2.0))+ cnoise);
}

void psd_drw_array(double *fk, double *arg, double *psd, int n)
{
  double A=exp(arg[0]), fknee=exp(arg[1]), cnoise = exp(arg[2]);
  int i;

  for(i=0; i<n; i++)
  {
    psd[i] = A/(1.0 + pow(fk[i]/fknee, 2.0)) + cnoise;
  }
  return;
}

void psd_drw_sqrt_array(double *fk, double *arg, double *psd_sqrt, int n)
{
  double A=exp(arg[0]), fknee=exp(arg[1]), cnoise = exp(arg[2]);
  int i;
  for(i=0; i<n; i++)
  {
    psd_sqrt[i] = sqrt(A/(1.0 + pow(fk[i]/fknee, 2.0))+ cnoise);
  }
}

//======================================================================

double psd_power_law(double fk, double *arg)
{
  double A=exp(arg[0]), alpha=arg[1], cnoise=exp(arg[2]);

  if(fk > parset.freq_limit)
    return A * pow(fk, -alpha) + cnoise;
  else
    return A*pow(parset.freq_limit, -alpha) + cnoise;
    
}

double psd_power_law_sqrt(double fk, double *arg)
{
  double A=exp(arg[0]), alpha=arg[1], cnoise=exp(arg[2]);

  if(fk > parset.freq_limit)
    return sqrt(A * pow(fk, -alpha) + cnoise);
  else
    return sqrt(A*pow(parset.freq_limit, -alpha) + cnoise);
    
}


void psd_power_law_sqrt_array(double *fk, double *arg, double *psd_sqrt, int n)
{
  double A=exp(arg[0]), alpha=arg[1], cnoise=exp(arg[2]);
  int i;
  
  for(i=0; i<n; i++)
  {
    if(fk[i] > parset.freq_limit)
      psd_sqrt[i] = sqrt(A * pow(fk[i], -alpha) + cnoise);
    else
      psd_sqrt[i] = sqrt(A*pow(parset.freq_limit, -alpha) + cnoise);
  }
  return;
}

void psd_power_law_array(double *fk, double *arg, double *psd, int n)
{
  int i;
  psd_power_law_sqrt_array(fk, arg, psd, n);
  
  for(i=0; i<n; i++)
  {
    psd[i] = psd[i]*psd[i];
  }
  return;
}

//======================================================================

double psd_bending_power_law(double fk, double *arg)
{
  double A=exp(arg[0]), alpha_hi=arg[1], alpha_lo=(arg[1] - arg[2]);
  double fc=exp(arg[3]), cnoise=exp(arg[4]);

  if(fk > fc)
    return A * pow(fk/fc, -alpha_hi) + cnoise;
  else if(fk > parset.freq_limit)
    return A * pow(fk/fc, -alpha_lo) + cnoise;
  else
    return A * pow(parset.freq_limit/fc, -alpha_lo) + cnoise;
}

double psd_bending_power_law_sqrt(double fk, double *arg)
{
  double A=exp(arg[0]), alpha_hi=arg[1], alpha_lo=(arg[1] - arg[2]);
  double fc=exp(arg[3]), cnoise=exp(arg[4]);

  if(fk > fc)
    return sqrt(A * pow(fk/fc, -alpha_hi) + cnoise);
  else if(fk > parset.freq_limit)
    return sqrt(A * pow(fk/fc, -alpha_lo) + cnoise);
  else
    return sqrt(A * pow(parset.freq_limit/fc, -alpha_lo) + cnoise);
}


void psd_bending_power_law_sqrt_array(double *fk, double *arg, double *psd_sqrt, int n)
{
  double A=exp(arg[0]), alpha_hi=arg[1], alpha_lo=(arg[1] - arg[2]);
  double fc=exp(arg[3]), cnoise=exp(arg[4]);
  int i;

  for(i=0; i<n; i++)
  {
    if(fk[i] > fc)
      psd_sqrt[i] = sqrt(A * pow(fk[i]/fc, -alpha_hi) + cnoise);
    else if(fk[i] > parset.freq_limit)
      psd_sqrt[i] = sqrt(A * pow(fk[i]/fc, -alpha_lo) + cnoise);
    else
      psd_sqrt[i] = sqrt(A * pow(parset.freq_limit/fc, -alpha_lo) + cnoise);
  }
  return;
}

void psd_bending_power_law_array(double *fk, double *arg, double *psd, int n)
{
  int i;
  psd_bending_power_law_sqrt_array(fk, arg, psd, n);
  for(i=0; i<n; i++)
  {
    psd[i] = psd[i] * psd[i];
  }
  return;
}


//======================================================================

double psd_gaussian(double fk, double *arg)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), sigma=exp(arg[2]);

  return Ap * 1.0/sqrt(2.0*PI)/sigma * exp(-0.5 * pow( (fk - center)/sigma, 2.0 ));
}

double psd_gaussian_sqrt(double fk, double *arg)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), sigma=exp(arg[2]);

  return sqrt(Ap * 1.0/sqrt(sqrt(2.0*PI))/sigma * exp(-0.5 * pow( (fk - center)/sigma, 2.0 )));
}

void psd_gaussian_array(double *fk, double *arg, double *psd, int n)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), sigma=exp(arg[2]);
  int i;
  
  for(i=0; i<n; i++)
  {
    psd[i] = Ap * 1.0/sqrt(2.0*PI)/sigma * exp(-0.5 * pow( (fk[i] - center)/sigma, 2.0 ));
  }
}

void psd_gaussian_sqrt_array(double *fk, double *arg, double *psd_sqrt, int n)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), sigma=exp(arg[2]);
  int i;
  for(i=0; i<n; i++)
  {
    psd_sqrt[i] = sqrt(Ap * 1.0/sqrt(sqrt(2.0*PI))/sigma * exp(-0.5 * pow( (fk[i] - center)/sigma, 2.0 )));
  }
}

//======================================================================

double psd_lorentz(double fk, double *arg)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), width=exp(arg[2]);

  return Ap/PI *width/(width*width + pow(fk-center, 2.0));
}

double psd_lorentz_sqrt(double fk, double *arg)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), width=exp(arg[2]);

  return sqrt( Ap/PI *width/(width*width + pow(fk-center, 2.0)) );
}

void psd_lorentz_array(double *fk, double *arg, double *psd, int n)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), width=exp(arg[2]);
  int i;
  for(i=0; i<n; i++)
  {
    psd[i] = Ap/PI *width/(width*width + pow(fk[i]-center, 2.0));
  }
}

void psd_lorentz_sqrt_array(double *fk, double *arg, double *psd_sqrt, int n)
{
  double Ap=exp(arg[0]), center=exp(arg[1]), width=exp(arg[2]);
  int i;

  for(i=0; i<n; i++)
  {
    psd_sqrt[i] = sqrt( Ap/PI *width/(width*width + pow(fk[i]-center, 2.0)) );
  }
}