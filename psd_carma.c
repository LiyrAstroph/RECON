#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "allvars.h"
#include "proto.h"

/*
 * PSD for a CARMA proces
 */

void get_ar_roots(double *theta, complex *roots)
{
  int i;
  double a, b;

  for(i=0; i<parset.carma_p/2; i++)
  {
    a = exp(theta[1+i*2]);
    b = exp(theta[1+i*2+1]);
    roots[i*2] = -2.0*PI* (a +  b * I);
    roots[i*2+1] = -2.0*PI* (a - b * I);
  }

  if(parset.carma_p % 2 == 1)
  {
    roots[parset.carma_p-1] = -2.0*PI * exp(theta[1 + parset.carma_p -1]);
  }
  return;
}

void get_ma_coefs(double *theta, double * ma_coefs)
{
  int i;
  ma_coefs[0] = 1.0;
  for(i=1; i<=parset.carma_q; i++)
  {
    ma_coefs[i] = exp(theta[1 + parset.carma_p + i-1]);
  }
  for(i=parset.carma_q+1; i<=parset.carma_p; i++)
  {
    ma_coefs[i] = 0.0;
  }
  return;
}

void get_poly_coefs(complex *roots, int n, double *coefs)
{
  int i, j;
  complex *pcoefs = workspace_complex + parset.carma_p;

  pcoefs[0] = 1.0;
  for(i=0; i<=n; i++)
  {
    pcoefs[i+1] = 0.0;
    for(j=i+1; j>=1; j--)
      pcoefs[j] = pcoefs[j] - roots[i] * pcoefs[j-1];
  }

  for(i=0; i <=n; i++)
  {
    coefs[i] = creal(pcoefs[n-i]);
  }
  return;
}

/*
 * orders of arg: sigma, Lorentizan width, Lorentizn centriod, moving-average coefficients
 *
 */

double psd_carma_sqrt(double fk, double *arg)
{
  double sigma, psd_sqrt, noise;
  double *ar_coefs, *ma_coefs;
  complex *ar_roots;
  int i; 
  complex ar_poly, ma_poly, tmp;

  ar_roots = workspace_complex;
  ar_coefs = workspace_psd;
  ma_coefs = ar_coefs + parset.carma_p+1;
  
  sigma = exp(arg[0]);
  noise = exp(arg[num_params_psd - 1]);

  get_ar_roots(arg, ar_roots);
  get_poly_coefs(ar_roots, parset.carma_p, ar_coefs);
  get_ma_coefs(arg, ma_coefs);

  ar_poly = 0.0;
  ma_poly = 0.0;
  for(i=0; i<=parset.carma_p; i++)
  {
    tmp = cpow(2.0*PI*fk*I, i);
    ar_poly += ar_coefs[i] * tmp;
    ma_poly += ma_coefs[i] * tmp;
  }

  psd_sqrt = sigma * cabs(ma_poly)/cabs(ar_poly);
  psd_sqrt = sqrt(psd_sqrt*psd_sqrt + noise);
  return psd_sqrt;
}

double psd_carma(double fk, double *arg)
{
  double psd_sqrt;
  psd_sqrt = psd_carma(fk, arg);

  return psd_sqrt*psd_sqrt;
}

void psd_carma_sqrt_array(double *fk, double *arg, double *psd_sqrt, int n)
{
  double sigma, noise, tmp_sqrt;
  double *ar_coefs, *ma_coefs;
  complex *ar_roots;
  int i, j; 
  complex ar_poly, ma_poly, tmp;

  ar_roots = workspace_complex;
  ar_coefs = workspace_psd;
  ma_coefs = ar_coefs + parset.carma_p+1;
  
  sigma = exp(arg[0]);
  noise = exp(arg[num_params_psd - 1]);

  get_ar_roots(arg, ar_roots);
  get_poly_coefs(ar_roots, parset.carma_p, ar_coefs);
  get_ma_coefs(arg, ma_coefs);

  for(j=0; j<n; j++)
  {
    ar_poly = ar_coefs[0];
    ma_poly = ma_coefs[0];
    for(i=1; i<=parset.carma_p; i++)
    {
      tmp = freq_array_pow[i-1][j];
      ar_poly += ar_coefs[i] * tmp;
      ma_poly += ma_coefs[i] * tmp;
    }
    tmp_sqrt = sigma* cabs(ma_poly)/cabs(ar_poly);
    psd_sqrt[j] = sqrt( tmp_sqrt*tmp_sqrt + noise );
  }

  return;
}

void psd_carma_array(double *fk, double *arg, double *psd, int n)
{
  int i;
  psd_carma_sqrt_array(fk, arg, psd, n);
  for(i=0; i<n; i++)
  {
    psd[i] = psd[i]*psd[i];
  }
  
  return;
}