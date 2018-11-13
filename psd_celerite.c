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

#include "allvars.h"
#include "proto.h"

/*
 *
 * PSD for damped simple harmonic oscillator
 * see Foreman-Mackey et al. 2017, ApJ, 154, 220
 *  
 */
double psd_harmonic_sqrt(double fk, double *arg)
{
  double S1, S2, omega1_sqr, omega2_sqr, Q_sqr;
  double psd, fk_sqr, noise;
  int i;

  fk_sqr = fk*fk;
  S1 = exp(arg[0]);
  omega1_sqr = exp(2.0 * arg[1]);

  psd = S1*omega1_sqr*omega1_sqr/ ( (fk_sqr - omega1_sqr)*(fk_sqr - omega1_sqr) + 2*omega1_sqr*fk_sqr );
  for(i=1; i<parset.harmonic_term_num; i++)
  {
    S2 = exp(arg[2 + (i-1)*3]);
    omega2_sqr = exp(2.0 * arg[2 + (i-1)*3 + 1]);
    Q_sqr = pow(exp(arg[2 + (i-1)*3 + 2]) - 1.0, 2.0)/4.0;
    psd += S2*omega2_sqr*omega2_sqr/ ( (fk_sqr - omega2_sqr)*(fk_sqr - omega2_sqr) +   omega2_sqr*fk_sqr/Q_sqr );
  }
  
  noise = exp(arg[num_params_psd_tot - 1]);

  return sqrt(sqrt(2.0/PI) * psd + noise);
}

double psd_harmonic(double fk, double *arg)
{
  double psd_sqrt;
  psd_sqrt = psd_harmonic_sqrt(fk, arg);

  return psd_sqrt*psd_sqrt;
}

/*
 *
 * PSD for damped simple harmonic oscillator
 * see Foreman-Mackey et al. 2017, ApJ, 154, 220
 *  
 */
void psd_harmonic_sqrt_array(double *fk, double *arg, double *psd_sqrt, int n)
{
  double S1, *S2, omega1_sqr, *omega2_sqr, *Q_sqr;
  double fk_sqr, noise;
  int i, j;

  S2 = workspace_psd;
  omega2_sqr = S2 + parset.harmonic_term_num-1;
  Q_sqr = omega2_sqr + parset.harmonic_term_num-1;

  S1 = exp(arg[0]);
  omega1_sqr = exp(2.0 * arg[1]);
  noise = exp(arg[num_params_psd_tot - 1]);

  for(i=0; i<parset.harmonic_term_num-1; i++)
  {
    S2[i] = exp(arg[2 + (i)*3]);
    omega2_sqr[i] = exp(2.0 * arg[2 + (i)*3 + 1]);
    Q_sqr[i] = pow(exp(arg[2 + (i)*3 + 2]) - 1.0, 2.0)/4.0;
  }

  for(j=0; j<n; j++)
  {
    fk_sqr = fk[j]*fk[j];

    psd_sqrt[j] = S1*omega1_sqr*omega1_sqr/ ( (fk_sqr - omega1_sqr)*(fk_sqr - omega1_sqr) + 2*omega1_sqr*fk_sqr );
    for(i=0; i<parset.harmonic_term_num-1; i++)
    {
      psd_sqrt[j] += S2[i]*omega2_sqr[i]*omega2_sqr[i]/ ( (fk_sqr - omega2_sqr[i])*(fk_sqr - omega2_sqr[i]) +   omega2_sqr[i]*fk_sqr/Q_sqr[i] );
    }

    psd_sqrt[j] = sqrt(psd_sqrt[j] * sqrt(2.0/PI) + noise);
  }

  return;
}

void psd_harmonic_array(double *fk, double *arg, double *psd, int n)
{
  int i;
  psd_harmonic_sqrt_array(fk, arg, psd, n);
  for(i=0; i<n; i++)
  {
    psd[i] = psd[i] * psd[i];
  }
  return;
}