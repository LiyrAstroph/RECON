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
double psd_harmonic(double fk, double *arg)
{
  double S1, S2, omega1_sqr, omega2_sqr, Q_sqr;
  double psd, fk_sqr, noise;
  int i;

  fk_sqr = fk*fk;
  S1 = exp(arg[0]);
  omega1_sqr = exp(2.0 * arg[1]);

  psd = S1*omega1_sqr*omega1_sqr/ ( (fk_sqr - omega1_sqr)*(fk_sqr - omega1_sqr) + 2*omega1_sqr*fk_sqr );
  for(i=1; i<harmonic_term_num; i++)
  {
    S2 = exp(arg[2 + (i-1)*3]);
    omega2_sqr = exp(2.0 * arg[2 + (i-1)*3 + 1]);
    Q_sqr = pow(exp(arg[2 + (i-1)*3 + 2]) - 1.0, 2.0)/4.0;
    psd += S2*omega2_sqr*omega2_sqr/ ( (fk_sqr - omega2_sqr)*(fk_sqr - omega2_sqr) +   omega2_sqr*fk_sqr/Q_sqr );
  }
  
  noise = exp(arg[2 + (harmonic_term_num-1) * 3]);

  return sqrt(2.0/PI) * psd + noise;
}

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
  for(i=1; i<harmonic_term_num; i++)
  {
    S2 = exp(arg[2 + (i-1)*3]);
    omega2_sqr = exp(2.0 * arg[2 + (i-1)*3 + 1]);
    Q_sqr = pow(exp(arg[2 + (i-1)*3 + 2]) - 1.0, 2.0)/4.0;
    psd += S2*omega2_sqr*omega2_sqr/ ( (fk_sqr - omega2_sqr)*(fk_sqr - omega2_sqr) +   omega2_sqr*fk_sqr/Q_sqr );
  }
  
  noise = exp(arg[2 + (harmonic_term_num-1) * 3]);

  return sqrt(sqrt(2.0/PI) * psd + noise);
}