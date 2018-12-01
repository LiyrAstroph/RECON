/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

#ifndef _PROTO_H
#define _PROTO_H
#include <complex.h>

int recon_init();
int recon_end();
double recon();
double recon_run(int argc, char **argv);
void genlc(const void *model);
void genlc_array(const void *model);
void genlc_array_initial(const void *model);
int get_line_number(char *fname);
void get_num_particles(char *fname);
void get_max_num_levels(char *fname);
void get_posterior_sample_file(char *fname, char *samplefile);
int read_data(char *fname, int n, double *t, double *f, double *e);
void set_par_range();
void set_par_fix();
void set_psd_functions();
double prob_recon(const void *model);
double prob_initial_recon(const void *model);
int recon_postproc();
void sim();
void print_help();
void read_parset();
void read_sim_arg();
void time_cad_cal();
int recon_cmp(const void *a, const void *b);

/* psd */
int psddata_cal();
int psd_fft(double *t, double *f, int n, double *freq, double *psd, int *nf);
int psd_fft_rebin(double *freq, double *psd, int nf, double *freq_rebin, double *psd_rebin, int *nf_rebin);
int resample(double *t, double *f, int n, double *ts, double *fs);
void psd_fit_check(double *freq, double *psd, int nf);
double psd_drw(double fk, double *arg);
double psd_drw_sqrt(double fk, double *arg);
double psd_power_law(double fk, double *arg);
double psd_power_law_sqrt(double fk, double *arg);
double psd_gaussian(double fk, double *arg);
double psd_gaussian_sqrt(double fk, double *arg);
double psd_lorentz(double fk, double *arg);
double psd_lorentz_sqrt(double fk, double *arg);
double psd_bending_power_law(double fk, double *arg);
double psd_bending_power_law_sqrt(double fk, double *arg);

void psd_drw_array(double *fk, double *arg, double *psd, int n);
void psd_drw_sqrt_array(double *fk, double *arg, double *psd, int n);
void psd_power_law_array(double *fk, double *arg, double *psd, int n);
void psd_power_law_sqrt_array(double *fk, double *arg, double *psd, int n);
void psd_gaussian_array(double *fk, double *arg, double *psd, int n);
void psd_gaussian_sqrt_array(double *fk, double *arg, double *psd, int n);
void psd_lorentz_array(double *fk, double *arg, double *psd, int n);
void psd_lorentz_sqrt_array(double *fk, double *arg, double *psd, int n);
void psd_bending_power_law_array(double *fk, double *arg, double *psd, int n);
void psd_bending_power_law_sqrt_array(double *fk, double *arg, double *psd, int n);

double psd_harmonic(double fk, double *arg);
double psd_harmonic_sqrt(double fk, double *arg);
void psd_harmonic_array(double *fk, double *arg, double *psd, int n);
void psd_harmonic_sqrt_array(double *fk, double *arg, double *psd, int n);

void get_ar_roots(double *theta, complex *roots);
void get_ma_coefs(double *theta, double * ma_coefs);
void get_poly_coefs(complex *roots, int n, double *coefs);
double psd_carma(double fk, double *arg);
double psd_carma_sqrt(double fk, double *arg);
void psd_carma_array(double *fk, double *arg, double *psd, int n);
void psd_carma_sqrt_array(double *fk, double *arg, double *psd, int n);


/* system */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);
int remove_restart_file();

/* functions */
void from_prior_recon(void *model);
void print_particle_recon(FILE *fp, const void *model);
void read_particle_recon(FILE *fp, void *model);
void print_particle_recon_saveoutput(FILE *fp, const void *model);
void read_particle_recon_saveoutput(FILE *fp, void *model);
double log_likelihoods_cal_recon(const void *model);
double log_likelihoods_cal_initial_recon(const void *model);
double log_likelihoods_cal_recon_exam(const void *model);
double perturb_recon(void *model);
double perturb_recon_limits(void *model);
void restart_action_recon(int iflag);
void accept_action_recon();
void kill_action_recon(int i, int i_copy);

double (*psdfunc)(double, double *);
double (*psdfunc_sqrt)(double, double *);

double (*psdfunc_period)(double, double *);
double (*psdfunc_period_sqrt)(double, double *);

void (*psdfunc_array)(double *, double *, double *, int );
void (*psdfunc_sqrt_array)(double *, double *, double *, int );

void (*psdfunc_period_array)(double *, double *, double *, int );
void (*psdfunc_period_sqrt_array)(double *, double *, double *, int );
#endif