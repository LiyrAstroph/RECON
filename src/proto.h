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

extern void init();
extern void end();
extern void set_var_range();

extern int recon_init();
extern int recon_end();
extern double recon();
extern double recon_run(int argc, char **argv);
extern void genlc(const void *model);
extern void genlc_array(const void *model);
extern void genlc_array_initial(const void *model);
extern void genlc_array_simple_period(const void *model);
extern void genlc_array_simple_delta(const void *model);
extern int get_line_number(char *fname);
extern void get_num_particles(char *fname);
extern void get_max_num_levels(char *fname);
extern int read_data(char *fname, int n, double *t, double *f, double *e);
extern void set_par_range();
extern void set_par_fix();
extern void set_psd_functions();
extern double prob_recon(const void *model);
extern double prob_initial_recon(const void *model);
extern int recon_postproc();
extern void sim();
extern void print_help();
extern void read_parset();
extern void read_sim_arg();
extern void time_cad_cal();
extern int recon_cmp(const void *a, const void *b);

/* psd fit*/
extern double psd_fit();
extern void psd_fit_init();
extern void psd_fit_end();
extern void from_prior_fit(void *model);
extern double log_likelihoods_cal_fit(const void *model);
extern double log_likelihoods_cal_initial_fit(const void *model);
extern double perturb_fit(void *model);
extern void set_par_fix_fit();
extern void set_par_range_fit();
extern double prob_fit(const void *model);
extern void psd_fit_postproc();
extern void restart_action_fit(int iflag);
extern void kill_action_fit(int i, int i_copy);
extern void accept_action_fit();

/* psd */
extern int psddata_cal();
extern int psd_fft(double *t, double *f, int n, double *freq, double *psd, int *nf);
extern int psd_fft_rebin(double *freq, double *psd, int nf, double *freq_rebin, double *psd_rebin, int *nf_rebin);
extern int resample(double *t, double *f, int n, double *ts, double *fs);
extern void psd_fit_check(double *freq, double *psd, int nf);
extern double psd_drw(double fk, double *arg);
extern double psd_drw_sqrt(double fk, double *arg);
extern double psd_power_law(double fk, double *arg);
extern double psd_power_law_sqrt(double fk, double *arg);
extern double psd_gaussian(double fk, double *arg);
extern double psd_gaussian_sqrt(double fk, double *arg);
extern double psd_lorentz(double fk, double *arg);
extern double psd_lorentz_sqrt(double fk, double *arg);
extern double psd_bending_power_law(double fk, double *arg);
extern double psd_bending_power_law_sqrt(double fk, double *arg);

extern void psd_drw_array(double *fk, double *arg, double *psd, int n);
extern void psd_drw_sqrt_array(double *fk, double *arg, double *psd, int n);
extern void psd_power_law_array(double *fk, double *arg, double *psd, int n);
extern void psd_power_law_sqrt_array(double *fk, double *arg, double *psd, int n);
extern void psd_gaussian_array(double *fk, double *arg, double *psd, int n);
extern void psd_gaussian_sqrt_array(double *fk, double *arg, double *psd, int n);
extern void psd_lorentz_array(double *fk, double *arg, double *psd, int n);
extern void psd_lorentz_sqrt_array(double *fk, double *arg, double *psd, int n);
extern void psd_bending_power_law_array(double *fk, double *arg, double *psd, int n);
extern void psd_bending_power_law_sqrt_array(double *fk, double *arg, double *psd, int n);

extern double psd_harmonic(double fk, double *arg);
extern double psd_harmonic_sqrt(double fk, double *arg);
extern void psd_harmonic_array(double *fk, double *arg, double *psd, int n);
extern void psd_harmonic_sqrt_array(double *fk, double *arg, double *psd, int n);

extern void get_ar_roots(double *theta, complex *roots);
extern void get_ma_coefs(double *theta, double * ma_coefs);
extern void get_poly_coefs(complex *roots, int n, double *coefs);
extern double psd_carma(double fk, double *arg);
extern double psd_carma_sqrt(double fk, double *arg);
extern void psd_carma_array(double *fk, double *arg, double *psd, int n);
extern void psd_carma_sqrt_array(double *fk, double *arg, double *psd, int n);


/* system */
extern double second();
extern double timediff(double t0, double t1);
extern void get_hms(double dt, int *h, int *m, double *s);
extern int remove_restart_file();

/* functions */
extern void from_prior_recon(void *model);
extern void print_particle_recon(FILE *fp, const void *model);
extern void read_particle_recon(FILE *fp, void *model);
extern void print_particle_recon_saveoutput(FILE *fp, const void *model);
extern void read_particle_recon_saveoutput(FILE *fp, void *model);
extern double log_likelihoods_cal_recon(const void *model);
extern double log_likelihoods_cal_initial_recon(const void *model);
extern double log_likelihoods_cal_recon_exam(const void *model);
extern double perturb_recon(void *model);
extern double perturb_recon_limits(void *model);
extern void restart_action_recon(int iflag);
extern void accept_action_recon_period();
extern void accept_action_recon_other();
extern void kill_action_recon_period(int i, int i_copy);
extern void kill_action_recon_other(int i, int i_copy);

extern double (*psdfunc)(double, double *);
extern double (*psdfunc_sqrt)(double, double *);

extern double (*psdfunc_period)(double, double *);
extern double (*psdfunc_period_sqrt)(double, double *);

extern void (*psdfunc_array)(double *, double *, double *, int );
extern void (*psdfunc_sqrt_array)(double *, double *, double *, int );

extern void (*psdfunc_period_array)(double *, double *, double *, int );
extern void (*psdfunc_period_sqrt_array)(double *, double *, double *, int );

extern void (*func_genlc_array)(const void *);
#endif