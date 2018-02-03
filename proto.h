/*
 * RECON Copyright (C) 2018 Yan-Rong Li
 * A package for measuring spectral power and reconstructing time series in AGN.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 */

#ifndef _PROTO_H
#define _PROTO_H

int recon_init();
int recon_end();
int recon();
int genlc(const void *model);
int get_line_number(char *fname);
void get_num_particles(char *fname);
void get_posterior_sample_file(char *fname, char *samplefile);
int read_data(char *fname, int n, double *t, double *f, double *e);
void set_par_range();
double prob_con(const void *model);
double psd_drw(double fk, double *arg);
double psd_drw_sqrt(double fk, double *arg);
double psd_power_law(double fk, double *arg);
double psd_power_law_sqrt(double fk, double *arg);
double psd_period(double fk, double *arg);
double psd_period_sqrt(double fk, double *arg);
double psd_bending_power_law(double fk, double *arg);
double psd_bending_power_law_sqrt(double fk, double *arg);
int recon_postproc();
void sim();
void print_help();
void read_parset();
void time_cad_cal();
int recon_cmp(const void *a, const void *b);

/* psd */
int psddata_cal();
int psd_fft(double *t, double *f, int n, double *freq, double *psd, int *nf);
int psd_fft_rebin(double *freq, double *psd, int nf, double *freq_rebin, double *psd_rebin, int *nf_rebin);
int resample(double *t, double *f, int n, double *ts, double *fs);
void psd_fit_check(double *freq, double *psd, int nf);

/* system */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);
int remove_restart_file();

/* functions */
void from_prior_recon(void *model);
void print_particle_recon(FILE *fp, const void *model);
double log_likelihoods_cal_recon(const void *model);
double log_likelihoods_cal_initial_recon(const void *model);
double log_likelihoods_cal_restart_recon(const void *model);
double log_likelihoods_cal_recon_exam(const void *model);
double perturb_recon(void *model);
int get_num_params_recon();
void restart_clouds_recon(int iflag);

void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*log_likelihoods_cal_initial)(const void *model);
double (*log_likelihoods_cal_restart)(const void *model);
double (*perturb)(void *model);
int (*get_num_params)();
void (*restart_clouds)(int iflag);

double (*psdfunc)(double, double *);
double (*psdfunc_sqrt)(double, double *);

#endif