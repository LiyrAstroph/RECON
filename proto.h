
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
double psd_power_law(double fk, double *arg);
int recon_postproc();
void sim();
void print_help();
void read_parset();

/* psd */
int psddata_cal();
int psd_fft(double *t, double *f, int n, double *freq, double *psd, int *nf);
int psd_fft_rebin(double *freq, double *psd, int nf, double *freq_rebin, double *psd_rebin, int *nf_rebin);
int resample(double *t, double *f, int n, double *ts, double *fs);

/* system */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);

/* functions */
void from_prior_recon(void *model);
void print_particle_recon(FILE *fp, const void *model);
double log_likelihoods_cal_recon(const void *model);
double log_likelihoods_cal_initial_recon(const void *model);
double log_likelihoods_cal_restart_recon(const void *model);
double perturb_recon(void *model);
void copy_model_recon(void *dest, const void *src);
void* create_model_recon();
int get_num_params_recon();
void restart_clouds_recon(int iflag);

void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*log_likelihoods_cal_initial)(const void *model);
double (*log_likelihoods_cal_restart)(const void *model);
double (*perturb)(void *model);
void (*copy_model)(void *dest, const void *src);
void* (*create_model)();
int (*get_num_params)();
void (*restart_clouds)(int iflag);

double (*psdfunc)(double, double *);

#endif