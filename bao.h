
#include <fftw3.h>

void set_mask_params(double params[],double kmin,double kmax,int Nlin, double kmin_data, double kmax_data);


void apply_mask(char *type_of_analysis,int modeP0, int modeP2, int modeP4, double *k_theo,double P_theo0[], double P_theo2[], double P_theo4[],int Nlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *sampling_mask, fftw_plan plan1, fftw_plan plan2,double kmin, double kmax, double kmin_data0, double kmax_data0, double kmin_data2, double kmax_data2, double kmin_data4, double kmax_data4, char *sampling_theo, double Sigma_smooth,int isbao);


void  make_a_bao_plot(char *type_BAO_fit, char *type_of_analysis,char *fit_BAO,double *parameters2, double chi2_min, double *k0, double *P0, double *errP0,int NeffP0, double *k0SGC, double *P0SGC, double *errP0SGC, int NeffP0SGC, double *k2, double *P2, double *errP2,int NeffP2, double *k2SGC, double *P2SGC, double *errP2SGC, int NeffP2SGC, double *k4, double *P4, double *errP4, int NeffP4, double *k4SGC, double *P4SGC, double *errP4SGC, int NeffP4SGC, double *k_Plin, double *Plin, int Nlin, double *k_Olin, double *Olin, int NOlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC,double *W4SGC, double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC,  char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks,int Npoints,int Ndof, char *path_output, char *identifier,fftw_plan plan1, fftw_plan plan2,int Nalphas,int Nsigmas_tot, int Nsigmas_free,double Sigma_smooth,int factor_sampling_mask,char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory);


double gauss(double x,double Sigma_nl_mean,double Sigma_nl_stddev);


double chi2_bao(char *type_BAO_fit,char *type_of_analysis,char *fit_BAO,double *parameters2, double *k_Plin, double *Plin,int Nlin,double *k_Olin,double *Olin,int NOlin, double *pos, double *W0,double *W2,double *W4,double *W6,double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC,double *W2SGC,double *W4SGC,double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0, double *P0,int NeffP0,double *k2, double *P2,int NeffP2,double *k4,double *P4, int NeffP4, double *k0SGC, double *P0SGC, int NeffP0SGC,double *k2SGC, double *P2SGC,int NeffP2SGC, double *k4SGC, double *P4SGC, int NeffP4SGC, double *cov, double *covSGC,  char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks, fftw_plan plan1, fftw_plan plan2, char *do_power_spectrum, char *do_bispectrum,int Nalphas,int Nsigmas_tot, int Nsigmas_free, double Sigma_smooth,int factor_sampling_mask,char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory);

void do_bao_analytic(char *type_BAO_fit,char *type_of_analysis,char *fit_BAO,double *k_Plin, double *Plin, int Nlin, double *k_Olin, double *Olin, int NOlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC,char *path_to_mask2, char *spacing_maskSGC,  double *k0, double *P0, double *errP0,  int NeffP0, double *k2, double *P2, double *errP2,  int NeffP2,  double *k4, double *P4, double *errP4,  int NeffP4, double *k0SGC, double *P0SGC, double *errP0SGC, int NeffP0SGC, double *k2SGC, double *P2SGC, double *errP2SGC, int NeffP2SGC,double *k4SGC, double *P4SGC, double *errP4SGC, int NeffP4SGC, double *cov, double *covSGC, double alpha_min, double alpha_max, double alpha_step, char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks, char *path_output, char *identifier, char *do_plot, char *do_power_spectrum, char *do_bispectrum,char *spacing_dataNGC,char *spacing_dataSGC, char *spacing_theory);


double hi0_mask(int index,int mode,int modeP0,int modeP2,int modeP4, double ki, double alpha_0,double alpha_00, double Sigma_nl,int mode_parameter, double *pos, double *W0,double *W2,double *W4,double *W6,double *W8, int Nmask, char *path_to_mask, double *klin, double *Plin, int Nlin, double *kolin, double *Polin, int Nolin,fftw_plan plan1, fftw_plan plan2,double kmin,double kmax,double kmin_data,double kmax_data);

double hi0(int mode,int modeP0,int modeP2,int modeP4, double ki, double alpha_0,double alpha_00, double Sigma_nl,int mode_parameter, double *klin, double *Plin, int Nlin, double *kolin, double *Polin, int Nolin);

void parameter_bestfit(int mode, int modeP0,int modeP2,int modeP4,double alpha_0,double alpha_00,double Sigma_nl, double *Knl, double *Pk,double *Cov,int Ncov, double A[],int Npoly, double *pos, double *W0,double *W2, double *W4, double *W6,double *W8, int Nmask, char *path_to_mask, double *klin, double *Plin, int Nlin, double *kolin, double *Polin, int Nolin,fftw_plan plan1, fftw_plan plan2);

void do_log_file2(int nthreads, char *path_output,char *identifier,char *likelihood_file,char *fit_BAO,int Nchunks,int Npolynomial,double time_run,int processors,double amin,double amax,double astep, int Nalphas,int Nsigmas);

