#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>


double Zspt(double k1, double k2, double k3, double knl, double *k_n, double *n,int N,double s8,double beta_F,double beta_G, double beta_S, double beta_mu,double C1,double C2);

double Zeff(double k1, double k2, double k3, double knl, double *k_n, double *n,int N,double s8,double beta_F,double beta_G, double beta_S, double beta_mu,double C1,double C2);

void do_Btheo_iso(double B_theo0[], int NeffB0, double *parameters1,double *k_Plin,double *Plin,int Nlin, double *k_Olin, double *Olin, int NOlin,double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *spacing_mask, char *path_to_mask1, double k_1[], double k_2[], double k_3[], int Npoly, fftw_plan plan1, fftw_plan plan2, double kmin_data0 , double kmax_data0, int wiggle,char *spacing_theory, double Sigma_smooth,double knl,double *k_n,double *n, double s8,char *bispectrum_BQ);

void do_Btheo_RSD(char *ptmodel, char *type_fog, char *type_bispectrummodel,double B_theo0[], int NeffB0,double *NoiseB, double *parameters1,double **Theory, int Nlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *spacing_mask, char *path_to_mask1, double k_1[], double k_2[], double k_3[], fftw_plan plan1, fftw_plan plan2, double kmin_data0 , double kmax_data0, char *spacing_theory,double knl,double *n, char *bispectrum_BQ, char *mask_matrix, double **MatrixFS_mask, int noise_option,double redshift_in);

void integrals8(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

double Wth(double x);

double get_fid_knl(double *klin, double *Plin, int Nlin);

double get_fid_sigma8(double *klin, double *Plin, int Nlin);

void generate_knl_array(double *klin, double *Plin, int Nlin, double *knl, double *knl_s8, int Nknl, double s8ref);

void generate_n_array(double *klin, double *Plin, int Nlin, double *n_func, double *k_n_func);

void smooth_n_array(double *k_in, int N_in, double *n_new, double *n_func, double *k_n_func, int Nfunc);

double M1(double beta, double k1, double k2, double k3);

double M2(double beta, double k1, double k2, double k3);

double M3(double beta, double k1, double k2, double k3);

double F2_eff(double k1, double k2, double k3, double knl, double **Theory, double *k_n, double *n,int N,double s8);

double G2_eff(double k1, double k2, double k3, double knl, double **Theory, double *k_n, double *n, int N, double s8);

double F2(double k1, double k2, double k3);

double G2(double k1, double k2, double k3);

double Q3(double n);

double aF(double n, double k, double knl,double s8);

double bF(double n, double k, double knl);

double cF(double n, double k, double knl);

double aG(double n, double k, double knl,double s8);

double bG(double n, double k, double knl);

double cG(double n, double k, double knl);

double S_ker(double k1, double k2, double k3);

double geo_factor(double z, double A, double ctheta_min, double ctheta_med, double ctheta_max);

double interpol_f(double z, double **f,int i);
