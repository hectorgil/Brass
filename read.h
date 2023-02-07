

int get_compression_file(char *path_to_compression, double **matrix_linear_compression, int input_elements_data_vector);


void get_Pk_theory(char *perturbation_theory_file,double **Theory,int N_Plin,int N_inputs);

void get_Pk_bao( char *path, double k[], double P[]);

void get_data_alphas(char *path, double *P0bao);

void get_data(char *path, double k0[],double kav0[],double P0[], double k2[],double kav2[], double P2[], double k4[], double kav4[], double P4[], double parameter_value[], double kminP0,double kmaxP0, double kminP2, double kmaxP2, double kminP4, double kmaxP4, char *type_BAORSD,int bao);

void get_data_bis(char *path, double k11[], double k22[], double k33[], double B0[], double Bnoise[], double kminB0, double kmaxB0, char *bispectrumBQ);

void get_mask(char *path,double posAV[], double pos[],double W0[], double W2[], double W4[], double W6[], double W8[], int Nmask, char *type_BAORSD, double params[],char *renormalize_window);

//void get_cov_from_file();

void write_cov_file(char *covfile, double Covariance[], int Ncov, int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao, int NeffP2rsd, int NeffP4bao, int NeffP4rsd, int NeffB0bao, int NeffB0rsd);

void get_cov_from_file(char *path_to_cov, char *path_to_cov_bis, double cov[], int Ncov,int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao,  int NeffP2rsd, int NeffP4bao, int NeffP4rsd , int NeffB0bao, int NeffB0rsd, double errP0bao[], double errP0rsd[],double errP2bao[], double errP2rsd[], double errP4bao[],  double errP4rsd[], double errB0bao[], double errB0rsd[], char *do_power_spectrum, char *do_bispectrum,double **matrix_compression,int input, int output, char *covariance_correction);

int get_cov_from_mocks(char *path_to_mocks_bao, char *path_to_mocks_rsd, char *path_to_mocks_bis_bao, char *path_to_mocks_bis_rsd,  char *covfile, double cov[], int Ncov,int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao,  int NeffP2rsd, int NeffP4bao, int NeffP4rsd , int NeffB0bao, int NeffB0rsd, double errP0bao[], double errP0rsd[],double errP2bao[], double errP2rsd[], double errP4bao[],  double errP4rsd[] , double errB0bao[], double errB0rsd[], double kminP0bao,  double kminP0rsd,  double kmaxP0bao,  double kmaxP0rsd,  double kminP2bao, double kminP2rsd ,double kmaxP2bao, double kmaxP2rsd,double kminP4bao, double kminP4rsd,double kmaxP4bao, double kmaxP4rsd,double kminB0bao, double kminB0rsd,double kmaxB0bao, double kmaxB0rsd, char *type_BAORSD, char *fit_BAO, char *fit_RSD, char *do_power_spectrum, char *do_bispectrum,char *bispectrumBQ, double **matrix_compression,int input, int output, char *covariance_correction, char *covarianceFSa_option, double *paramsBAO);

//void get_cov_from_mocks(char *path_to_mocks_bao, char *path_to_mocks_rsd, char *path_to_mocks_bis_bao, char *path_to_mocks_bis_rsd, double cov[], int Ncov,int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao,  int NeffP2rsd, int NeffP4bao, int NeffP4rsd , int NeffB0bao, int NeffB0rsd, double errP0bao[], double errP0rsd[],double errP2bao[], double errP2rsd[], double errP4bao[],  double errP4rsd[] , double errB0bao[], double errB0rsd[], double kminP0bao,  double kminP0rsd,  double kmaxP0bao,  double kmaxP0rsd,  double kminP2bao, double kminP2rsd ,double kmaxP2bao, double kmaxP2rsd,double kminP4bao, double kminP4rsd,double kmaxP4bao, double kmaxP4rsd,double kminB0bao, double kminB0rsd,double kmaxB0bao, double kmaxB0rsd, char *type_BAORSD, char *fit_BAO, char *fit_RSD, char *do_power_spectrum, char *do_bispectrum);


int countPk(int mode, char *path ,double kmin,double kmax);
