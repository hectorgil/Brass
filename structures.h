

typedef struct f_params {

char *type_fog;
char *type_ptmodel;
char *type_rsdmodel;
char *type_bispectrummodel;

int mode;
//FOR RSD
int N;
double sigmaP;
double avir;
double sigmaB;
double b1;
double b2;
double bs2;
double b3nl;
double A;//same for B and P
double Pnoise;
double f;
double mBGV;
double m2BGV;
double sigma8;
double sigma8_value;
double a_parallel;
double a_perpendicular;
double **theory;
double kinput;
int noise_option;
//for Bispectrum
double knl;
double *k_n,*n;
double k1input,k2input,k3input;

//FOR BAO
int Nlin,NOlin;
double *parameters1;
int wiggle;
double Sigma_smooth;
double *k_Plin,*Plin,*k_Olin,*Olin;

int Npoly,modeP0,modeP2,modeP4;

char *spacing;

}f_params;


