

typedef struct f_params {

char *type_fog;
char *type_ptmodel;
char *type_rsdmodel;

int mode;
//FOR RSD
int N;
double sigmaP;
double b1;
double b2;
double bs2;
double b3nl;
double A;
double Pnoise;
double f;
double sigma8;
double a_parallel;
double a_perpendicular;
double **theory;
double kinput;


//FOR BAO
int Nlin,NOlin;
double *parameters1;
int wiggle;
double Sigma_smooth;
double *k_Plin,*Plin,*k_Olin,*Olin;

int Npoly,modeP0,modeP2,modeP4;

char *spacing;

}f_params;


