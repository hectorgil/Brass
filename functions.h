
//int smallest_common_multiple( int a, int b ); 

void determine_spacing_theo2(char *spacing,double **k0,int N);

void determine_spacing_theo(char *spacing,double *k0,int N);

void determine_spacing(char *spacing,double *k0, double *k2, double *k4,int NeffP0,int NeffP2,int NeffP4);

double get_error1(double *alpha, double *chi2, double chi2min, int imin, int N);

double get_error2(double *alpha, double *chi2, double chi2min, int imin, int N);

void string_copy(char *from, char *to);

double P_interpol(double k0, double *k, double *P, int N);

double P_interpol_doublearray(double k0, int index, double **P, int N);

double P_interpolLOG(double k0, double *k, double *P, int N);

void skipheader(FILE *myfile);

void skiplines(FILE *myfile, long int nlines);

void concatenate_files(int n_files, char *filename_wo_extension);

int countlines(char *filename);

long int countlinesLI(char *filename);

void freeTokens(double** tokens, int N);

void freeTokensInt(int** tokens, int N);

void freeTokens2(double ***tokens, int N1, int N2);

double P_interpol_fast(double k0, double *P, int N, char *spacing,int interpolation_order, int Ninterpol,double w0,double w1,double w2 );

double P_interpol_fast_doublearray(double k0, double **P, int index, int N, char *spacing,int interpolation_order, int Ninterpol,double w0,double w1,double w2 );

long int getMin(long int arr[], int n);

long int getMax(long int arr[], int n);


//

 int determine_N_doublearray(double **Theory,double kinput,int Nlin, char *spacing);

 int determine_N_singlearray(double *Theory,double kinput,int Nlin, char *spacing);


///
double determine_w1_doublearray(double **Theory,double kinput,int N1,char *spacing);

double determine_w1_singlearray(double *k,double kinput,int N1, char *spacing);

double P_interpol_w1_singlearray(double *Theory, int N1,double w1, char *spacing);

double P_interpol_w1_doublearray(double **Theory,int index, int N1,double w1,char *spacing);

//


double determine_w0_2ndorder_singlearray(double *k,double kinput,int N1, char *spacing);

double determine_w1_2ndorder_singlearray(double *k,double kinput,int N1,char *spacing);

double determine_w2_2ndorder_singlearray(double *k,double kinput,int N1,char *spacing);

double determine_w0_2ndorder_doublearray(double **k,double kinput,int N1,char *spacing);

double determine_w1_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing);

double determine_w2_2ndorder_doublearray(double **k,double kinput,int N1, char *spacing);

double P_interpol_w012_singlearray(double *Theory, int N1,double w0, double w1, double w2, char *spacing);

double P_interpol_w012_doublearray(double **Theory,int index, int N1,double w0, double w1, double w2,char *spacing);

