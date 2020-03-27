
void integralP(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void integralPaniso(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

//double P_interpol2_aniso_2ndorder(double *Theory, int N1,double w0, double w1, double w2);

//double P_interpol2_2ndorder(double **Theory,int index, int N1,double w0, double w1, double w2);

//double determine_w0_2ndorder_aniso(double *k,double kinput,int N1);

//double determine_w1_2ndorder_aniso(double *k,double kinput,int N1);

//double determine_w2_2ndorder_aniso(double *k,double kinput,int N1);

//double determine_w0_2ndorder(double *k,double kinput,int N1);

//double determine_w1_2ndorder(double *k,double kinput,int N1);

//double determine_w2_2ndorder(double *k,double kinput,int N1);

//int determine_N1(double **Theory,double kinput,int Nlin, char *spacing);

//double determine_w1(double **Theory,double kinput,int N1);

//double determine_w1_aniso(double *k,double kinput,int N1);

//double P_interpol2(double **Theory,int index, int N1,double w1);

//double P_interpol2_aniso(double *Theory, int N1,double w1);

double N2(double pl,double p13,double p15);

int determine_N1_aniso(double *Theory,double kinput,int Nlin, char *spacing);
