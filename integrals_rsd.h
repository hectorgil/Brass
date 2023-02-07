
void integralP(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void integralPaniso(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void integralB(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Preal(double *kin,double *Pin,int N, void *fdata);

double N2(double pl,double p13,double p15);

int determine_N1_aniso(double *Theory,double kinput,int Nlin, char *spacing);
