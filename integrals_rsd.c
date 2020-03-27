#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>
#include "functions.h"
#include "structures.h"
#include "rsd.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "fftlog.h"
#include "mcmc_bao.h"
#include "cubature.h"
#define Pi (4.*atan(1.))

double N2(double pl,double p13,double p15)
{
        double f;
    if (p15>=0)
    {
                f=cosh(sqrt(p15/pl))+0.5*p13/pl*sqrt(pl/p15)*sinh(sqrt(p15/pl));
    }
        else
        {
        f=cos(sqrt(-p15/pl))+0.5*p13/pl*sqrt(-pl/p15)*sin(sqrt(-p15/pl));
        }

        f=f*f;
        return f;
}





void integralPaniso(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

int i,j;
double func,olin;
int mode;
int Nlin,NOlin,Nlin1,NOlin1;
double w1lin,w1olin;
double w2lin,w2olin;
double w0lin,w0olin;
double *parameters1;
int wiggle;
double Sigma_smooth;
double *k_Plin,*Plin,*k_Olin,*Olin;
int Npoly;
int modeP0,modeP2,modeP4;
double alpha_perpendicular,alpha_parallel,kinput,mu,Fap,mup,k1p;
double sigma_para,sigma_perp,beta,R,Legendre;
int interpolation_order;
int shiftN;
double interP;
char *spacing;
interpolation_order=1;//1 or 2

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}


spacing=params_function.spacing;

    mode=params_function.mode;
    Nlin=params_function.Nlin;
    NOlin=params_function.NOlin;
    parameters1=params_function.parameters1;
    wiggle=params_function.wiggle;
    Sigma_smooth=params_function.Sigma_smooth;

    k_Plin=params_function.k_Plin;
    k_Olin=params_function.k_Olin;

    Plin=params_function.Plin;
    Olin=params_function.Olin;

    modeP0=params_function.modeP0;
    modeP2=params_function.modeP2;
    modeP4=params_function.modeP4;

    Npoly=params_function.Npoly;

    alpha_parallel=parameters1[0];
    alpha_perpendicular=parameters1[1];
    sigma_para=parameters1[2];
    sigma_perp=parameters1[3];
    beta=parameters1[4];

//printf("%lf %lf %lf %lf %lf\n",alpha_perpendicular,alpha_parallel, sigma_para,sigma_perp, beta);
//printf("%lf %lf %lf %lf %lf %lf %lf\n",parameters1[5],parameters1[6],parameters1[7],parameters1[8],parameters1[9],parameters1[10],parameters1[11]);
//printf("%d %d %d %d\n",mode,modeP0,modeP2,modeP4);
//get sigmas
//exit(0);
    kinput=params_function.kinput;
    mu=x[0];
    Fap=alpha_parallel*1./alpha_perpendicular*1.;
    mup=mu/Fap*pow(1.+mu*mu*(1./(Fap*Fap)-1),-0.5);
    k1p=kinput/alpha_perpendicular*pow( 1+mu*mu*(1./(Fap*Fap)-1) ,0.5);
//printf("%lf %lf\n",kinput,k1p);
//printf("N=%d %d\n",Nlin,NOlin);
    Nlin1=determine_N_singlearray(k_Plin,k1p,Nlin,spacing);
    NOlin1=determine_N_singlearray(k_Olin,k1p,NOlin,spacing);

//printf("N=%d %d\n",Nlin1,NOlin1);

if(Nlin1>=Nlin-shiftN || Nlin1<0 || NOlin1>=NOlin-shiftN || NOlin1<0 || k1p<=0 || kinput<=0)
{
fval[0]=0;//integral returns 0 out of range of Theory
//printf("Error: %d,%d (%d),%d\n",Nlin1,NOlin1,Nlin,NOlin);
//exit(0);
}
else{

if(interpolation_order==1){
w1lin=determine_w1_singlearray(k_Plin,k1p,Nlin1,spacing);
w1olin=determine_w1_singlearray(k_Olin,k1p,NOlin1,spacing);
}
if(interpolation_order==2){
w0lin=determine_w0_2ndorder_singlearray(k_Plin,k1p,Nlin1,spacing);
w1lin=determine_w1_2ndorder_singlearray(k_Plin,k1p,Nlin1,spacing);
w2lin=determine_w2_2ndorder_singlearray(k_Plin,k1p,Nlin1,spacing);

w0olin=determine_w0_2ndorder_singlearray(k_Olin,k1p,NOlin1,spacing);
w1olin=determine_w1_2ndorder_singlearray(k_Olin,k1p,NOlin1,spacing);
w2olin=determine_w2_2ndorder_singlearray(k_Olin,k1p,NOlin1,spacing);
}


if(mode==0){Legendre=1.0;}
if(mode==2){Legendre=3./2.*mu*mu-0.5;}
if(mode==4){Legendre=35./8.*mu*mu*mu*mu-30./8.*mu*mu+3./8.;}

//Legendre=1;
if(Sigma_smooth==0){//pre-recon case
R=1;
}
else{//post-recon
R=1-exp(-k1p*k1p*Sigma_smooth*Sigma_smooth*0.5);
}


interP=P_interpol_fast(k1p,Plin,Nlin,spacing,interpolation_order,Nlin1,w0lin,w1lin,w2lin);
olin=P_interpol_fast(k1p,Olin,NOlin,spacing,interpolation_order,NOlin1,w0olin,w1olin,w2olin);



              func=interP*parameters1[5]*pow(1+parameters1[4]*mup*mup*R,2)*Legendre; 

if(wiggle==0){olin=1;}

//fval[0]=0.5*(2.*mode+1.)/(alpha_parallel*alpha_perpendicular*alpha_perpendicular)*func*(1+(olin-1)*exp(-0.5*(k1p*k1p*mup*mup*sigma_para*sigma_para+k1p*k1p*(1.-mup*mup)*sigma_perp*sigma_perp)));

fval[0]=0.5*(2.*mode+1.)*func*(1+(olin-1)*exp(-0.5*(k1p*k1p*mup*mup*sigma_para*sigma_para+k1p*k1p*(1.-mup*mup)*sigma_perp*sigma_perp)));
//printf("%lf %lf, %d %lf, %lf,%lf, %lf %lf\n",k1p,interP,Nlin1,w1lin, k_Plin[Nlin1],k_Plin[Nlin1+1],Plin[Nlin1],Plin[Nlin1+1]);
//exit(0);
}


}

void integralP(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;
    
    int N, nomask;
    double b1,b2,bs2,b3nl,Anoise,Pnoise;
    double alpha_perpendicular,alpha_parallel,F,sigma8_scaling,sigma_ps;
    double **Theory;
    double a11,a12,a22,a23,a33;
    double b1_11,b1_12,b1_21,b1_22;
    double b2_11,b2_12,b2_21,b2_22;
    double b3_12,b3_21,b3_22,b4_22;
    double Legendre, Pdd, Ptt, Pdt, Pb2_d, Pb2_t, Pbs2_d, Pbs2_t, sigma3, Pb22, Pb2s2, Pbs22;
    double fog;
    char *fog_model_ps, *ptmodel, *rsdmodel_ps;
    int mode;
    double kinput, mu, Fap, mup;
    double k1p;
    int N1;
    double w1,w0,w2;
    int interpolation_order,shiftN;
    char *spacing;

    interpolation_order=1;//1 or 2
    if(interpolation_order==1){shiftN=1;}
    if(interpolation_order==2){shiftN=2;}

    spacing=params_function.spacing;
    N=params_function.N;
    b1=params_function.b1;
    b2=params_function.b2;
    bs2=params_function.bs2;
    b3nl=params_function.b3nl;
    Anoise=params_function.A;
    Pnoise=params_function.Pnoise;
    alpha_perpendicular=params_function.a_perpendicular;
    alpha_parallel=params_function.a_parallel;
    F=params_function.f;
    sigma8_scaling=params_function.sigma8;
    sigma_ps=params_function.sigmaP;
    Theory=params_function.theory;

    fog_model_ps=params_function.type_fog;
    ptmodel=params_function.type_ptmodel;
    rsdmodel_ps=params_function.type_rsdmodel;

    mode=params_function.mode;
    kinput=params_function.kinput;
    mu=x[0];
    Fap=alpha_parallel*1./alpha_perpendicular*1.;
    mup=mu/Fap*pow(1.+mu*mu*(1./(Fap*Fap)-1),-0.5);
    k1p=kinput/alpha_perpendicular*pow( 1+mu*mu*(1./(Fap*Fap)-1) ,0.5);

    N1=determine_N_doublearray(Theory,k1p,N,spacing);

if(N1>=N-shiftN || N1<0 || k1p<=0 || kinput<=0)
{
fval[0]=0;//integral returns 0 out of range of Theory
}
else{
    //w1=determine_w1(Theory,k1p,N1);

if(interpolation_order==1){
    w1=determine_w1_doublearray(Theory,k1p,N1,spacing);
}
if(interpolation_order==2){

w0=determine_w0_2ndorder_doublearray(Theory,k1p,N1,spacing);
w1=determine_w1_2ndorder_doublearray(Theory,k1p,N1,spacing);
w2=determine_w2_2ndorder_doublearray(Theory,k1p,N1,spacing);


}

if( strcmp(rsdmodel_ps, "TNS10") ==0){

    a11=P_interpol_fast_doublearray(k1p,Theory,24,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*b1*F;
    a12=P_interpol_fast_doublearray(k1p,Theory,25,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F;
    a22=P_interpol_fast_doublearray(k1p,Theory,26,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F;
    a23=P_interpol_fast_doublearray(k1p,Theory,27,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F;
    a33=P_interpol_fast_doublearray(k1p,Theory,28,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F;

    b1_11=P_interpol_fast_doublearray(k1p,Theory,29,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*b1*F*F;
    b1_12=P_interpol_fast_doublearray(k1p,Theory,30,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b1_21=P_interpol_fast_doublearray(k1p,Theory,31,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b1_22=P_interpol_fast_doublearray(k1p,Theory,32,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;

    b2_11=P_interpol_fast_doublearray(k1p,Theory,33,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*b1*F*F;
    b2_12=P_interpol_fast_doublearray(k1p,Theory,34,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b2_21=P_interpol_fast_doublearray(k1p,Theory,35,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b2_22=P_interpol_fast_doublearray(k1p,Theory,36,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;

    b3_12=P_interpol_fast_doublearray(k1p,Theory,37,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b3_21=P_interpol_fast_doublearray(k1p,Theory,38,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b3_22=P_interpol_fast_doublearray(k1p,Theory,39,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;
    b4_22=P_interpol_fast_doublearray(k1p,Theory,40,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;

}


if( strcmp(ptmodel, "linear") == 0)
{

Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2);
Pdt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2);
Ptt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2);

}

if( strcmp(ptmodel, "1L-SPT") == 0)
{

Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2));
Pdt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2));
Ptt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2));

}


if( strcmp(ptmodel, "2L-SPT") == 0)
{

Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2))+sigma8_scaling*sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,5,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,12,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)/(4.*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)));

Pdt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))+sigma8_scaling*sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,6,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,14,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,15,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2)+0.25*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))/(4.*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)));

Ptt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))+sigma8_scaling*sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,7,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,13,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)/(4.*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)));


}

if( strcmp(ptmodel, "1L-RPT") == 0)
{

Pdd = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2))*exp(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)/P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2));
Pdt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2))*exp(sigma8_scaling*0.5*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))/P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2));
Ptt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2))*exp(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)/P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2));

}


if( strcmp(ptmodel, "2L-RPT") == 0)
{

Pdd =(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,5,N,spacing,interpolation_order,N1,w0,w1,w2))*N2(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2));
Pdt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,6,N,spacing,interpolation_order,N1,w0,w1,w2))*N2(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*0.5*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)),sigma8_scaling*sigma8_scaling*sigma8_scaling*0.5*(P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2)));
Ptt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,7,N,spacing,interpolation_order,N1,w0,w1,w2))*N2(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2));


}



     Pb2_d=P_interpol_fast_doublearray(k1p,Theory,16,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pb2_t=P_interpol_fast_doublearray(k1p,Theory,17,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pbs2_d=P_interpol_fast_doublearray(k1p,Theory,18,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pbs2_t=P_interpol_fast_doublearray(k1p,Theory,19,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     sigma3=P_interpol_fast_doublearray(k1p,Theory,23,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pb22=P_interpol_fast_doublearray(k1p,Theory,22,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pb2s2=P_interpol_fast_doublearray(k1p,Theory,20,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pbs22=P_interpol_fast_doublearray(k1p,Theory,21,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;


if( strcmp(fog_model_ps, "Lorentzian") ==0){fog=pow( 1+sigma_ps*sigma_ps/2.*k1p*k1p*mup*mup ,-2);}
if( strcmp(fog_model_ps, "Exponential") ==0){fog=exp( -sigma_ps*sigma_ps/2.*k1p*k1p*mup*mup);}

if(fabs(sigma_ps)<0.001){fog=1.;}

if(mode==0){Legendre=1.0;nomask=1;}
if(mode==2){Legendre=3./2.*mu*mu-0.5;nomask=0;}
if(mode==4){Legendre=35./8.*mu*mu*mu*mu-30./8.*mu*mu+3./8.;nomask=0;}


if( strcmp(rsdmodel_ps, "TNS10") ==0){

              fval[0]=(2*mode+1)/2.*((b1*b1*Pdd+2.*b2*b1*Pb2_d+2.*b1*bs2*Pbs2_d+2.*sigma3*b1*b3nl+b2*b2*Pb22+2.*b2*bs2*Pb2s2+bs2*bs2*Pbs22)+2.*F*mup*mup*( b1*Pdt+b2*Pb2_t+bs2*Pbs2_t+sigma3*b3nl )+pow(mup,4)*F*F*Ptt+a11*mup*mup+a12*mup*mup+a22*mup*mup*mup*mup+a23*mup*mup*mup*mup+a33*mup*mup*mup*mup*mup*mup+b1_11*mup*mup+b1_12*mup*mup+b1_21*mup*mup+b1_22*mup*mup+b2_11*mup*mup*mup*mup+b2_12*mup*mup*mup*mup+b2_21*mup*mup*mup*mup+b2_22*mup*mup*mup*mup+b3_12*mup*mup*mup*mup*mup*mup+b3_21*mup*mup*mup*mup*mup*mup+b3_22*mup*mup*mup*mup*mup*mup+b4_22*mup*mup*mup*mup*mup*mup*mup*mup)*fog*Legendre/(alpha_parallel*alpha_perpendicular*alpha_perpendicular)+nomask*Anoise*Pnoise*(2*mode+1)/(2*alpha_parallel*alpha_perpendicular*alpha_perpendicular)*Legendre;

}


if( strcmp(rsdmodel_ps, "Scoccimarro04") ==0){
              fval[0]=(2*mode+1)/2.*((b1*b1*Pdd+2.*b2*b1*Pb2_d+2.*b1*bs2*Pbs2_d+2.*sigma3*b1*b3nl+b2*b2*Pb22+2.*b2*bs2*Pb2s2+bs2*bs2*Pbs22)+2.*F*mup*mup*( b1*Pdt+b2*Pb2_t+bs2*Pbs2_t+sigma3*b3nl )+pow(mup,4)*F*F*Ptt)*fog*Legendre/(alpha_parallel*alpha_perpendicular*alpha_perpendicular)+nomask*Anoise*Pnoise*(2*mode+1)/(2*alpha_parallel*alpha_perpendicular*alpha_perpendicular)*Legendre;

}

if( strcmp(rsdmodel_ps, "Kaiser87") ==0){

 fval[0]=(2*mode+1)/2.*((b1*b1*Pdd+2.*b2*b1*Pb2_d+2.*b1*bs2*Pbs2_d+2.*sigma3*b1*b3nl+b2*b2*Pb22+2.*b2*bs2*Pb2s2+bs2*bs2*Pbs22)*(1+F/b1*mup*mup)*(1+F/b1*mup*mup))*fog*Legendre/(alpha_parallel*alpha_perpendicular*alpha_perpendicular)+nomask*Anoise*Pnoise*(2*mode+1)/(2*alpha_parallel*alpha_perpendicular*alpha_perpendicular)*Legendre;

}


}


}

