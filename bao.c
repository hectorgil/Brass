#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>
#include "functions.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "fftlog.h"
#include "mcmc_bao.h"
#define Pi (4.*atan(1.))

void set_mask_params(double params[],double kmin,double kmax,int Nlin,double kmin_data, double kmax_data)
{

//these parameters have been tested with functions between 0.005 and 0.5 and around 5000 points

if(kmin_data<kmin || kmax_data>kmax){printf("Error, the data k-range (%lf,%lf) exceeds the Plin k-range (%lf,%lf). Exiting now...\n",kmin_data,kmax_data,kmin,kmax);exit(0);}

double kpad_up;//Everything above is padded
double kpad_down;//Everything below is padded
double kini;//first bin
double kfin;//last bin
int N=128;//Number of points sampled (including padding areas)

//this seems to work fine for kmax=0.30 and kmin=0.02
kpad_up=0.5;
kpad_down=0.001;
kini=0.001/20.;
kfin=0.5*2;

params[0]=kpad_up;
params[1]=kpad_down;;
params[2]=N;
params[3]=kini;
params[4]=log(kfin/kini)/(N*1.-1.);

}

void apply_mask(char *type_of_analysis, int modeP0, int modeP2, int modeP4, double *k_theo,double P_theo0[], double P_theo2[], double P_theo4[],int Nlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *sampling_mask, fftw_plan plan1, fftw_plan plan2,double kmin, double kmax, double kmin_data0, double kmax_data0, double kmin_data2, double kmax_data2, double kmin_data4, double kmax_data4, char *sampling_theo, double Sigma_smooth, int isbao)
{
//isbao=0 for rsd fit, isbao=1 for baofit. Only relevant if FSBAOISO
int i,j;
double kmax_pad,kmin_pad,k0,slope,k,W0_int,W2_int,W4_int,W6_int,W8_int;
int Nlog;
double *klog,*Plog0, *Plog2, *Plog4;
double *Xi_mono,*Xi_mono_mask,*R,*P_mono_mask,*Knl0logmask;
double *Xi_quadru,*Xi_quadru_mask,*P_quadru_mask,*Knl2logmask;
double *Xi_hexadeca,*Xi_hexadeca_mask,*P_hexadeca_mask,*Knl4logmask;
double abs_kmin,abs_kmax;
double params[5];
int index0,index2,index4;
char sampling_log[2000];
int interpolation_order,shiftN,Ninterpol;
double w1,w2,w0;

interpolation_order=1;

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}

sprintf(sampling_log, "log");

if(modeP0==1){abs_kmin=kmin_data0;abs_kmax=kmax_data0;}
if(modeP2==1){abs_kmin=kmin_data2;abs_kmax=kmax_data2;}
if(modeP4==1){abs_kmin=kmin_data4;abs_kmax=kmax_data4;}

if(modeP0==1 && abs_kmin>kmin_data0){abs_kmin=kmin_data0;}
if(modeP2==1 && abs_kmin>kmin_data2){abs_kmin=kmin_data2;}
if(modeP4==1 && abs_kmin>kmin_data4){abs_kmin=kmin_data4;}

if(modeP0==1 && abs_kmax<kmax_data0){abs_kmax=kmax_data0;}
if(modeP2==1 && abs_kmax<kmax_data2){abs_kmax=kmax_data2;}
if(modeP4==1 && abs_kmax<kmax_data4){abs_kmax=kmax_data4;}

set_mask_params(params,kmin,kmax,Nlin,abs_kmin,abs_kmax);

kmax_pad=params[0];
kmin_pad=params[1];
Nlog=(int)(params[2]);
k0=params[3];
slope=params[4];

klog =  (double*) calloc( Nlog, sizeof(double));
if(modeP0==1){Plog0 =  (double*) calloc( Nlog, sizeof(double));}
if(modeP2==1){Plog2 =  (double*) calloc( Nlog, sizeof(double));}
if(modeP4==1){Plog4 =  (double*) calloc( Nlog, sizeof(double));}

index0=0;
index2=0;
index4=0;
for(i=0;i<Nlog;i++)
{
   k=k0*exp(i*slope);
   klog[i]=k;
   if(k<kmin_pad || k>kmax_pad)//padding
   {
      if(modeP0==1){Plog0[i]=0;}
      if(modeP2==1){Plog2[i]=0;}
      if(modeP4==1){Plog4[i]=0;}

   }
   else
   {

Ninterpol=determine_N_singlearray(k_theo,k,Nlin,sampling_theo);

if(Ninterpol>=Nlin-shiftN || Ninterpol<0  || k<=0)
{
      if(modeP0==1){Plog0[i]=0;}
      if(modeP2==1){Plog2[i]=0;}
      if(modeP4==1){Plog4[i]=0;}
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k,Ninterpol,sampling_theo);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k,Ninterpol,sampling_theo);
w1=determine_w1_2ndorder_singlearray(k_theo,k,Ninterpol,sampling_theo);
w2=determine_w2_2ndorder_singlearray(k_theo,k,Ninterpol,sampling_theo);
}

      if(modeP0==1){Plog0[i]=P_interpol_fast(k,P_theo0,Nlin,sampling_theo,interpolation_order,Ninterpol,w0,w1,w2);}
      if(modeP2==1){Plog2[i]=P_interpol_fast(k,P_theo2,Nlin,sampling_theo,interpolation_order,Ninterpol,w0,w1,w2);}
      if(modeP4==1){Plog4[i]=P_interpol_fast(k,P_theo4,Nlin,sampling_theo,interpolation_order,Ninterpol,w0,w1,w2);}
}

/*
      if(modeP0==1){Plog0[i]=P_interpol(k,k_theo,P_theo0,Nlin);}
      if(modeP2==1){Plog2[i]=P_interpol(k,k_theo,P_theo2,Nlin);}
      if(modeP4==1){Plog4[i]=P_interpol(k,k_theo,P_theo4,Nlin);}
*/
      if(strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0)
      {
       if(modeP0==1 && modeP2==1 && isbao==1){Plog2[i]=5./2.*(Plog2[i]-Plog0[i]);index2=2;}
       if(modeP0==1 && modeP2==1 && modeP4==1 && isbao==1){Plog4[i]=63./8.*(Plog4[i]-Plog0[i]-4./7.*Plog2[i]);index4=4;}
      }
   }
}


      if(strcmp(type_of_analysis, "FS") == 0 || strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){index2=2;index4=4;}
      if(strcmp(type_of_analysis, "FSBAOISO") == 0 && isbao==0){index2=2;index4=4;}

//for(i=0;i<Nlin;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i], P_theo2[i]);}exit(0);

//k_theo,double P_theo0[], double P_theo2[], double P_theo4[],int Nlin
//printf("Nlog=%d\n",Nlog);
R = (double*) calloc( Nlog, sizeof(double));

if(modeP0==1){
Xi_mono = (double*) calloc( Nlog, sizeof(double));
Xi_mono_mask = (double*) calloc( Nlog, sizeof(double));

fftlog_ComputeXiLM_threadsafe(index0, 2, Nlog, klog, Plog0, R, Xi_mono,plan1,plan2);
}
//exit(0);


if(modeP2==1){
Xi_quadru = (double*) calloc( Nlog, sizeof(double));
Xi_quadru_mask = (double*) calloc( Nlog, sizeof(double));
fftlog_ComputeXiLM_threadsafe(index2, 2, Nlog, klog, Plog2, R, Xi_quadru,plan1,plan2);

}

if(modeP4==1){
Xi_hexadeca = (double*) calloc( Nlog, sizeof(double));
Xi_hexadeca_mask = (double*) calloc( Nlog, sizeof(double));
fftlog_ComputeXiLM_threadsafe(index4, 2, Nlog, klog, Plog4, R, Xi_hexadeca,plan1,plan2);

}

free(klog);
if(modeP0==1){free(Plog0);}
if(modeP2==1){free(Plog2);}
if(modeP4==1){free(Plog4);}

for(i=0;i<Nlog;i++)
{

if(R[i]<=pos[0]){

W0_int=W0[0];
W2_int=W2[0];
W4_int=W4[0];
W6_int=W6[0];
W8_int=W8[0];

}
else{

/*
Ninterpol=determine_N_singlearray(pos,R[i],Nmask,sampling_mask);
if(Ninterpol>=Nmask-shiftN || Ninterpol<0  || R[i]<=0)
{
W0_int=0;
W2_int=0;
W4_int=0;
W6_int=0;
W8_int=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(pos,R[i],Ninterpol,sampling_mask);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(pos,R[j],Ninterpol,sampling_mask);
w1=determine_w1_2ndorder_singlearray(pos,R[j],Ninterpol,sampling_mask);
w2=determine_w2_2ndorder_singlearray(pos,R[j],Ninterpol,sampling_mask);
}


W0_int=P_interpol_fast(R[i],W0,Nmask,sampling_mask,interpolation_order,Ninterpol,w0,w1,w2);
W2_int=P_interpol_fast(R[i],W2,Nmask,sampling_mask,interpolation_order,Ninterpol,w0,w1,w2);
W4_int=P_interpol_fast(R[i],W4,Nmask,sampling_mask,interpolation_order,Ninterpol,w0,w1,w2);
W6_int=P_interpol_fast(R[i],W6,Nmask,sampling_mask,interpolation_order,Ninterpol,w0,w1,w2);
W8_int=P_interpol_fast(R[i],W8,Nmask,sampling_mask,interpolation_order,Ninterpol,w0,w1,w2);

}
*/

//printf("A %lf %lf %lf %lf %lf\n",W0_int,W2_int,W4_int,W6_int,W8_int);
W0_int=P_interpol( R[i] , pos, W0, Nmask);
W2_int=P_interpol( R[i] , pos, W2, Nmask);
W4_int=P_interpol( R[i] , pos, W4, Nmask);
W6_int=P_interpol( R[i] , pos, W6, Nmask);
W8_int=P_interpol( R[i] , pos, W8, Nmask);
//printf("B %lf %lf %lf %lf %lf\n",W0_int,W2_int,W4_int,W6_int,W8_int);


}



if(modeP0==1 && modeP2==1 && modeP4==1)//cas 111
{
Xi_mono_mask[i]=W0_int*Xi_mono[i]+1./5.*W2_int*Xi_quadru[i]+1./9.*W4_int*Xi_hexadeca[i];
Xi_quadru_mask[i]=W2_int*Xi_mono[i]+Xi_quadru[i]*(W0_int+2./7.*W2_int+2./7.*W4_int)+Xi_hexadeca[i]*(2./7.*W2_int+100./693.*W4_int+25./143.*W6_int);
Xi_hexadeca_mask[i]=W4_int*Xi_mono[i]+Xi_quadru[i]*(18./35.*W2_int+20./77.*W4_int+45./143.*W6_int)+Xi_hexadeca[i]*(W0_int+20./77.*W2_int+162./1001.*W4_int+20./143.*W6_int+490./2431.*W8_int);
}

if(modeP0==1 && modeP2==0 && modeP4==0)//100
{
Xi_mono_mask[i]=W0_int*Xi_mono[i];
}

if(modeP0==0 && modeP2==1 && modeP4==0)//010
{
Xi_quadru_mask[i]=Xi_quadru[i]*(W0_int+2./7.*W2_int+2./7.*W4_int);
}

if(modeP0==0 && modeP2==0 && modeP4==1)//001
{
Xi_hexadeca_mask[i]=Xi_hexadeca[i]*(W0_int+20./77.*W2_int+162./1001.*W4_int+20./143.*W6_int+490./2431.*W8_int);
}

if(modeP0==1 && modeP2==1 && modeP4==0)//110
{
Xi_mono_mask[i]=W0_int*Xi_mono[i]+1./5.*W2_int*Xi_quadru[i];
Xi_quadru_mask[i]=W2_int*Xi_mono[i]+Xi_quadru[i]*(W0_int+2./7.*W2_int+2./7.*W4_int);
}

if(modeP0==1 && modeP2==0 && modeP4==1)//101
{
Xi_mono_mask[i]=W0_int*Xi_mono[i]+1./9.*W4_int*Xi_hexadeca[i];
Xi_hexadeca_mask[i]=W4_int*Xi_mono[i]+Xi_hexadeca[i]*(W0_int+20./77.*W2_int+162./1001.*W4_int+20./143.*W6_int+490./2431.*W8_int);
}

if(modeP0==0 && modeP2==1 && modeP4==1)//011
{
Xi_quadru_mask[i]=Xi_quadru[i]*(W0_int+2./7.*W2_int+2./7.*W4_int)+Xi_hexadeca[i]*(2./7.*W2_int+100./693.*W4_int+25./143.*W6_int);
Xi_hexadeca_mask[i]=Xi_quadru[i]*(18./35.*W2_int+20./77.*W4_int+45./143.*W6_int)+Xi_hexadeca[i]*(W0_int+20./77.*W2_int+162./1001.*W4_int+20./143.*W6_int+490./2431.*W8_int);
}



}

if(modeP0==1){free(Xi_mono);}
if(modeP2==1){free(Xi_quadru);}
if(modeP4==1){free(Xi_hexadeca);}

if(modeP0==1){
P_mono_mask = (double*) calloc( Nlog, sizeof(double));
Knl0logmask = (double*) calloc( Nlog, sizeof(double));
fftlog_ComputeXiLM_threadsafe(index0, 2, Nlog, R, Xi_mono_mask, Knl0logmask, P_mono_mask,plan1,plan2);
free(Xi_mono_mask);
}

if(modeP2==1){
P_quadru_mask = (double*) calloc( Nlog, sizeof(double));
Knl2logmask = (double*) calloc( Nlog, sizeof(double));

fftlog_ComputeXiLM_threadsafe(index2, 2, Nlog, R, Xi_quadru_mask, Knl2logmask, P_quadru_mask,plan1,plan2);

free(Xi_quadru_mask);
}


if(modeP4==1){
P_hexadeca_mask = (double*) calloc( Nlog, sizeof(double));
Knl4logmask = (double*) calloc( Nlog, sizeof(double));


fftlog_ComputeXiLM_threadsafe(index4, 2, Nlog, R, Xi_hexadeca_mask, Knl4logmask, P_hexadeca_mask,plan1,plan2);

free(Xi_hexadeca_mask);
}

free(R);
for(i=0;i<Nlin;i++)
{

if(modeP0==1){

Ninterpol=determine_N_singlearray(Knl0logmask,k_theo[i],Nlog,sampling_log);

if(Ninterpol>=Nlog-shiftN || Ninterpol<0  || k_theo[i]<=0)
{
P_theo0[i]=0;
}else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(Knl0logmask,k_theo[i],Ninterpol,sampling_log);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(Knl0logmask,k_theo[i],Ninterpol,sampling_log);
w1=determine_w1_2ndorder_singlearray(Knl0logmask,k_theo[i],Ninterpol,sampling_log);
w2=determine_w2_2ndorder_singlearray(Knl0logmask,k_theo[i],Ninterpol,sampling_log);
}

P_theo0[i]=P_interpol_fast(k_theo[i],P_mono_mask,Nlog,sampling_log,interpolation_order,Ninterpol,w0,w1,w2)*pow(2.*Pi,3);

}

//P_theo0[i]=P_interpol(k_theo[i],Knl0logmask,P_mono_mask,Nlog)*pow(2.*Pi,3);
//P_theo0[i]=P_interpolLOG(k_theo[i],Knl0logmask,P_mono_mask,Nlog)*pow(2.*Pi,3);

}

if(modeP2==1){

Ninterpol=determine_N_singlearray(Knl2logmask,k_theo[i],Nlog,sampling_log);
if(Ninterpol>=Nlog-shiftN || Ninterpol<0  || k_theo[i]<=0)
{
P_theo2[i]=0;
}else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(Knl2logmask,k_theo[i],Ninterpol,sampling_log);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(Knl2logmask,k_theo[i],Ninterpol,sampling_log);
w1=determine_w1_2ndorder_singlearray(Knl2logmask,k_theo[i],Ninterpol,sampling_log);
w2=determine_w2_2ndorder_singlearray(Knl2logmask,k_theo[i],Ninterpol,sampling_log);
}

P_theo2[i]=P_interpol_fast(k_theo[i],P_quadru_mask,Nlog,sampling_log,interpolation_order,Ninterpol,w0,w1,w2)*pow(2.*Pi,3);

}

//P_theo2[i]=P_interpol(k_theo[i],Knl2logmask,P_quadru_mask,Nlog)*pow(2.*Pi,3);


if(strcmp(type_of_analysis, "BAOISO") == 0 && modeP0==1){P_theo2[i]=P_theo0[i]+2./5.*P_theo2[i];}
if(strcmp(type_of_analysis, "FSBAOISO") == 0 && modeP0==1 && isbao==1){P_theo2[i]=P_theo0[i]+2./5.*P_theo2[i];}
}

if(modeP4==1){

Ninterpol=determine_N_singlearray(Knl4logmask,k_theo[i],Nlog,sampling_log);
if(Ninterpol>=Nlog-shiftN || Ninterpol<0  || k_theo[i]<=0)
{
P_theo4[i]=0;
}else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(Knl4logmask,k_theo[i],Ninterpol,sampling_log);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(Knl4logmask,k_theo[i],Ninterpol,sampling_log);
w1=determine_w1_2ndorder_singlearray(Knl4logmask,k_theo[i],Ninterpol,sampling_log);
w2=determine_w2_2ndorder_singlearray(Knl4logmask,k_theo[i],Ninterpol,sampling_log);
}
P_theo4[i]=P_interpol_fast(k_theo[i],P_hexadeca_mask,Nlog,sampling_log,interpolation_order,Ninterpol,w0,w1,w2)*pow(2.*Pi,3);



}

//P_theo4[i]=P_interpol(k_theo[i],Knl4logmask,P_hexadeca_mask,Nlog)*pow(2.*Pi,3);

if(strcmp(type_of_analysis, "BAOISO") == 0 && modeP0==1 && modeP2==1){P_theo4[i]=P_theo0[i]+(4./7.)*(5./2.)*(P_theo2[i]-P_theo0[i])+8./63.*P_theo4[i];}
if(strcmp(type_of_analysis, "FSBAOISO") == 0 && modeP0==1 && modeP2==1 && isbao==1){P_theo4[i]=P_theo0[i]+(4./7.)*(5./2.)*(P_theo2[i]-P_theo0[i])+8./63.*P_theo4[i];}
}

}

if(modeP0==1){
free(P_mono_mask);
free(Knl0logmask);
}

if(modeP2==1){
free(P_quadru_mask);
free(Knl2logmask);
}

if(modeP4==1){
free(P_hexadeca_mask);
free(Knl4logmask);
}


//exit(0);
}


double gauss(double x, double Sigma_nl_mean,double Sigma_nl_stddev)
{
double f;

f=pow(x-Sigma_nl_mean,2)/(Sigma_nl_stddev*Sigma_nl_stddev);
return f;

}

void do_log_file2(char *path_output,char *identifier,char *likelihood_file,char *fit_BAO,int Nchunks,int Npolynomial,double time_run,int processors,double amin,double amax,double astep,int Nalphas,int Nsigmas)
{

FILE *f;
long int Nlines,i,j,Nlines1,l1,l2,l1min,l2min,itrial;
int junk;
char log_file[200];
double *alpha1,*alpha2,*chi2_value,errorA_alpha1,errorB_alpha1,errorA_alpha2,errorB_alpha2;
double *alpha11,*alpha22,*chi2_value11,*chi2_value22;
int modeP0,modeP2,modeP4;
long int imin;
double chi2min,chi2min1,chi2min2;
double alpha1mean,alpha2mean,error1,error2;
int detection_alpha1;
int detection_alpha2;
sprintf(log_file,"%s/logfileanalytic_%s.txt",path_output,identifier);

modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){modeP4=1;}

Nlines=countlinesLI(likelihood_file)-2;

alpha1=  (double*) calloc( Nlines, sizeof(double));
if(modeP0+modeP2+modeP4>1){
alpha2= (double*) calloc( Nlines, sizeof(double));

Nlines1=(long int)(sqrt(Nlines*1.));
alpha11= (double*) calloc( Nlines1, sizeof(double));
alpha22= (double*) calloc( Nlines1, sizeof(double));
chi2_value11= (double*) calloc( Nlines1, sizeof(double));
chi2_value22= (double*) calloc( Nlines1, sizeof(double));


}
chi2_value= (double*) calloc( Nlines, sizeof(double));

junk=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+Nalphas+Nsigmas+1;
//if(modeP0+modeP2+modeP4>1){junk++;}
//printf("junk=%d\n",junk);
f=fopen(likelihood_file,"r");
for(j=0;j<junk;j++){
if(j==junk-1){fscanf(f,"%*s\n");}
else{fscanf(f,"%*s ");}
}

fscanf(f,"%*s\n");
l1=0;l2=0;
for(i=0;i<Nlines;i++)
{

//if(modeP0+modeP2+modeP4>1){for(j=0;j<junk-3;j++){fscanf(f,"%*f ");}}
//if(modeP0+modeP2+modeP4==1){for(j=0;j<junk-2;j++){fscanf(f,"%*f ");}}
for(j=0;j<junk-Nalphas-1;j++){fscanf(f,"%*f ");}

if(modeP0+modeP2+modeP4>1){fscanf(f,"%lf %lf %lf\n",&alpha1[i],&alpha2[i],&chi2_value[i]);}
if(modeP0+modeP2+modeP4==1){fscanf(f,"%lf %lf\n",&alpha1[i],&chi2_value[i]);}

if(i==0){chi2min=chi2_value[i];imin=i;}
else{if(chi2min>chi2_value[i]){chi2min=chi2_value[i];imin=i;}}

if(modeP0+modeP2+modeP4>1)
{

   if(i==0){alpha11[l1]=alpha1[i];chi2_value11[l1]=chi2_value[i];}
   else{

         if(alpha1[i]==alpha11[l1] && chi2_value11[l1]>chi2_value[i]){chi2_value11[l1]=chi2_value[i];}
         if(alpha1[i]>alpha11[l1]){alpha11[l1+1]=alpha1[i];chi2_value11[l1+1]=chi2_value[i];l1++;}
       }

         if(i<Nlines1)
         {
             alpha22[i]=alpha2[i];
             chi2_value22[i]=chi2_value[i];

         }
         else
         {
              itrial=(long int)(i*1./Nlines1*1.);
              l2=i-Nlines1*itrial;
              if(chi2_value[i]<chi2_value22[l2] && alpha2[i]==alpha22[l2]){chi2_value22[l2]=chi2_value[i];}

         }

     

}

}


fclose(f);


if(modeP0+modeP2+modeP4==1){

errorA_alpha1=get_error1(alpha1,chi2_value,chi2min,imin,Nlines);
errorB_alpha1=get_error2(alpha1,chi2_value,chi2min,imin,Nlines);
}

if(modeP0+modeP2+modeP4>1){

for(i=0;i<Nlines1;i++)
{

if(i==0){chi2min1=chi2_value11[i];l1min=i;}
else{if(chi2min1>chi2_value11[i]){chi2min1=chi2_value11[i];l1min=i;}}

if(i==0){chi2min2=chi2_value22[i];l2min=i;}
else{if(chi2min2>chi2_value22[i]){chi2min2=chi2_value22[i];l2min=i;}}
}

//printf("%lf %d %d\n",chi2min,l1min,Nlines1);
errorA_alpha1=get_error1(alpha11,chi2_value11,chi2min,l1min,Nlines1);
errorB_alpha1=get_error2(alpha11,chi2_value11,chi2min,l1min,Nlines1);
//exit(0);
errorA_alpha2=get_error1(alpha22,chi2_value22,chi2min,l2min,Nlines1);
errorB_alpha2=get_error2(alpha22,chi2_value22,chi2min,l2min,Nlines1);

}

f=fopen(log_file,"w");

fprintf(f,"Range of alphas %lf - %lf\n",amin,amax);
fprintf(f,"Step in alpha %lf\n\n",astep);

detection_alpha1=1;
if(errorB_alpha1==0 || errorA_alpha1==0){detection_alpha1=0;}
//printf("%lf %lf\n",errorB_alpha1,errorA_alpha1);

//for(i=0;i<Nlines1;i++){printf("%d %lf %lf\n",i,alpha11[i],chi2_value11[i]);}


fprintf(f,"Mean parameters\n");
if(detection_alpha1==1){

alpha1mean=alpha1[imin]-errorA_alpha1+0.5*(errorA_alpha1+errorB_alpha1);
error1=0.5*(errorA_alpha1+errorB_alpha1);

if(modeP0==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha0= ");}
if(modeP2==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha2= ");}
if(modeP4==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha4= ");}
if(modeP0+modeP2+modeP4>1){fprintf(f,"alpha_para= ");}

fprintf(f,"%lf pm %lf\n",alpha1mean,error1);
}
else
{

if(modeP0==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha0= ");}
if(modeP2==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha2= ");}
if(modeP4==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha4= ");}
if(modeP0+modeP2+modeP4>1){fprintf(f,"alpha_para= ");}

fprintf(f,"not detected in this range\n");

}


if(modeP0+modeP2+modeP4>1){

detection_alpha2=1;
if(errorB_alpha2==0 || errorA_alpha2==0){detection_alpha2=0;}

if(detection_alpha2==1){

alpha2mean=alpha2[imin]-errorA_alpha2+0.5*(errorA_alpha2+errorB_alpha2);
error2=0.5*(errorA_alpha2+errorB_alpha2);

fprintf(f,"alpha_perp= %lf pm %lf\n",alpha2mean,error2);
}
else
{
fprintf(f,"alpha_perp= not detected in this range\n");
}

}




fprintf(f,"\nBest-fitting parameters\n");
if(modeP0==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha0= ");}
if(modeP2==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha2= ");}
if(modeP4==1 && modeP0+modeP2+modeP4==1){fprintf(f,"alpha4= ");}
if(modeP0+modeP2+modeP4>1){fprintf(f,"alpha_para= ");}
fprintf(f,"%lf\n",alpha1[imin]);
if(modeP0+modeP2+modeP4>1){fprintf(f,"alpha_perp= %lf\n",alpha2[imin]);}
if(modeP0+modeP2+modeP4>1){fprintf(f,"rho= XX\n");} 
fprintf(f,"chi2min= %lf\n",chi2min);


fprintf(f,"\nTime take for analysis to run");

if(time_run>24.)
{
time_run=time_run/24.;
fprintf(f," %lf days",time_run);
}
else
{

if(time_run>1){
fprintf(f," %lf hours",time_run);
}
else
{

if(time_run*60<1)
{
fprintf(f," %lf seconds",time_run*60*60);
}
else{
fprintf(f," %lf minutes",time_run*60);
}

}

}
fprintf(f," using %d processors.",processors);

fclose(f);


free(alpha1);
if(modeP0+modeP2+modeP4>1){
free(alpha2);

free(alpha11);
free(alpha22);
free(chi2_value11);
free(chi2_value22);
}
free(chi2_value);
}

double hi0_mask(int index,int mode,int modeP0,int modeP2,int modeP4, double ki, double alpha_0,double alpha_00, double Sigma_nl,int mode_parameter, double *pos, double *W0,double *W2,double *W4,double *W6,double *W8, int Nmask, char *path_to_mask, double *klin, double *Plin, int Nlin, double *kolin, double *Polin, int Nolin,fftw_plan plan1, fftw_plan plan2,double kmin,double kmax,double kmin_data,double kmax_data)
{
double f;

int i,j;
double kmax_pad,kmin_pad,k0,slope,k,W0_int,W2_int,W4_int,W6_int,W8_int;
int Nlog;
double *klog,*Plog,*Xi_mono,*Xi_mono_mask,*R0,*P_mono_mask,*Knl0logmask;
double alpha_in;
double olin_eff;
double func;

double params[5];
set_mask_params(params,kmin,kmax,Nlin,kmin_data,kmax_data);

kmax_pad=params[0];
kmin_pad=params[1];
Nlog=(int)(params[2]);
k0=params[3];
slope=params[4];

if(modeP0+modeP2+modeP4>1){
if(mode==0){alpha_in=pow(alpha_0,1./3.)*pow(alpha_00,2./3.);}
if(mode==2){alpha_in=pow(alpha_0,3./5.)*pow(alpha_00,2./5.);}
if(mode==4){alpha_in=pow(alpha_0,5./7.)*pow(alpha_00,2./7.);}
}
else{
alpha_in=alpha_0;
}

klog =  (double*) calloc( Nlog, sizeof(double));
Plog =  (double*) calloc( Nlog, sizeof(double));

for(i=0;i<Nlog;i++)
{
   k=k0*exp(i*slope);
   klog[i]=k;
   if(k<kmin_pad || k>kmax_pad)//padding
   {
      Plog[i]=0;
   }
   else
   {

       if(mode_parameter==1)
       {
          func=P_interpol(k,klin,Plin,Nlin);
       }
       else
       {
             func=pow(k,3-mode_parameter);
       }

          olin_eff=P_interpol(k/alpha_in,kolin,Polin,Nolin);
          Plog[i]=func*(1.+(olin_eff-1)*exp(-0.5*Sigma_nl*Sigma_nl*k*k));
   }
}

Xi_mono = (double*) calloc( Nlog, sizeof(double));
Xi_mono_mask = (double*) calloc( Nlog, sizeof(double));
R0 = (double*) calloc( Nlog, sizeof(double));

fftlog_ComputeXiLM_threadsafe(index, 2, Nlog, klog, Plog, R0, Xi_mono,plan1,plan2);
free(klog);
free(Plog);

for(i=0;i<Nlog;i++)
{
if(R0[i]<=pos[0]){
W0_int=W0[0];
W2_int=W2[0];
W4_int=W4[0];
W6_int=W6[0];
W8_int=W8[0];
}
else{
W0_int=P_interpol( R0[i] , pos, W0, Nmask);
W2_int=P_interpol( R0[i] , pos, W2, Nmask);
W4_int=P_interpol( R0[i] , pos, W4, Nmask);
W6_int=P_interpol( R0[i] , pos, W6, Nmask);
W8_int=P_interpol( R0[i] , pos, W8, Nmask);
}

if(mode==0){Xi_mono_mask[i]=W0_int*Xi_mono[i];}
if(mode==2){Xi_mono_mask[i]=(W0_int+2./7.*W2_int+2./7.*W4_int)*Xi_mono[i];}
if(mode==4){Xi_mono_mask[i]=(W0_int+20./77.*W2_int+162./1001.*W4_int+20./143.*W6_int+490./2431.*W8_int)*Xi_mono[i];}

}

free(Xi_mono);

P_mono_mask = (double*) calloc( Nlog, sizeof(double));
Knl0logmask = (double*) calloc( Nlog, sizeof(double));
fftlog_ComputeXiLM_threadsafe(index, 2, Nlog, R0, Xi_mono_mask, Knl0logmask, P_mono_mask,plan1,plan2);
free(Xi_mono_mask);
free(R0);

f=P_interpol(ki,Knl0logmask,P_mono_mask,Nlog)*pow(2.*Pi,3);
if(f!=f){printf("Error with in hi0_mask function, f=%lf\n. Exiting now...\n",f);exit(0);}
free(P_mono_mask);
free(Knl0logmask);


return f;
}


double hi0(int mode,int modeP0,int modeP2,int modeP4, double ki, double alpha_0,double alpha_00, double Sigma_nl,int mode_parameter, double *klin, double *Plin, int Nlin, double *kolin, double *Polin, int Nolin)
{
int i;
double f;
double olin_eff;
double func;
double alpha_in;
if(modeP0+modeP2+modeP4>1){
if(mode==0){alpha_in=pow(alpha_0,1./3.)*pow(alpha_00,2./3.);}
if(mode==2){alpha_in=pow(alpha_0,3./5.)*pow(alpha_00,2./5.);}
if(mode==4){alpha_in=pow(alpha_0,5./7.)*pow(alpha_00,2./7.);}
}
else{
alpha_in=alpha_0;
}


    if(mode_parameter==1)
    {
      func=P_interpol(ki,klin,Plin,Nlin);
    }
    else
    {
       func=pow(ki,3-mode_parameter);
    }

olin_eff=P_interpol(ki/alpha_in,kolin,Polin,Nolin);
f=func*(1.+(olin_eff-1)*exp(-0.5*Sigma_nl*Sigma_nl*ki*ki));

return f;
}




void parameter_bestfit(int mode, int modeP0,int modeP2,int modeP4,double alpha_0,double alpha_00,double Sigma_nl, double *Knl, double *Pk,double *Cov,int Ncov, double A[],int Npoly, double *pos, double *W0,double *W2, double *W4, double *W6,double *W8, int Nmask, char *path_to_mask, double *klin, double *Plin, int Nlin, double *kolin, double *Polin, int Nolin,fftw_plan plan1, fftw_plan plan2)
{
        
        int i,j;
        int l,lp;
        double cov;
        int n_params=Npoly+1;
        int s;
        double input;
        double hl,hlp;
        double pk;
        double *B,*H,*Hp;
        int mask;//0 no mask, 1 mask;

//in the case of P2<<P0 the Pmu-windowed is dominated by the signal of the monopole

if(strcmp(path_to_mask, "none") == 0){mask=0;}else{mask=1;}
 
        gsl_matrix * m = gsl_matrix_alloc (n_params, n_params);
        gsl_matrix * inverse = gsl_matrix_alloc (n_params, n_params);
        gsl_matrix * identity = gsl_matrix_alloc (n_params, n_params);

        gsl_permutation * perm = gsl_permutation_alloc (n_params);



         for(i=0;i<n_params;i++)
         {
         for(j=i;j<n_params;j++)
         {
            input=0;
            for(l=0;l<Ncov;l++)
            {
                   if(mask==0){ hl=hi0(mode,modeP0,modeP2,modeP4,Knl[l], alpha_0,alpha_00, Sigma_nl,i+1, klin, Plin, Nlin, kolin, Polin, Nolin);}
                   if(mask==1){ hl=hi0_mask(0,mode,modeP0,modeP2,modeP4,Knl[l], alpha_0,alpha_00, Sigma_nl,i+1,pos,W0,W2,W4,W6,W8,Nmask,path_to_mask, klin, Plin, Nlin, kolin, Polin, Nolin,plan1,plan2,klin[0],klin[Nlin-1],Knl[0],Knl[Ncov-1]);}


                    for(lp=0;lp<Ncov;lp++)
                    {
                            if(mask==0){hlp=hi0(mode,modeP0,modeP2,modeP4,Knl[lp], alpha_0,alpha_00, Sigma_nl,j+1,  klin, Plin, Nlin, kolin, Polin, Nolin);}
                            if(mask==1){hlp=hi0_mask(0,mode,modeP0,modeP2,modeP4,Knl[lp], alpha_0,alpha_00, Sigma_nl,j+1,pos,W0,W2,W4,W6,W8,Nmask,path_to_mask,  klin, Plin, Nlin, kolin, Polin, Nolin,plan1,plan2,klin[0],klin[Nlin-1],Knl[0],Knl[Ncov-1]);}

                            cov=Cov[l+lp*(Ncov)];                             
                            input=input+hl*(1./cov)*hlp;
                   }
              }

              gsl_matrix_set (m, i, j, input);
              if(i!=j){gsl_matrix_set (m, j, i, input);}
        }
        }
           gsl_linalg_LU_decomp (m, perm, &s);
           gsl_linalg_LU_invert (m, perm, inverse);

              B = (double*) calloc( n_params, sizeof(double));
              
              H = (double*) calloc( n_params, sizeof(double));
          
              Hp = (double*) calloc( n_params, sizeof(double));


            for(l=0;l<Ncov;l++)
            {

              for(i=0;i<n_params;i++){
if(mask==0){Hp[i]=hi0(mode,modeP0,modeP2,modeP4,Knl[l], alpha_0,alpha_00, Sigma_nl,i+1, klin, Plin, Nlin, kolin, Polin, Nolin);}
if(mask==1){Hp[i]=hi0_mask(0,mode,modeP0,modeP2,modeP4,Knl[l], alpha_0,alpha_00, Sigma_nl,i+1, pos,W0,W2,W4,W6,W8,Nmask,path_to_mask, klin, Plin, Nlin, kolin, Polin, Nolin,plan1,plan2,klin[0],klin[Nlin-1],Knl[0],Knl[Ncov-1]);}

}

              for(lp=0;lp<Ncov;lp++)
              {
                 pk=Pk[lp];
                 cov=Cov[l+lp*(Ncov)];

                 for(j=0;j<n_params;j++)
                 {
                     for(i=0;i<n_params;i++){H[i]=Hp[i]*gsl_matrix_get(inverse, j, i)*pk;}
                     for(i=0;i<n_params;i++){B[j]=B[j]+(1./cov)*H[i];}
                 }
              }

            }

for(i=0;i<Npoly+1;i++){A[i]=B[i];}

free(B);
free(H);
free(Hp);
//free matrices
gsl_matrix_free(m);
gsl_matrix_free(inverse);
gsl_matrix_free(identity);
gsl_permutation_free(perm);

}



void  make_a_bao_plot(char *type_BAO_fit, char *type_of_analysis, char *fit_BAO,double *parameters2, double chi2_min, double *k0, double *P0, double *errP0,int NeffP0, double *k0SGC, double *P0SGC, double *errP0SGC, int NeffP0SGC, double *k2, double *P2, double *errP2,int NeffP2, double *k2SGC, double *P2SGC, double *errP2SGC, int NeffP2SGC, double *k4, double *P4, double *errP4, int NeffP4, double *k4SGC, double *P4SGC, double *errP4SGC, int NeffP4SGC, double *k_Plin, double *Plin, int Nlin, double *k_Olin, double *Olin, int NOlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC,double *W4SGC, double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC,  char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks,int Npoints,int Ndof, char *path_output, char *identifier,fftw_plan plan1, fftw_plan plan2,int Nalphas,int Nsigmas_tot, int Nsigmas_free, double Sigma_smooth,int factor_sampling_mask_in,char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory)
{
FILE *f1,*f2,*f12;
int i,j,l,open,NeffP;
char name1[2000],name2[2000],name12[2000];
char nameP0_1[2000],nameP0_2[2000],nameP0_12[2000];
char nameP2_1[2000],nameP2_2[2000],nameP2_12[2000];
char nameP4_1[2000],nameP4_2[2000],nameP4_12[2000];

int combine;
double ptheo,ptheoSGC,ptheotot;
double pnw,pnwSGC,pnwtot;
double P0tot,errP0tot;
combine=0;//0 yes, otherwise no
int Ntheo;
int modeP0,modeP2,modeP4;
//int dimension;
double *ktheo,*ksm,*Ptheo0,*ksm0,*Psm0,*ksm2,*Psm2,*ksm4,*Psm4;
double *Ptheo2;
double *Ptheo4;

double *ktheo0,*ktheo2,*ktheo4;

double *ktheoSGC,*ksmSGC,*Ptheo0SGC, *Psm0SGC;
double *Ptheo2SGC,*Psm2SGC;
double *Ptheo4SGC, *Psm4SGC;
double *ktheo0SGC,*ktheo2SGC,*ktheo4SGC,*ksm0SGC,*ksm2SGC,*ksm4SGC;

double *parameters1;
int ik01_P0,ik01_P2,ik01_P4;
double difference;
double difference_min;Ntheo=Nlin;
double kmin0,kmin2,kmin4;
double kmax0,kmax2,kmax4;
double kmin0SGC,kmin2SGC,kmin4SGC;
double kmax0SGC,kmax2SGC,kmax4SGC;
int offset,offset_ini;
double Sigmanl0,Sigmanl2,Sigmanl4;
int Nalphas_for_param1,Nsigmas_for_param1;
int Neffmax,NeffmaxSGC;
int factor_sampling_mask;
int interpolation_order,shiftN;
double w1,w2,w0;
int Ninterpol;
double w1sm,w2sm,w0sm;
int Ninterpolsm;
interpolation_order=1;

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}

//This part was missing, but in the plotting....
modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}


factor_sampling_mask=factor_sampling_mask_in;

kmin0=0;kmax0=0;
kmin2=0;kmax2=0;
kmin4=0;kmax4=0;

if(strcmp(path_to_mask1, "none") != 0 ){

/*
if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
*/

if( strcmp(spacing_dataNGC,"linear") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}

}

if( strcmp(spacing_dataNGC,"log") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(log(k0[0])/(log(k0[NeffP0-1])-log(k0[0]))*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(log(k0[0])/(log(k0[NeffP0-1])-log(k0[0]))*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(log(k2[0])/(log(k2[NeffP2-1])-log(k2[0]))*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(log(k2[0])/(log(k2[NeffP2-1])-log(k2[0]))*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(log(k4[0])/(log(k4[NeffP4-1])-log(k4[0]))*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(log(k4[0])/(log(k4[NeffP4-1])-log(k4[0]))*(NeffP4-1));}

}

if( strcmp(spacing_dataNGC,"log10") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(log10(k0[0])/(log10(k0[NeffP0-1])-log10(k0[0]))*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(log10(k0[0])/(log10(k0[NeffP0-1])-log10(k0[0]))*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(log10(k2[0])/(log10(k2[NeffP2-1])-log10(k2[0]))*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(log10(k2[0])/(log10(k2[NeffP2-1])-log10(k2[0]))*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(log10(k4[0])/(log10(k4[NeffP4-1])-log10(k4[0]))*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(log10(k4[0])/(log10(k4[NeffP4-1])-log10(k4[0]))*(NeffP4-1));}

}

if( strcmp(spacing_dataNGC,"irregular") == 0  ){
Neffmax=Ntheo;
factor_sampling_mask=1;
}



ktheo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
ksm = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}

modeP0=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
if(strcmp(path_to_mask1, "none") == 0 ){
ksm0 = (double*) calloc( NeffP0, sizeof(double));
ktheo0 = (double*) calloc( NeffP0, sizeof(double));
Ptheo0 = (double*) calloc( NeffP0, sizeof(double));
Psm0 = (double*) calloc( NeffP0, sizeof(double));
}
else{
Ptheo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
Psm0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}

modeP0=1;
kmin0=k0[0];
kmax0= k0[NeffP0-1];
}

modeP2=0;
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
if(strcmp(path_to_mask1, "none") == 0 ){
Ptheo2 = (double*) calloc( NeffP2, sizeof(double));
ktheo2 = (double*) calloc( NeffP2, sizeof(double));
Psm2 = (double*) calloc( NeffP2, sizeof(double));
ksm2 = (double*) calloc( NeffP2, sizeof(double));
}
else{
Ptheo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
Psm2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}
modeP2=1;
kmin2=k2[0];
kmax2=k2[NeffP2-1];
}

modeP4=0;
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
if(strcmp(path_to_mask1, "none") == 0 ){
Ptheo4 = (double*) calloc( NeffP4, sizeof(double));
ktheo4 = (double*) calloc( NeffP4, sizeof(double));
Psm4 = (double*) calloc( NeffP4, sizeof(double));
ksm4 = (double*) calloc( NeffP4, sizeof(double));
}
else{
Ptheo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
Psm4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}
modeP4=1;
kmin4=k4[0];
kmax4=k4[NeffP4-1];
}


//if(modeP0+modeP2+modeP4>1){dimension=3+(Npolynomial+1)*(modeP0+modeP2+modeP4);}
//else{dimension=2+(Npolynomial+1);}
//parameters1 = (double*) calloc( dimension, sizeof(double));
/*
if(modeP0+modeP2+modeP4==1){parameters1 =  (double*) calloc( Npolynomial+1+2, sizeof(double));}
if(modeP0+modeP2+modeP4>1){

      if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 ){   parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+5, sizeof(double));}
      else{

               if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0){parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+4, sizeof(double));}
               else{

                      if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0){parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+4, sizeof(double));}
                      else{
                             parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+3, sizeof(double));
                          }
                   }
         }
}
*/

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){

if(modeP0+modeP2+modeP4==1){Nalphas_for_param1=1;}
else{Nalphas_for_param1=2;}

Nsigmas_for_param1=modeP0+modeP2+modeP4;

parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+Nalphas_for_param1+Nsigmas_for_param1, sizeof(double));

//if(strcmp(type_BAO_fit, "analytic") == 0){
offset=Nalphas+Nsigmas_tot;
}


if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){

Nalphas_for_param1=2;
Nsigmas_for_param1=2;;
parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1, sizeof(double));

offset=Nalphas+Nsigmas_tot+1;

}



if(strcmp(Sigma_def_type, "effective") == 0)
{
if(modeP0+modeP2+modeP4==1){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];offset_ini=2;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];parameters1[4]=parameters2[2];offset_ini=5;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];parameters1[4]=parameters2[4];offset_ini=5;}


}
else
{
if(modeP0+modeP2+modeP4==1)
{

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[2],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[2],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[2],4./14.);}
}
else
{
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[1]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[1]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[1]/(1.+ffactor),4./14.);}
}

     parameters1[0]=parameters2[0];
if(modeP0==1){parameters1[1]=Sigmanl0;}
if(modeP2==1){parameters1[1]=Sigmanl2;}
if(modeP4==1){parameters1[1]=Sigmanl4;}
offset_ini=2;
}
else{
     parameters1[0]=parameters2[0];
     parameters1[1]=parameters2[1];

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[3],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[3],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[3],4./14.);}
}
else
{
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[2]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[2]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[2]/(1.+ffactor),4./14.);}
}


if(modeP0==0){parameters1[2]=Sigmanl2;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP2==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP4==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;offset_ini=4;}
if(modeP0+modeP2+modeP4==3){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;parameters1[4]=Sigmanl4;offset_ini=5;}

}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){

parameters1[2]=parameters2[2];
if(strcmp(Sigma_independent, "yes") == 0)
{
parameters1[3]=parameters2[3];
parameters1[4]=parameters2[4];
}
else
{
parameters1[3]=parameters2[2]/(1+ffactor);
parameters1[4]=parameters2[3];
}
offset_ini=5;

}


}

}

//}

//if(strcmp(type_BAO_fit, "mcmc") == 0)
//{


//}
//printf("%d %d\n",offset_ini,(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini);
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){for(i=offset_ini;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){parameters1[i]=parameters2[i-offset_ini+offset];}}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){for(i=offset_ini;i<1+(Npolynomial)*(modeP0+modeP2+modeP4)+offset_ini;i++){parameters1[i]=parameters2[i-offset_ini+offset];}}


//for(i=0;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){printf("%d %lf\n",i,parameters1[i]);}


//for(i=0;i<dimension;i++){parameters1[i]=parameters2[i];}

//for(i=0;i<20;i++){printf("NGC+SGC %d %lf\n",i,parameters2[i]);}

//for(i=0;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){printf("NGC %d %lf\n",i,parameters1[i]);}


if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){
do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ktheo,ktheo0,ktheo2,ktheo4, Ptheo0, Ptheo2, Ptheo4,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0,k2,k4, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],kmin0, kmax0,kmin2, kmax2,kmin4, kmax4, 1,spacing_dataNGC,spacing_theory,Sigma_smooth);

do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ksm,ksm0,ksm2,ksm4, Psm0, Psm2, Psm4,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC, path_to_mask1,k0,k2,k4, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1] ,kmin0 , kmax0,kmin2, kmax2,kmin4, kmax4, 0,spacing_dataNGC,spacing_theory,Sigma_smooth);
}


if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ktheo,ktheo0,ktheo2,ktheo4, Ptheo0, Ptheo2, Ptheo4,NeffP0,NeffP2,NeffP4, factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,Sigma_smooth,pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0,k2,k4, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],kmin0, kmax0,kmin2, kmax2,kmin4, kmax4, 1,spacing_dataNGC,spacing_theory);

do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ksm,ksm0,ksm2,ksm4, Psm0, Psm2, Psm4,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,Sigma_smooth,pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0,k2,k4, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1] ,kmin0 , kmax0,kmin2, kmax2,kmin4, kmax4, 0,spacing_dataNGC,spacing_theory);
}

if(Nchunks==2)
{
factor_sampling_mask=factor_sampling_mask_in;
kmin0SGC=0;kmax0SGC=0;
kmin2SGC=0;kmax2SGC=0;
kmin4SGC=0;kmax4SGC=0;

if(strcmp(path_to_mask2, "none") != 0 ){

/*
if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
*/


if( strcmp(spacing_dataSGC,"linear") == 0  ){

if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}

}

if( strcmp(spacing_dataSGC,"log") == 0  ){

if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP0SGC-2+25+(int)(log(k0SGC[0])/(log(k0SGC[NeffP0SGC-1])-log(k0SGC[0]))*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP0SGC-2+25+(int)(log(k0SGC[0])/(log(k0SGC[NeffP0SGC-1])-log(k0SGC[0]))*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP2SGC-2+25+(int)(log(k2SGC[0])/(log(k2SGC[NeffP2SGC-1])-log(k2SGC[0]))*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP2SGC-2+25+(int)(log(k2SGC[0])/(log(k2SGC[NeffP2SGC-1])-log(k2SGC[0]))*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP4SGC-2+25+(int)(log(k4SGC[0])/(log(k4SGC[NeffP4SGC-1])-log(k4SGC[0]))*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP4SGC-2+25+(int)(log(k4SGC[0])/(log(k4SGC[NeffP4SGC-1])-log(k4SGC[0]))*(NeffP4SGC-1));}

}

if( strcmp(spacing_dataSGC,"log10") == 0  ){

if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP0SGC-2+25+(int)(log10(k0SGC[0])/(log10(k0SGC[NeffP0SGC-1])-log10(k0SGC[0]))*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP0SGC-2+25+(int)(log10(k0SGC[0])/(log10(k0SGC[NeffP0SGC-1])-log10(k0SGC[0]))*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP2SGC-2+25+(int)(log10(k2SGC[0])/(log10(k2SGC[NeffP2SGC-1])-log10(k2SGC[0]))*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){NeffmaxSGC=NeffP2SGC-2+25+(int)(log10(k2SGC[0])/(log10(k2SGC[NeffP2SGC-1])-log10(k2SGC[0]))*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){NeffmaxSGC=NeffP4SGC-2+25+(int)(log10(k4SGC[0])/(log10(k4SGC[NeffP4SGC-1])-log10(k4SGC[0]))*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){NeffmaxSGC=NeffP4SGC-2+25+(int)(log10(k4SGC[0])/(log10(k4SGC[NeffP4SGC-1])-log10(k4SGC[0]))*(NeffP4SGC-1));}

}

if( strcmp(spacing_dataSGC,"irregular") == 0  ){
NeffmaxSGC=Ntheo;
factor_sampling_mask=1;
}


ktheoSGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
ksmSGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
}

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
if(strcmp(path_to_mask2, "none") == 0 ){
Ptheo0SGC = (double*) calloc( NeffP0SGC, sizeof(double));
Psm0SGC = (double*) calloc( NeffP0SGC, sizeof(double));
ktheo0SGC = (double*) calloc( NeffP0SGC, sizeof(double));
ksm0SGC = (double*) calloc( NeffP0SGC, sizeof(double));

}
else{
Ptheo0SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
Psm0SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
}

kmin0SGC=k0SGC[0];
kmax0SGC=k0SGC[NeffP0SGC-1];
}

if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
if(strcmp(path_to_mask2, "none") == 0 ){
Ptheo2SGC = (double*) calloc( NeffP2SGC, sizeof(double));
Psm2SGC = (double*) calloc( NeffP2SGC, sizeof(double));
ktheo2SGC = (double*) calloc( NeffP2SGC, sizeof(double));
ksm2SGC = (double*) calloc( NeffP2SGC, sizeof(double));
}
else{
Ptheo2SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
Psm2SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
}
kmin2SGC=k2SGC[0];
kmax2SGC=k2SGC[NeffP2SGC-1];

}

if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
if(strcmp(path_to_mask2, "none") == 0 ){
Ptheo4SGC = (double*) calloc( NeffP4SGC, sizeof(double));
Psm4SGC = (double*) calloc( NeffP4SGC, sizeof(double));
ktheo4SGC = (double*) calloc( NeffP4SGC, sizeof(double));
ksm4SGC = (double*) calloc( NeffP4SGC, sizeof(double));
}
else{
Ptheo4SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
Psm4SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
}
kmin4SGC=k4SGC[0];
kmax4SGC=k4SGC[NeffP4SGC-1];
}

/*
if(modeP0+modeP2+modeP4==1){parameters1 =  (double*) calloc( Npolynomial+1+2, sizeof(double));}
if(modeP0+modeP2+modeP4>1){

      if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 ){   parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+5, sizeof(double));}
      else{

               if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0){parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+4, sizeof(double));}
               else{

                      if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0){parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+4, sizeof(double));}
                      else{
                             parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+3, sizeof(double));
                          }
                   }
         }
}
*/
//if(modeP0+modeP2+modeP4==1){Nalphas_for_param1=1;}
//else{Nalphas_for_param1=2;}

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){

if(modeP0+modeP2+modeP4==1){Nalphas_for_param1=1;}
else{Nalphas_for_param1=2;}

Nsigmas_for_param1=modeP0+modeP2+modeP4;
parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+Nalphas_for_param1+Nsigmas_for_param1, sizeof(double));
offset=Nalphas+Nsigmas_tot;
}

if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
Nalphas_for_param1=2;
Nsigmas_for_param1=2;
parameters1 =  (double*) calloc( 1+(modeP0+modeP2+modeP4)*(Npolynomial)+Nalphas_for_param1+Nsigmas_for_param1+1, sizeof(double));
offset=Nalphas+Nsigmas_tot+1;
}

if(strcmp(Sigma_def_type, "effective") == 0)
{

if(modeP0+modeP2+modeP4==1){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];offset_ini=2;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];parameters1[4]=parameters2[2];offset_ini=5;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];parameters1[4]=parameters2[4];offset_ini=5;}




}
else
{


if(modeP0+modeP2+modeP4==1)
{

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[2],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[2],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[2],4./14.);}
}
else
{
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[1]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[1]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[1]/(1.+ffactor),4./14.);}
}

     parameters1[0]=parameters2[0];
if(modeP0==1){parameters1[1]=Sigmanl0;}
if(modeP2==1){parameters1[1]=Sigmanl2;}
if(modeP4==1){parameters1[1]=Sigmanl4;}
offset_ini=2;
}
else{
     parameters1[0]=parameters2[0];
     parameters1[1]=parameters2[1];
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[3],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[3],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[3],4./14.);}
}else
{
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[2]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[2]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[2]/(1.+ffactor),4./14.);}
}


if(modeP0==0){parameters1[2]=Sigmanl2;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP2==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP4==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;offset_ini=4;}
if(modeP0+modeP2+modeP4==3){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;parameters1[4]=Sigmanl4;offset_ini=5;}
}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){

parameters1[2]=parameters2[2];
if(strcmp(Sigma_independent, "yes") == 0)
{
parameters1[3]=parameters2[3];
parameters1[4]=parameters2[4];
}
else
{
parameters1[3]=parameters2[2]/(1+ffactor);
parameters1[4]=parameters2[3];
}
offset_ini=5;


}

}

}

//}

//if(strcmp(type_BAO_fit, "mcmc") == 0)
//{
//}

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){for(i=offset_ini;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){parameters1[i]=parameters2[i+(Npolynomial+1)*(modeP0+modeP2+modeP4)-offset_ini+offset];}}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){for(i=offset_ini;i<1+(Npolynomial)*(modeP0+modeP2+modeP4)+offset_ini;i++){parameters1[i]=parameters2[i+1+(Npolynomial)*(modeP0+modeP2+modeP4)-offset_ini+offset];}}

//for(i=0;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){printf("SGC %d %lf\n",i,parameters1[i]);}

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){

do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ktheoSGC,ktheo0SGC,ktheo2SGC,ktheo4SGC, Ptheo0SGC, Ptheo2SGC, Ptheo4SGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC, path_to_mask2,k0SGC,k2SGC,k4SGC, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],kmin0SGC, kmax0SGC, kmin2SGC,kmax2SGC,kmin4SGC,kmax4SGC, 1,spacing_dataSGC,spacing_theory,Sigma_smooth);//full power spectrum

do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ksmSGC,ksm0SGC,ksm2SGC,ksm4SGC, Psm0SGC, Psm2SGC, Psm4SGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC, path_to_mask2,k0SGC,k2SGC,k4SGC, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],kmin0SGC, kmax0SGC,kmin2SGC, kmax2SGC,kmin4SGC,kmax4SGC, 0,spacing_dataSGC,spacing_theory,Sigma_smooth);//broadband only
}

if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){

do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ktheoSGC,ktheo0SGC,ktheo2SGC,ktheo4SGC, Ptheo0SGC, Ptheo2SGC, Ptheo4SGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,Sigma_smooth, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC, path_to_mask2,k0SGC,k2SGC,k4SGC, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],kmin0SGC, kmax0SGC, kmin2SGC,kmax2SGC,kmin4SGC,kmax4SGC, 1,spacing_dataSGC,spacing_theory);//full power spectrum

do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, ksmSGC,ksm0SGC,ksm2SGC,ksm4SGC, Psm0SGC, Psm2SGC, Psm4SGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,Sigma_smooth, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC, path_to_mask2,k0SGC,k2SGC,k4SGC, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],kmin0SGC, kmax0SGC,kmin2SGC, kmax2SGC,kmin4SGC,kmax4SGC, 0,spacing_dataSGC,spacing_theory);//broadband only

}

}

if(Nchunks==2)
{
combine=0;
//   Decide if combine or not here.if combine==0 then combine, otherwise no
if(NeffP0!=NeffP0SGC || NeffP2!=NeffP2SGC || NeffP4!=NeffP4SGC)
{ 
combine=1;
}
else{

  for(i=0;i<NeffP0;i++){if(k0[i]!=k0SGC[i]){combine++;}}
}

}
if(combine==0)//decide if they are going to be combined according to Area~1/err^2 for k=0.1 bin
{
//in this point NeffP0=NeffP0SGC, NeffP2=NeffP2SGC, NeffP4=NeffP4SGC
   ik01_P0=-1;
   for(i=0;i<NeffP0;i++)
   {
     if(i==0){difference_min=fabs(k0[i]-0.1);ik01_P0=i;}
     else
     {
        difference=fabs(k0[i]-0.1);
        if(difference<difference_min){difference_min=difference;ik01_P0=i;}

     }
  }
   ik01_P2=-1;
   for(i=0;i<NeffP2;i++)
   {
     if(i==0){difference_min=fabs(k2[i]-0.1);ik01_P2=i;}
     else
     {
        difference=fabs(k2[i]-0.1);
        if(difference<difference_min){difference_min=difference;ik01_P2=i;}

     }
  }
   ik01_P4=-1;
   for(i=0;i<NeffP4;i++)
   {
     if(i==0){difference_min=fabs(k4[i]-0.1);ik01_P4=i;}
     else
     {
        difference=fabs(k4[i]-0.1);
        if(difference<difference_min){difference_min=difference;ik01_P4=i;}

     }
  }



}


sprintf(nameP0_1,"%s/Monopole1_bao_%s.txt",path_output,identifier);
sprintf(nameP0_2,"%s/Monopole2_bao_%s.txt",path_output,identifier);
sprintf(nameP0_12,"%s/Monopole12_bao_%s.txt",path_output,identifier);

sprintf(nameP2_1,"%s/Quadrupole1_bao_%s.txt",path_output,identifier);
sprintf(nameP2_2,"%s/Quadrupole2_bao_%s.txt",path_output,identifier);
sprintf(nameP2_12,"%s/Quadrupole12_bao_%s.txt",path_output,identifier);

sprintf(nameP4_1,"%s/Hexadecapole1_bao_%s.txt",path_output,identifier);
sprintf(nameP4_2,"%s/Hexadecapole2_bao_%s.txt",path_output,identifier);
sprintf(nameP4_12,"%s/Hexadecapole12_bao_%s.txt",path_output,identifier);



for(l=0;l<=4;l=l+2)
{
open=0;
if(modeP0==1 && l==0){string_copy(nameP0_1,name1);string_copy(nameP0_2,name2);string_copy(nameP0_12,name12);open=1;}
if(modeP2==1 && l==2){string_copy(nameP2_1,name1);string_copy(nameP2_2,name2);string_copy(nameP2_12,name12);open=1;}
if(modeP4==1 && l==4){string_copy(nameP4_1,name1);string_copy(nameP4_2,name2);string_copy(nameP4_12,name12);open=1;}

if(open==1){

f1=fopen(name1,"w");
if(Nchunks==2)
{
f2=fopen(name2,"w");
if(combine==0){f12=fopen(name12,"w");}
}

fprintf(f1,"#k \t Pdata \t err \t Pmodel \t Pnowiggle\n");

if(Nchunks==2)
{
fprintf(f2,"#k \t Pdata \t err \t Pmodel \t Pnowiggle\n");
if(combine==0){fprintf(f12,"#k \t Pdata \t err \t Pmodel \t Pnowiggle\n");}
}

if(l==0){NeffP=NeffP0;}
if(l==2){NeffP=NeffP2;}
if(l==4){NeffP=NeffP4;}

for(i=0;i<NeffP;i++)
{

if(l==0)
{
if(strcmp(path_to_mask1, "none") == 0 ){
ptheo=Ptheo0[i];
pnw=Psm0[i];
}
else{

Ninterpol=determine_N_singlearray(ktheo,k0[i],Neffmax*factor_sampling_mask,spacing_dataNGC);
Ninterpolsm=determine_N_singlearray(ksm,k0[i],Neffmax*factor_sampling_mask,spacing_dataNGC);

if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0[i]<=0 || Ninterpolsm>=Neffmax*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheo=0;
pnw=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_singlearray(ksm,k0[i],Ninterpol,spacing_dataNGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);

w0sm=determine_w0_2ndorder_singlearray(ksm,k0[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_2ndorder_singlearray(ksm,k0[i],Ninterpol,spacing_dataNGC);
w2sm=determine_w2_2ndorder_singlearray(ksm,k0[i],Ninterpol,spacing_dataNGC);
}
ptheo=P_interpol_fast(k0[i],Ptheo0,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
pnw=P_interpol_fast(k0[i],Psm0,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);

}


//ptheo=P_interpol(k0[i],ktheo,Ptheo0,Neffmax*factor_sampling_mask);
//pnw=P_interpol(k0[i],ksm,Psm0,Neffmax*factor_sampling_mask);


}


fprintf(f1,"%e %e %e %e %e\n",k0[i],P0[i],errP0[i],ptheo,pnw);
}

if(l==2)
{
if(strcmp(path_to_mask1, "none") == 0 ){
ptheo=Ptheo2[i];
pnw=Psm2[i];
}
else{

Ninterpol=determine_N_singlearray(ktheo,k2[i],Neffmax*factor_sampling_mask,spacing_dataNGC);
Ninterpolsm=determine_N_singlearray(ksm,k2[i],Neffmax*factor_sampling_mask,spacing_dataNGC);

if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2[i]<=0 || Ninterpolsm>=Neffmax*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheo=0;
pnw=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_singlearray(ksm,k2[i],Ninterpol,spacing_dataNGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);

w0sm=determine_w0_2ndorder_singlearray(ksm,k2[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_2ndorder_singlearray(ksm,k2[i],Ninterpol,spacing_dataNGC);
w2sm=determine_w2_2ndorder_singlearray(ksm,k2[i],Ninterpol,spacing_dataNGC);
}
ptheo=P_interpol_fast(k2[i],Ptheo2,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
pnw=P_interpol_fast(k2[i],Psm2,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);

}


//ptheo=P_interpol(k2[i],ktheo,Ptheo2,Neffmax*factor_sampling_mask);
//pnw=P_interpol(k2[i],ksm,Psm2,Neffmax*factor_sampling_mask);
}
fprintf(f1,"%e %e %e %e %e\n",k2[i],P2[i],errP2[i],ptheo,pnw);
}

if(l==4)
{
if(strcmp(path_to_mask1, "none") == 0 ){
ptheo=Ptheo4[i];
pnw=Psm4[i];
}
else{

Ninterpol=determine_N_singlearray(ktheo,k4[i],Neffmax*factor_sampling_mask,spacing_dataNGC);
Ninterpolsm=determine_N_singlearray(ksm,k4[i],Neffmax*factor_sampling_mask,spacing_dataNGC);

if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4[i]<=0 || Ninterpolsm>=Neffmax*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheo=0;
pnw=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_singlearray(ksm,k4[i],Ninterpol,spacing_dataNGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);

w0sm=determine_w0_2ndorder_singlearray(ksm,k4[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_2ndorder_singlearray(ksm,k4[i],Ninterpol,spacing_dataNGC);
w2sm=determine_w2_2ndorder_singlearray(ksm,k4[i],Ninterpol,spacing_dataNGC);
}
ptheo=P_interpol_fast(k4[i],Ptheo4,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
pnw=P_interpol_fast(k4[i],Psm4,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);

}


//ptheo=P_interpol(k4[i],ktheo,Ptheo4,Neffmax*factor_sampling_mask);
//pnw=P_interpol(k4[i],ksm,Psm4,Neffmax*factor_sampling_mask);
}
fprintf(f1,"%e %e %e %e %e\n",k4[i],P4[i],errP4[i],ptheo,pnw);
}



if(Nchunks==2)
{
//This part is written after, separated from this loop in case Neff is not NeffSGC

if(combine==0)
{


if(l==0)
{

if(strcmp(path_to_mask2, "none") == 0 ){
ptheoSGC=Ptheo0SGC[i];
pnwSGC=Psm0SGC[i];
}
else{

Ninterpol=determine_N_singlearray(ktheoSGC,k0SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
Ninterpolsm=determine_N_singlearray(ksmSGC,k0SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);

if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k0SGC[i]<=0 || Ninterpolsm>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheoSGC=0;
pnwSGC=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);

w0sm=determine_w0_2ndorder_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_2ndorder_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w2sm=determine_w2_2ndorder_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k0SGC[i],Ptheo0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
pnwSGC=P_interpol_fast(k0SGC[i],Psm0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);

}


//ptheoSGC=P_interpol(k0SGC[i],ktheoSGC,Ptheo0SGC,NeffmaxSGC*factor_sampling_mask);
//pnwSGC=P_interpol(k0SGC[i],ksmSGC,Psm0SGC,NeffmaxSGC*factor_sampling_mask);
}

}

if(l==2)
{

if(strcmp(path_to_mask2, "none") == 0 ){
ptheoSGC=Ptheo2SGC[i];
pnwSGC=Psm2SGC[i];
}else{

Ninterpol=determine_N_singlearray(ktheoSGC,k2SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
Ninterpolsm=determine_N_singlearray(ksmSGC,k2SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);

if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k2SGC[i]<=0 || Ninterpolsm>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheoSGC=0;
pnwSGC=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);

w0sm=determine_w0_2ndorder_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_2ndorder_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w2sm=determine_w2_2ndorder_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k2SGC[i],Ptheo2SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
pnwSGC=P_interpol_fast(k2SGC[i],Psm2SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);

}



//ptheoSGC=P_interpol(k2SGC[i],ktheoSGC,Ptheo2SGC,NeffmaxSGC*factor_sampling_mask);
//pnwSGC=P_interpol(k2SGC[i],ksmSGC,Psm2SGC,NeffmaxSGC*factor_sampling_mask);
}

}

if(l==4)
{

if(strcmp(path_to_mask2, "none") == 0 ){
ptheoSGC=Ptheo4SGC[i];
pnwSGC=Psm4SGC[i];
}else{

Ninterpol=determine_N_singlearray(ktheo,k4SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
Ninterpolsm=determine_N_singlearray(ksm,k4SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);

if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k4SGC[i]<=0 || Ninterpolsm>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheoSGC=0;
pnwSGC=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);

w0sm=determine_w0_2ndorder_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_2ndorder_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w2sm=determine_w2_2ndorder_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k4SGC[i],Ptheo4SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
pnwSGC=P_interpol_fast(k4SGC[i],Psm4SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);

}


//ptheoSGC=P_interpol(k4SGC[i],ktheoSGC,Ptheo4SGC,NeffmaxSGC*factor_sampling_mask);
//pnwSGC=P_interpol(k4SGC[i],ksmSGC,Psm4SGC,NeffmaxSGC*factor_sampling_mask);
}

}

if(l==0)
{

P0tot=(P0[i]*pow(errP0[ik01_P0],-2)+P0SGC[i]*pow(errP0SGC[ik01_P0],-2))/( pow(errP0[ik01_P0],-2)+pow(errP0SGC[ik01_P0],-2) );

ptheotot=(ptheo*pow(errP0[ik01_P0],-2)+ptheoSGC*pow(errP0SGC[ik01_P0],-2))/( pow(errP0[ik01_P0],-2)+pow(errP0SGC[ik01_P0],-2) );

pnwtot=(pnw*pow(errP0[ik01_P0],-2)+pnwSGC*pow(errP0SGC[ik01_P0],-2))/( pow(errP0[ik01_P0],-2)+pow(errP0SGC[ik01_P0],-2) );

errP0tot=errP0[i]*sqrt( pow(errP0[ik01_P0],-2)/(pow(errP0[ik01_P0],-2)+pow(errP0SGC[ik01_P0],-2)) );//errNS=errNGC*sqrt(ANGC/(ANGC+ASGC))

fprintf(f12,"%e %e %e %e %e\n",k0[i],P0tot,errP0tot,ptheotot,pnwtot);
}

if(l==2)
{
P0tot=(P2[i]*pow(errP2[ik01_P2],-2)+P2SGC[i]*pow(errP2SGC[ik01_P2],-2))/( pow(errP2[ik01_P2],-2)+pow(errP2SGC[ik01_P2],-2) );

ptheotot=(ptheo*pow(errP2[ik01_P2],-2)+ptheoSGC*pow(errP2SGC[ik01_P2],-2))/( pow(errP2[ik01_P2],-2)+pow(errP2SGC[ik01_P2],-2) );

pnwtot=(pnw*pow(errP2[ik01_P2],-2)+pnwSGC*pow(errP2SGC[ik01_P2],-2))/( pow(errP2[ik01_P2],-2)+pow(errP2SGC[ik01_P2],-2) );

errP0tot=errP2[i]*sqrt( pow(errP2[ik01_P2],-2)/(pow(errP2[ik01_P2],-2)+pow(errP2SGC[ik01_P2],-2)) );//errNS=errNGC*sqrt(ANGC/(ANGC+ASGC))

fprintf(f12,"%e %e %e %e %e\n",k2[i],P0tot,errP0tot,ptheotot,pnwtot);

}

if(l==4)
{
P0tot=(P4[i]*pow(errP4[ik01_P4],-2)+P4SGC[i]*pow(errP4SGC[ik01_P4],-2))/( pow(errP4[ik01_P4],-2)+pow(errP4SGC[ik01_P4],-2) );

ptheotot=(ptheo*pow(errP4[ik01_P4],-2)+ptheoSGC*pow(errP4SGC[ik01_P4],-2))/( pow(errP4[ik01_P4],-2)+pow(errP4SGC[ik01_P4],-2) );

pnwtot=(pnw*pow(errP4[ik01_P4],-2)+pnwSGC*pow(errP4SGC[ik01_P4],-2))/( pow(errP4[ik01_P4],-2)+pow(errP4SGC[ik01_P4],-2) );

errP0tot=errP4[i]*sqrt( pow(errP4[ik01_P4],-2)/(pow(errP4[ik01_P4],-2)+pow(errP4SGC[ik01_P4],-2)) );//errNS=errNGC*sqrt(ANGC/(ANGC+ASGC))

fprintf(f12,"%e %e %e %e %e\n",k4[i],P0tot,errP0tot,ptheotot,pnwtot);

}

}

}//Chunks=2

}//for-NeffP

if(Nchunks==2)
{

if(l==0){NeffP=NeffP0SGC;}
if(l==2){NeffP=NeffP2SGC;}
if(l==4){NeffP=NeffP4SGC;}

for(i=0;i<NeffP;i++)
{

if(l==0)
{

if(strcmp(path_to_mask2, "none") == 0 ){
ptheoSGC=Ptheo0SGC[i];
pnwSGC=Psm0SGC[i];
}
else{

Ninterpol=determine_N_singlearray(ktheoSGC,k0SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
Ninterpolsm=determine_N_singlearray(ksmSGC,k0SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);

if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k0SGC[i]<=0 || Ninterpolsm>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheoSGC=0;
pnwSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);

w0sm=determine_w0_2ndorder_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_2ndorder_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataNGC);
w2sm=determine_w2_2ndorder_singlearray(ksmSGC,k0SGC[i],Ninterpol,spacing_dataNGC);
}
ptheoSGC=P_interpol_fast(k0SGC[i],Ptheo0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
pnwSGC=P_interpol_fast(k0SGC[i],Psm0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);
}



//ptheoSGC=P_interpol(k0SGC[i],ktheoSGC,Ptheo0SGC,NeffmaxSGC*factor_sampling_mask);
//pnwSGC=P_interpol(k0SGC[i],ksmSGC,Psm0SGC,NeffmaxSGC*factor_sampling_mask);


}

fprintf(f2,"%e %e %e %e %e\n",k0SGC[i],P0SGC[i],errP0SGC[i],ptheoSGC,pnwSGC);
}

if(l==2)
{

if(strcmp(path_to_mask2, "none") == 0 ){
ptheoSGC=Ptheo2SGC[i];
pnwSGC=Psm2SGC[i];
}
else{

Ninterpol=determine_N_singlearray(ktheoSGC,k2SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
Ninterpolsm=determine_N_singlearray(ksmSGC,k2SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);

if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k2SGC[i]<=0 || Ninterpolsm>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheoSGC=0;
pnwSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);

w0sm=determine_w0_2ndorder_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_2ndorder_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataNGC);
w2sm=determine_w2_2ndorder_singlearray(ksmSGC,k2SGC[i],Ninterpol,spacing_dataNGC);
}
ptheoSGC=P_interpol_fast(k2SGC[i],Ptheo2SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
pnwSGC=P_interpol_fast(k2SGC[i],Psm2SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);
}

//ptheoSGC=P_interpol(k2SGC[i],ktheoSGC,Ptheo2SGC,NeffmaxSGC*factor_sampling_mask);
//pnwSGC=P_interpol(k2SGC[i],ksmSGC,Psm2SGC,NeffmaxSGC*factor_sampling_mask);
}

fprintf(f2,"%e %e %e %e %e\n",k2SGC[i],P2SGC[i],errP2SGC[i],ptheoSGC,pnwSGC);
}

if(l==4)
{
if(strcmp(path_to_mask2, "none") == 0 ){
ptheoSGC=Ptheo4SGC[i];
pnwSGC=Psm4SGC[i];
}
else{

Ninterpol=determine_N_singlearray(ktheoSGC,k4SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
Ninterpolsm=determine_N_singlearray(ksmSGC,k4SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);

if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k4SGC[i]<=0 || Ninterpolsm>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpolsm<0)
{
ptheoSGC=0;
pnwSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w1sm=determine_w1_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);

w0sm=determine_w0_2ndorder_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataNGC);
w1sm=determine_w1_2ndorder_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataNGC);
w2sm=determine_w2_2ndorder_singlearray(ksmSGC,k4SGC[i],Ninterpol,spacing_dataNGC);
}
ptheoSGC=P_interpol_fast(k4SGC[i],Ptheo4SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
pnwSGC=P_interpol_fast(k4SGC[i],Psm4SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0sm,w1sm,w2sm);
}


//ptheoSGC=P_interpol(k4SGC[i],ktheoSGC,Ptheo4SGC,NeffmaxSGC*factor_sampling_mask);
//pnwSGC=P_interpol(k4SGC[i],ksmSGC,Psm4SGC,NeffmaxSGC*factor_sampling_mask);
}

fprintf(f2,"%e %e %e %e %e\n",k4SGC[i],P4SGC[i],errP4SGC[i],ptheoSGC,pnwSGC);
}


}//for NeffP

}//if Nchunks=2

fprintf(f1,"\n#");
fprintf(f1,"B\t");
if(Nchunks==2)
{
fprintf(f2,"\n#");
if(combine==0){fprintf(f12,"\n#B\t");
}
}

for(j=1;j<=Npolynomial;j++){fprintf(f1,"A%d\t",j);}
if(Nchunks==2)
{
if(combine==0){for(j=1;j<=Npolynomial;j++){fprintf(f12,"A%d\t",j);}
}
}


if(Nchunks==2)
{
fprintf(f2,"BSGC\t");
if(combine==0){fprintf(f12,"BSGC\t");
}
}


if(Nchunks==2)
{
for(j=1;j<=Npolynomial;j++){fprintf(f2,"ASGC%d\t",j);}
if(combine==0){for(j=1;j<=Npolynomial;j++){fprintf(f12,"ASGC%d\t",j);}
}
}


if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f1,"Sigmanl_eff \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f1,"Sigmanl_eff \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f1,"Sigmanl_eff \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){fprintf(f1,"Sigmanl_eff \t alpha_para \t alpha_perp \t chi2\t");}
}
if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f1,"Sigmanl0_eff \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f1,"Sigmanl2_eff \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f1,"Sigmanl4_eff \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 ){fprintf(f1,"Sigmanl0_eff \t Sigmanl2_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P04") == 0 ){fprintf(f1,"Sigmanl0_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P24") == 0 ){fprintf(f1,"Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P024") == 0 ){fprintf(f1,"Sigmanl0_eff \t Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f1,"Sigmanl_para \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f1,"Sigmanl_para \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f1,"Sigmanl_para \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f1,"Sigmanl_para \t alpha_para \t alpha_perp \t chi2\t");}

if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f1,"Sigmanl_para \t beta \t alpha_para \t alpha_perp \t chi2\t");}

}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f1,"Sigmanl_para \t Sigma_perp \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f1,"Sigmanl_para \t Sigma_perp \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f1,"Sigmanl_para \t Sigma_perp \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 ||  strcmp(fit_BAO, "P04") == 0 ||  strcmp(fit_BAO, "P24") == 0 ||  strcmp(fit_BAO, "P024") == 0){

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f1,"Sigmanl_para \t Sigmanl_perp \t alpha_para \t alpha_perp \t chi2\t");}

if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f1,"Sigmanl_para \t Sigmanl_perp \t beta \t alpha_para \t alpha_perp \t chi2\t");}

}

}


//if(strcmp(fit_BAO, "P0") == 0){fprintf(f1,"Sigmanl \t alpha0 \t chi2\t");}
//if(strcmp(fit_BAO, "P2") == 0){fprintf(f1,"Sigmanl \t alpha2 \t chi2\t");}
//if(strcmp(fit_BAO, "P4") == 0){fprintf(f1,"Sigmanl \t alpha4 \t chi2\t");}
//if(strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){fprintf(f1,"Sigmanl \t alpha_para \t alpha_perp \t chi2\t");}

fprintf(f1," Npoints-dof\n");
if(Nchunks==2)
{
//if(strcmp(fit_BAO, "P0") == 0){fprintf(f2,"Sigmanl \t alpha0 \t chi2\t");}
//if(strcmp(fit_BAO, "P2") == 0){fprintf(f2,"Sigmanl \t alpha2 \t chi2\t");}
//if(strcmp(fit_BAO, "P4") == 0){fprintf(f2,"Sigmanl \t alpha4 \t chi2\t");}
//if(strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){fprintf(f2,"Sigmanl \t alpha_para \t alpha_perp \t chi2\t");}

if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f2,"Sigmanl_eff \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f2,"Sigmanl_eff \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f2,"Sigmanl_eff \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){fprintf(f2,"Sigmanl_eff \t alpha_para \t alpha_perp \t chi2\t");}
}
if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f2,"Sigmanl0_eff \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f2,"Sigmanl2_eff \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f2,"Sigmanl4_eff \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 ){fprintf(f2,"Sigmanl0_eff \t Sigmanl2_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P04") == 0 ){fprintf(f2,"Sigmanl0_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P24") == 0 ){fprintf(f2,"Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P024") == 0 ){fprintf(f2,"Sigmanl0_eff \t Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f2,"Sigmanl_para \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f2,"Sigmanl_para \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f2,"Sigmanl_para \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){

if( strcmp(type_of_analysis, "BAOISO") == 0 ){fprintf(f2,"Sigmanl_para \t alpha_para \t alpha_perp \t chi2\t");}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){fprintf(f2,"Sigmanl_para \t beta \t alpha_para \t alpha_perp \t chi2\t");}


}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f2,"Sigmanl_para \t Sigma_perp \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f2,"Sigmanl_para \t Sigma_perp \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f2,"Sigmanl_para \t Sigma_perp \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 ||  strcmp(fit_BAO, "P04") == 0 ||  strcmp(fit_BAO, "P24") == 0 ||  strcmp(fit_BAO, "P024") == 0){

if( strcmp(type_of_analysis, "BAOISO") == 0 ){fprintf(f2,"Sigmanl_para \t Sigmanl_perp \t alpha_para \t alpha_perp \t chi2\t");}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){fprintf(f2,"Sigmanl_para \t Sigmanl_perp \t beta \t alpha_para \t alpha_perp \t chi2\t");}

}
}

fprintf(f2," Npoints-dof\n");
if(combine==0){

if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f12,"Sigmanl_eff \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f12,"Sigmanl_eff \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f12,"Sigmanl_eff \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){fprintf(f12,"Sigmanl_eff \t alpha_para \t alpha_perp \t chi2\t");}
}
if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f12,"Sigmanl0_eff \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f12,"Sigmanl2_eff \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f12,"Sigmanl4_eff \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 ){fprintf(f12,"Sigmanl0_eff \t Sigmanl2_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P04") == 0 ){fprintf(f12,"Sigmanl0_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P24") == 0 ){fprintf(f12,"Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
if( strcmp(fit_BAO, "P024") == 0 ){fprintf(f12,"Sigmanl0_eff \t Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\t");}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f12,"Sigmanl_para \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f12,"Sigmanl_para \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f12,"Sigmanl_para \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){

if( strcmp(type_of_analysis, "BAOISO") == 0 ){fprintf(f12,"Sigmanl_para \t alpha_para \t alpha_perp \t chi2\t");}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){fprintf(f12,"Sigmanl_para \t alpha_para \t beta \t alpha_perp \t chi2\t");}



}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f12,"Sigmanl_para \t Sigma_perp \t alpha0 \t chi2\t");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f12,"Sigmanl_para \t Sigma_perp \t alpha2 \t chi2\t");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f12,"Sigmanl_para \t Sigma_perp \t alpha4 \t chi2\t");}
if( strcmp(fit_BAO, "P02") == 0 ||  strcmp(fit_BAO, "P04") == 0 ||  strcmp(fit_BAO, "P24") == 0 ||  strcmp(fit_BAO, "P024") == 0){

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f12,"Sigmanl_para \t Sigmanl_perp \t alpha_para \t alpha_perp \t chi2\t");}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ||  strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f12,"Sigmanl_para \t Sigmanl_perp \t beta \t alpha_para \t alpha_perp \t chi2\t");}


}
}


//if(strcmp(fit_BAO, "P0") == 0){fprintf(f12,"Sigmanl \t alpha0 \t chi2\t");}
//if(strcmp(fit_BAO, "P2") == 0){fprintf(f12,"Sigmanl \t alpha2 \t chi2\t");}
//if(strcmp(fit_BAO, "P4") == 0){fprintf(f12,"Sigmanl \t alpha4 \t chi2\t");}
//if(strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){fprintf(f12,"Sigmanl \t alpha_para \t alpha_perp \t chi2\t");}

fprintf(f12," Npoints-dof\n");
}
}
///write best-fit parameters


for(j=0;j<=Npolynomial;j++){
if(j==0){fprintf(f1,"#");}

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f1,"%e\t",parameters2[j+Nalphas+Nsigmas_tot]);}//either P0, P2 or P4
else{fprintf(f1,"%e\t",parameters2[ j+(Npolynomial+1)*(l/2) + Nalphas+Nsigmas_tot]);}
}

if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f1,"%e\t",parameters2[j+Nalphas+Nsigmas_tot]);}//either P0, P2 or P4
else{

if(l==0){fprintf(f1,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot]);}
else{

if(j==0){fprintf(f1,"%e\t",parameters2[ j+(1+Npolynomial*(0/2)) + Nalphas+Nsigmas_tot]);}
if(j>0){fprintf(f1,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot]);}

}

}

}

}

if(Nchunks==2){


for(j=0;j<=Npolynomial;j++){
if(j==0){fprintf(f2,"#");}
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f2,"%e\t",parameters2[j+Nalphas+Nsigmas_tot+ (Npolynomial+1)]);}
else{fprintf(f2,"%e\t",parameters2[ j+(Npolynomial+1)*(l/2) + Nalphas+Nsigmas_tot + (Npolynomial+1)*(modeP0+modeP2+modeP4) ]);}
}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f2,"%e\t",parameters2[j+Nalphas+Nsigmas_tot + (Npolynomial+1)]);}
else{

//fprintf(f2,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot + 1+(Npolynomial)*(modeP0+modeP2+modeP4) ]);

if(l==0){fprintf(f2,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot+ 1+(Npolynomial)*(modeP0+modeP2+modeP4)]);}
else{

if(j==0){fprintf(f2,"%e\t",parameters2[ j+(1+Npolynomial*(0/2)) + Nalphas+Nsigmas_tot+ 1+(Npolynomial)*(modeP0+modeP2+modeP4)]);}
if(j>0){fprintf(f2,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot+ 1+(Npolynomial)*(modeP0+modeP2+modeP4)]);}

}

}

}

}


//if(modeP0+modeP2+modeP4==1){fprintf(f2,"%e\t %e \t %e \t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
//else{fprintf(f2,"%e\t %e\t %e \t %e \t %d-%d\n",parameters2[2],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}
if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "no") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f2,"%e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[2],parameters2[0], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f2,"%e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[1], parameters2[2], parameters2[0], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[1], parameters2[2],parameters2[3], parameters2[0], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "no") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[2],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2],parameters2[3],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}

}
if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[0],parameters2[1],chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3],parameters2[4], parameters2[0],parameters2[1],chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}
}

if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f2,"%e\t %e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[4], parameters2[0], parameters2[1], chi2_min,Npoints,Ndof);}
}


if(combine==0){

for(j=0;j<=Npolynomial;j++){
if(j==0){fprintf(f12,"#");}
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f12,"%e\t",parameters2[j+Nalphas+Nsigmas_tot]);}//NGC
else{fprintf(f12,"%e\t",parameters2[ j+(Npolynomial+1)*(l/2) + Nalphas+Nsigmas_tot ]);}//NGC
}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f12,"%e\t",parameters2[j+Nalphas+Nsigmas_tot]);}//NGC
else{

//fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot ]);

if(l==0){fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot]);}
else{

if(j==0){fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(0/2)) + Nalphas+Nsigmas_tot]);}
if(j>0){fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot]);}

}

}//NGC

}


}

for(j=0;j<=Npolynomial;j++){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f12,"%e\t",parameters2[Nalphas+Nsigmas_tot+Npolynomial+1+j]);}//SGC
else{fprintf(f12,"%e\t",parameters2[ j+(Npolynomial+1)*(l/2) + Nalphas+Nsigmas_tot + (Npolynomial+1)*(modeP0+modeP2+modeP4)]);}
}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
if(modeP0+modeP2+modeP4==1){fprintf(f12,"%e\t",parameters2[Nalphas+Nsigmas_tot+Npolynomial+1+j]);}//SGC
else{

//fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot + (Npolynomial*(modeP0+modeP2+modeP4)+1)]);

if(l==0){fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot+ 1+(Npolynomial)*(modeP0+modeP2+modeP4)]);}
else{

if(j==0){fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(0/2)) + Nalphas+Nsigmas_tot+ 1+(Npolynomial)*(modeP0+modeP2+modeP4)]);}
if(j>0){fprintf(f12,"%e\t",parameters2[ j+(1+Npolynomial*(l/2)) + Nalphas+Nsigmas_tot+ 1+(Npolynomial)*(modeP0+modeP2+modeP4)]);}

}


}

}

}

//if(modeP0+modeP2+modeP4==1){fprintf(f12,"%e\t %e \t %e \t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
//else{fprintf(f12,"%e\t %e\t %e \t %e \t %d-%d\n",parameters2[2],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}
if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "no") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f12,"%e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f12,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[3],parameters2[0], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f12,"%e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f12,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[1], parameters2[2], parameters2[0], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f12,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[2],parameters2[3], parameters2[0], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "no") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f12,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[2],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 ){fprintf(f12,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2],parameters2[3],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}

}
if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f12,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[0],parameters2[1],chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 ||  strcmp(type_of_analysis, "FSBAOANISO") == 0){fprintf(f12,"%e\t %e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3],parameters2[4], parameters2[0],parameters2[1],chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f12,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}
}

if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f12,"%e\t %e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[4], parameters2[0], parameters2[1], chi2_min,Npoints,Ndof);}
}


}

}

//if(modeP0+modeP2+modeP4==1){fprintf(f1,"%e\t %e \t %e \t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
//else{fprintf(f1,"%e\t %e\t %e \t %e \t %d-%d\n",parameters2[2],parameters2[0], parameters2[1], chi2_min,Npoints,Ndof);}
if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "no") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){fprintf(f1,"%e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[2],parameters2[0], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %d-%d\n",parameters2[1],parameters2[0], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[1], parameters2[2], parameters2[0], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[1], parameters2[2],parameters2[3],parameters2[0], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "no") == 0 ){

if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %d-%d\n",parameters2[2],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2],parameters2[3],parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}

}
if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[0],parameters2[1],chi2_min,Npoints,Ndof);}
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3],parameters2[4], parameters2[0],parameters2[1],chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[0],parameters2[1], chi2_min,Npoints,Ndof);}

}

if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
if( strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 ){fprintf(f1,"%e\t %e\t %e\t %e\t %e\t %e\t %d-%d\n",parameters2[2], parameters2[3], parameters2[4], parameters2[0], parameters2[1], chi2_min,Npoints,Ndof);}

}




//close files
fclose(f1);
if(Nchunks==2)
{
fclose(f2);
if(combine==0){fclose(f12);}
}

}//open

}//for l=0,2,4

if(modeP0==1){
if(strcmp(path_to_mask1, "none") == 0 ){free(ktheo0);free(ksm0);}
free(Ptheo0);
free(Psm0);}

if(modeP2==1){
if(strcmp(path_to_mask1, "none") == 0 ){free(ktheo2);free(ksm2);}
free(Ptheo2);
free(Psm2);}

if(modeP4==1){
if(strcmp(path_to_mask1, "none") == 0 ){free(ktheo4);free(ksm4);}
free(Ptheo4);
free(Psm4);}

if(strcmp(path_to_mask1, "none") != 0 ){
free(ktheo);
free(ksm);
}

free(parameters1);

if(Nchunks==2)
{
if(modeP0==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo0SGC);free(ksm0SGC);}
free(Ptheo0SGC);
free(Psm0SGC);}

if(modeP2==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo2SGC);free(ksm2SGC);}
free(Ptheo2SGC);
free(Psm2SGC);}

if(modeP4==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo4SGC);free(ksm4SGC);}
free(Ptheo4SGC);
free(Psm4SGC);}

if(strcmp(path_to_mask2, "none") != 0 ){
free(ktheoSGC);
free(ksmSGC);
}
}


}


double chi2_bao(char *type_BAO_fit,char *type_of_analysis,char *fit_BAO,double *parameters2, double *k_Plin, double *Plin,int Nlin,double *k_Olin,double *Olin,int NOlin, double *pos, double *W0,double *W2,double *W4,double *W6,double *W8, int Nmask, char *path_to_mask1,char *spacing_maskNGC, double *posSGC, double *W0SGC,double *W2SGC,double *W4SGC,double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC,double *k0, double *P0,int NeffP0,double *k2, double *P2,int NeffP2,double *k4,double *P4, int NeffP4, double *k0SGC, double *P0SGC, int NeffP0SGC,double *k2SGC, double *P2SGC,int NeffP2SGC, double *k4SGC, double *P4SGC, int NeffP4SGC, double *cov, double *covSGC,  char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks, fftw_plan plan1, fftw_plan plan2, char *do_power_spectrum, char *do_bispectrum,int Nalphas,int Nsigmas_tot, int Nsigmas_free, double Sigma_smooth,int factor_sampling_mask_in,char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory)
{
double prior_chi2;
double ch2,ch2SGC;
//double ch2_diag;
double *difference;
double *parameters1;
int i,j;
double ptheo,pobs;
int Ntheo=Nlin;
double *k_theo,*k_theo0,*k_theo2,*k_theo4;
double *P_theo0,*P_theo2,*P_theo4;
int modeP0,modeP2,modeP4;
int Ncov;
int points;
int offset;
int offset_ini;
double Sigmanl0,Sigmanl2,Sigmanl4;
int Nsigmas_for_param1;
int Nalphas_for_param1;
int Neffmax;
int factor_sampling_mask;
int interpolation_order,shiftN;
double w1,w2,w0;
int Ninterpol;

interpolation_order=1;

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}

factor_sampling_mask=factor_sampling_mask_in;



ch2=0;
ch2SGC=0;
modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

Ncov=NeffP0*modeP0+NeffP2*modeP2+NeffP4*modeP4;
points=Ncov;

if(strcmp(path_to_mask1, "none") == 0 )//no mask
{
if(NeffP0*modeP0>0){k_theo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2*modeP2>0){k_theo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4*modeP4>0){k_theo4 = (double*) calloc( NeffP4, sizeof(double));}
if(NeffP0*modeP0>0){P_theo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2*modeP2>0){P_theo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4*modeP4>0){P_theo4 = (double*) calloc( NeffP4, sizeof(double));}
}
else
{

/*
if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
*/

if( strcmp(spacing_dataNGC,"linear") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
//printf("Neffmax=%d\n",Neffmax);
}

if( strcmp(spacing_dataNGC,"log") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(log(k0[0])/(log(k0[NeffP0-1])-log(k0[0]))*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(log(k0[0])/(log(k0[NeffP0-1])-log(k0[0]))*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(log(k2[0])/(log(k2[NeffP2-1])-log(k2[0]))*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(log(k2[0])/(log(k2[NeffP2-1])-log(k2[0]))*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(log(k4[0])/(log(k4[NeffP4-1])-log(k4[0]))*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(log(k4[0])/(log(k4[NeffP4-1])-log(k4[0]))*(NeffP4-1));}

}

if( strcmp(spacing_dataNGC,"log10") == 0  ){

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(log10(k0[0])/(log10(k0[NeffP0-1])-log10(k0[0]))*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(log10(k0[0])/(log10(k0[NeffP0-1])-log10(k0[0]))*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(log10(k2[0])/(log10(k2[NeffP2-1])-log10(k2[0]))*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(log10(k2[0])/(log10(k2[NeffP2-1])-log10(k2[0]))*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(log10(k4[0])/(log10(k4[NeffP4-1])-log10(k4[0]))*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(log10(k4[0])/(log10(k4[NeffP4-1])-log10(k4[0]))*(NeffP4-1));}

}

if( strcmp(spacing_dataNGC,"irregular") == 0  ){
Neffmax=Ntheo;
factor_sampling_mask=1;
}


k_theo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0==1){P_theo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2==1){P_theo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4==1){P_theo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
}

difference = (double*) calloc( Ncov, sizeof(double));


if( strcmp(type_of_analysis, "BAOISO") == 0 ){

if(modeP0+modeP2+modeP4==1){Nalphas_for_param1=1;}
else{Nalphas_for_param1=2;}

Nsigmas_for_param1=modeP0+modeP2+modeP4;

parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+Nalphas_for_param1+Nsigmas_for_param1, sizeof(double));

offset=Nalphas+Nsigmas_tot;

}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){

Nalphas_for_param1=2;
Nsigmas_for_param1=2;
parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1, sizeof(double));

offset=Nalphas+Nsigmas_tot+1;

}

if(strcmp(Sigma_def_type, "effective") == 0)
{

if(modeP0+modeP2+modeP4==1){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];offset_ini=2;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];parameters1[4]=parameters2[2];offset_ini=5;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];parameters1[4]=parameters2[4];offset_ini=5;}


}
else//para-perp
{

if(modeP0+modeP2+modeP4==1)
{

if(strcmp(Sigma_independent, "yes") == 0)
{     
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[2],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[2],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[2],4./14.);}
}
else
{
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[1]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[1]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[1]/(1.+ffactor),4./14.);}
}

     parameters1[0]=parameters2[0];    
if(modeP0==1){parameters1[1]=Sigmanl0;}
if(modeP2==1){parameters1[1]=Sigmanl2;}
if(modeP4==1){parameters1[1]=Sigmanl4;}
offset_ini=2;
}else{

     parameters1[0]=parameters2[0];    
     parameters1[1]=parameters2[1];    

if( strcmp(type_of_analysis, "BAOISO") == 0 ){

if(strcmp(Sigma_independent, "yes") == 0)
{ 
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[3],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[3],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[3],4./14.);}
}
else
{
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[2]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[2]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[2]/(1.+ffactor),4./14.);}
}


if(modeP0==0){parameters1[2]=Sigmanl2;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP2==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP4==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;offset_ini=4;}
if(modeP0+modeP2+modeP4==3){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;parameters1[4]=Sigmanl4;offset_ini=5;}

}
if( strcmp(type_of_analysis, "BAOANISO") == 0 ){

parameters1[2]=parameters2[2];
if(strcmp(Sigma_independent, "yes") == 0)
{
parameters1[3]=parameters2[3];
parameters1[4]=parameters2[4];
}
else
{
parameters1[3]=parameters2[2]/(1+ffactor);
parameters1[4]=parameters2[3];
}
offset_ini=5;

}
}
}


if( strcmp(type_of_analysis, "BAOISO") == 0 ){for(i=offset_ini;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){parameters1[i]=parameters2[i-offset_ini+offset];}}
if( strcmp(type_of_analysis, "BAOANISO") == 0 ){for(i=offset_ini;i<1+(Npolynomial)*(modeP0+modeP2+modeP4)+offset_ini;i++){parameters1[i]=parameters2[i-offset_ini+offset];}}



if( strcmp(type_of_analysis, "BAOISO") == 0 ){


do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0,k2,k4, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],k0[0],k0[NeffP0-1],k2[0],k2[NeffP2-1], k4[0],k4[NeffP4-1], 1,spacing_dataNGC,spacing_theory,Sigma_smooth);

//exit(0);
}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){

//for(i=0;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){printf("%d %lf\n",i,parameters1[i]);}

do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin, Sigma_smooth, pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0,k2,k4, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],k0[0],k0[NeffP0-1],k2[0],k2[NeffP2-1], k4[0],k4[NeffP4-1], 1,spacing_dataNGC,spacing_theory);
//printf("%lf %lf %lf\n",k_theo[300],P_theo0[300],P_theo2[300]);
//printf("%lf %lf\n",P_theo0[10],P_theo2[10]);
//exit(0);
}

free(parameters1);

if(strcmp(path_to_mask1, "none") == 0 )
{
i=-1;

        if(modeP0==1){
        for(j=0;j<NeffP0;j++){
        ptheo=P_theo0[j];
        pobs=P0[j];
        i++;
        difference[i]=ptheo-pobs;//printf("P0 %lf %lf %lf %lf\n",k0[j],k_theo0[j],pobs,ptheo);
        }
        }

        if(modeP2==1){
        for(j=0;j<NeffP2;j++){
        ptheo=P_theo2[j];
        pobs=P2[j];
        i++;
        difference[i]=ptheo-pobs;//printf("P2 %lf %lf %lf %lf\n",k2[j],k_theo2[j],pobs,ptheo);
        }
        }

        if(modeP4==1){
        for(j=0;j<NeffP4;j++){
        ptheo=P_theo4[j];
        pobs=P4[j];
        i++;
        difference[i]=ptheo-pobs;//printf("P4 %lf %lf %lf %lf\n",k4[j],k_theo4[j],pobs,ptheo);
        }
        }
}
else{

j=0;
    for(i=0;i<Ncov;i++)
    {

        if(modeP4==1 && modeP0==0 && modeP2==0){
//        ptheo=P_interpol(k4[j],k_theo,P_theo4,Neffmax*factor_sampling_mask);//printf("P4 %lf %lf\n",P4[j],ptheo);
Ninterpol=determine_N_singlearray(k_theo,k4[j],Neffmax*factor_sampling_mask,spacing_dataNGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k4[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso(P_theo4,Ninterpol,w1);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k4[j],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k4[j],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k4[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo4,Ninterpol,w0,w1,w2);
}
ptheo=P_interpol_fast(k4[j],P_theo4,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
}



        pobs=P4[j];
        j++;
        if(j==NeffP4){j=0;modeP4=0;}
        }

        if(modeP2==1 && modeP0==0){
//        ptheo=P_interpol(k2[j],k_theo,P_theo2,Neffmax*factor_sampling_mask);//printf("P2 %lf %lf\n",P2[j],ptheo);
Ninterpol=determine_N_singlearray(k_theo,k2[j],Neffmax*factor_sampling_mask,spacing_dataNGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k2[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso(P_theo2,Ninterpol,w1);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k2[j],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k2[j],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k2[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo2,Ninterpol,w0,w1,w2);
}
ptheo=P_interpol_fast(k2[j],P_theo2,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
}



        pobs=P2[j];
        j++;
        if(j==NeffP2){j=0;modeP2=0;}
        }

        if(modeP0==1){
//        ptheo=P_interpol(k0[j],k_theo,P_theo0,Neffmax*factor_sampling_mask);//printf("P0 %lf %lf\n",P0[j],ptheo);
Ninterpol=determine_N_singlearray(k_theo,k0[j],Neffmax*factor_sampling_mask,spacing_dataNGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso(P_theo0,Ninterpol,w1);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo0,Ninterpol,w0,w1,w2);
}
ptheo=P_interpol_fast(k0[j],P_theo0,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
}


        pobs=P0[j];
        j++;
        if(j==NeffP0){j=0;modeP0=0;}
        }

        difference[i]=ptheo-pobs;
   }

}
//exit(0);

ch2=0;
//ch2_diag=0;
for(i=0;i<Ncov;i++)
{

   for(j=0;j<Ncov;j++)
      {
       
              ch2=ch2+difference[i]*1./cov[i+Ncov*j]*difference[j];
      }

}
//printf("chi2 bao=%lf\n",ch2);
//exit(0);

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

if(strcmp(path_to_mask1, "none") != 0 ){free(k_theo);}
free(difference);
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo0);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo0);}
}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo2);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo2);}
}   
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo4);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo4);}
}

if(Nchunks==2)
{
factor_sampling_mask=factor_sampling_mask_in;

modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

Ncov=NeffP0SGC*modeP0+NeffP2SGC*modeP2+NeffP4SGC*modeP4;
points=points+Ncov;

if(strcmp(path_to_mask2, "none") == 0 )
{
if(NeffP0SGC*modeP0>0){k_theo0 = (double*) calloc( NeffP0SGC, sizeof(double));}
if(NeffP2SGC*modeP2>0){k_theo2 = (double*) calloc( NeffP2SGC, sizeof(double));}
if(NeffP4SGC*modeP4>0){k_theo4 = (double*) calloc( NeffP4SGC, sizeof(double));}
if(NeffP0SGC*modeP0>0){P_theo0 = (double*) calloc( NeffP0SGC, sizeof(double));}
if(NeffP2SGC*modeP2>0){P_theo2 = (double*) calloc( NeffP2SGC, sizeof(double));}
if(NeffP4SGC*modeP4>0){P_theo4 = (double*) calloc( NeffP4SGC, sizeof(double));}
}
else{

/*
if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){Neffmax=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){Neffmax=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){Neffmax=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){Neffmax=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){Neffmax=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){Neffmax=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
*/

if( strcmp(spacing_dataSGC,"linear") == 0  ){

if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){Neffmax=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){Neffmax=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){Neffmax=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){Neffmax=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){Neffmax=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){Neffmax=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}

}

if( strcmp(spacing_dataSGC,"log") == 0  ){

if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){Neffmax=NeffP0SGC-2+25+(int)(log(k0SGC[0])/(log(k0SGC[NeffP0SGC-1])-log(k0SGC[0]))*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){Neffmax=NeffP0SGC-2+25+(int)(log(k0SGC[0])/(log(k0SGC[NeffP0SGC-1])-log(k0SGC[0]))*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){Neffmax=NeffP2SGC-2+25+(int)(log(k2SGC[0])/(log(k2SGC[NeffP2SGC-1])-log(k2SGC[0]))*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){Neffmax=NeffP2SGC-2+25+(int)(log(k2SGC[0])/(log(k2SGC[NeffP2SGC-1])-log(k2SGC[0]))*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){Neffmax=NeffP4SGC-2+25+(int)(log(k4SGC[0])/(log(k4SGC[NeffP4SGC-1])-log(k4SGC[0]))*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){Neffmax=NeffP4SGC-2+25+(int)(log(k4SGC[0])/(log(k4SGC[NeffP4SGC-1])-log(k4SGC[0]))*(NeffP4SGC-1));}

}

if( strcmp(spacing_dataSGC,"log10") == 0  ){

if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){Neffmax=NeffP0SGC-2+25+(int)(log10(k0SGC[0])/(log10(k0SGC[NeffP0SGC-1])-log10(k0SGC[0]))*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){Neffmax=NeffP0SGC-2+25+(int)(log10(k0SGC[0])/(log10(k0SGC[NeffP0SGC-1])-log10(k0SGC[0]))*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){Neffmax=NeffP2SGC-2+25+(int)(log10(k2SGC[0])/(log10(k2SGC[NeffP2SGC-1])-log10(k2SGC[0]))*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){Neffmax=NeffP2SGC-2+25+(int)(log10(k2SGC[0])/(log10(k2SGC[NeffP2SGC-1])-log10(k2SGC[0]))*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){Neffmax=NeffP4SGC-2+25+(int)(log10(k4SGC[0])/(log10(k4SGC[NeffP4SGC-1])-log10(k4SGC[0]))*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){Neffmax=NeffP4SGC-2+25+(int)(log10(k4SGC[0])/(log10(k4SGC[NeffP4SGC-1])-log10(k4SGC[0]))*(NeffP4SGC-1));}

}

if( strcmp(spacing_dataSGC,"irregular") == 0  ){
Neffmax=Ntheo;
factor_sampling_mask=1;
}


k_theo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0==1){P_theo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2==1){P_theo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4==1){P_theo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
}

difference = (double*) calloc( Ncov, sizeof(double));

if( strcmp(type_of_analysis, "BAOISO") == 0 ){
if(modeP0+modeP2+modeP4==1){Nalphas_for_param1=1;}
else{Nalphas_for_param1=2;}
Nsigmas_for_param1=modeP0+modeP2+modeP4;
parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial+1)+Nalphas_for_param1+Nsigmas_for_param1, sizeof(double));
offset=Nalphas+Nsigmas_tot;
}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){
Nalphas_for_param1=2;
Nsigmas_for_param1=2;;
parameters1 =  (double*) calloc( (modeP0+modeP2+modeP4)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1, sizeof(double));
offset=Nalphas+Nsigmas_tot+1;
}


if(strcmp(Sigma_def_type, "effective") == 0)
{

if(modeP0+modeP2+modeP4==1){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];offset_ini=2;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[2];parameters1[4]=parameters2[2];offset_ini=5;}
if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];offset_ini=4;}
if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2[0];parameters1[1]=parameters2[1];parameters1[2]=parameters2[2];parameters1[3]=parameters2[3];parameters1[4]=parameters2[4];offset_ini=5;}


}
else
{


if(modeP0+modeP2+modeP4==1)
{

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[2],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[2],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[2],4./14.);}
}
else
{
if(modeP0==1){Sigmanl0=pow(parameters2[1],2./6.)*pow(parameters2[1]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[1],6./10.)*pow(parameters2[1]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[1],10./14.)*pow(parameters2[1]/(1.+ffactor),4./14.);}
}

     parameters1[0]=parameters2[0];
if(modeP0==1){parameters1[1]=Sigmanl0;}
if(modeP2==1){parameters1[1]=Sigmanl2;}
if(modeP4==1){parameters1[1]=Sigmanl4;}
offset_ini=2;
}
else{

     parameters1[0]=parameters2[0];
     parameters1[1]=parameters2[1];

if( strcmp(type_of_analysis, "BAOISO") == 0 ){

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[3],4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[3],4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[3],4./14.);}
}else
{
if(modeP0==1){Sigmanl0=pow(parameters2[2],2./6.)*pow(parameters2[2]/(1.+ffactor),4./6.);}
if(modeP2==1){Sigmanl2=pow(parameters2[2],6./10.)*pow(parameters2[2]/(1.+ffactor),4./10.);}
if(modeP4==1){Sigmanl4=pow(parameters2[2],10./14.)*pow(parameters2[2]/(1.+ffactor),4./14.);}
}


if(modeP0==0){parameters1[2]=Sigmanl2;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP2==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl4;offset_ini=4;}
if(modeP4==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;offset_ini=4;}
if(modeP0+modeP2+modeP4==3){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;parameters1[4]=Sigmanl4;offset_ini=5;}


}
if( strcmp(type_of_analysis, "BAOANISO") == 0 ){

parameters1[2]=parameters2[2];
if(strcmp(Sigma_independent, "yes") == 0)
{
parameters1[3]=parameters2[3];
parameters1[4]=parameters2[4];
}
else
{
parameters1[3]=parameters2[2]/(1+ffactor);
parameters1[4]=parameters2[3];
}
offset_ini=5;


}


}

}



if( strcmp(type_of_analysis, "BAOISO") == 0 ){for(i=offset_ini;i<(Npolynomial+1)*(modeP0+modeP2+modeP4)+offset_ini;i++){parameters1[i]=parameters2[i+(Npolynomial+1)*(modeP0+modeP2+modeP4)-offset_ini+offset];}}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){for(i=offset_ini;i<(Npolynomial)*(modeP0+modeP2+modeP4)+1+offset_ini;i++){parameters1[i]=parameters2[i+1+(Npolynomial)*(modeP0+modeP2+modeP4)-offset_ini+offset];}}



if( strcmp(type_of_analysis, "BAOISO") == 0 ){
do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, spacing_maskSGC,path_to_mask2,k0SGC,k2SGC,k4SGC, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],k0SGC[0],k0SGC[NeffP0SGC-1],k2SGC[0],k2SGC[NeffP2SGC-1], k4SGC[0],k4SGC[NeffP4SGC-1], 1,spacing_dataSGC,spacing_theory,Sigma_smooth);
}

if( strcmp(type_of_analysis, "BAOANISO") == 0 ){
//for(i=offset_ini;i<(Npolynomial)*(modeP0+modeP2+modeP4)+1+offset_ini;i++){printf("%d %lf\n",i,parameters1[i]);}

do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,k_Plin,Plin,Nlin, k_Olin, Olin, NOlin, Sigma_smooth,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, spacing_maskSGC, path_to_mask2,k0SGC,k2SGC,k4SGC, Npolynomial, plan1, plan2, k_Plin[0], k_Plin[Nlin-1],k0SGC[0],k0SGC[NeffP0SGC-1],k2SGC[0],k2SGC[NeffP2SGC-1], k4SGC[0],k4SGC[NeffP4SGC-1], 1,spacing_dataSGC,spacing_theory);
//exit(0);
}

free(parameters1);

if(strcmp(path_to_mask2, "none") == 0 )
{
i=-1;

        if(modeP0==1){
        for(j=0;j<NeffP0SGC;j++){
        ptheo=P_theo0[j];
        pobs=P0SGC[j];
        i++;
        difference[i]=ptheo-pobs;
        }
        }

        if(modeP2==1){
        for(j=0;j<NeffP2SGC;j++){
        ptheo=P_theo2[j];
        pobs=P2SGC[j];
        i++;
        difference[i]=ptheo-pobs;
        }
        }

        if(modeP4==1){
        for(j=0;j<NeffP4SGC;j++){
        ptheo=P_theo4[j];
        pobs=P4SGC[j];
        i++;
        difference[i]=ptheo-pobs;
        }
        }
}
else{


j=0;
    for(i=0;i<Ncov;i++)
    {

        if(modeP4==1 && modeP0==0 && modeP2==0){
//        ptheo=P_interpol(k4SGC[j],k_theo,P_theo4,Neffmax*factor_sampling_mask);
Ninterpol=determine_N_singlearray(k_theo,k4SGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4SGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k4SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso(P_theo4,Ninterpol,w1);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k4SGC[j],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k4SGC[j],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k4SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo4,Ninterpol,w0,w1,w2);
}
ptheo=P_interpol_fast(k4SGC[j],P_theo4,Neffmax*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
}


        pobs=P4SGC[j];
        j++;
        if(j==NeffP4SGC){j=0;modeP4=0;}
        }

        if(modeP2==1 && modeP0==0){
//        ptheo=P_interpol(k2SGC[j],k_theo,P_theo2,Neffmax*factor_sampling_mask);
Ninterpol=determine_N_singlearray(k_theo,k2SGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2SGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k2SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso(P_theo2,Ninterpol,w1);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k2SGC[j],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k2SGC[j],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k2SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo2,Ninterpol,w0,w1,w2);
}
ptheo=P_interpol_fast(k2SGC[j],P_theo2,Neffmax*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
}

        pobs=P2SGC[j];//printf("%lf %lf %lf\n",k2SGC[j],P2SGC[j],ptheo);
        j++;
        if(j==NeffP2SGC){j=0;modeP2=0;}
        }

        if(modeP0==1){
//        ptheo=P_interpol(k0SGC[j],k_theo,P_theo0,Neffmax*factor_sampling_mask);
Ninterpol=determine_N_singlearray(k_theo,k0SGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0SGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso(P_theo0,Ninterpol,w1);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo0,Ninterpol,w0,w1,w2);
}
ptheo=P_interpol_fast(k0SGC[j],P_theo0,Neffmax*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
}


        pobs=P0SGC[j];//printf("%lf %lf %lf\n",k0SGC[j],P0SGC[j],ptheo);
        j++;
        if(j==NeffP0SGC){j=0;modeP0=0;}
        }

        difference[i]=ptheo-pobs;
   }


}

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

ch2SGC=0;
for(i=0;i<Ncov;i++)
{

   for(j=0;j<Ncov;j++)
      {
                  ch2SGC=ch2SGC+difference[i]*1./covSGC[i+Ncov*j]*difference[j];
//printf("%d %d, %lf %lf %e\n",i,j,difference[i],difference[j],1./covSGC[i+Ncov*j]);

      }

}

if(strcmp(path_to_mask2, "none") != 0 ){free(k_theo);}
free(difference);
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo0);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo0);}
}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo2);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo2);}
}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo4);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo4);}
}



//printf("%lf, %lf\n",ch2,ch2SGC);
ch2=ch2+ch2SGC;
}


for(i=0;i<Nsigmas_tot;i++)
{

if(Sigma_type[i]==1){prior_chi2=gauss(parameters2[Nalphas+i],Sigma_nl_mean[i],Sigma_nl_stddev[i]);ch2=ch2+prior_chi2;}

}

//printf("%lf %lf %d\n",ch2,ch2SGC,Ncov);
//exit(0);

return ch2;
}



void do_bao_analytic(char *type_BAO_fit,char *type_of_analysis,char *fit_BAO,double *k_Plin, double *Plin, int Nlin, double *k_Olin, double *Olin, int NOlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC,char *path_to_mask2, char *spacing_maskSGC,  double *k0, double *P0, double *errP0,  int NeffP0, double *k2, double *P2, double *errP2,  int NeffP2,  double *k4, double *P4, double *errP4,  int NeffP4, double *k0SGC, double *P0SGC, double *errP0SGC, int NeffP0SGC, double *k2SGC, double *P2SGC, double *errP2SGC, int NeffP2SGC,double *k4SGC, double *P4SGC, double *errP4SGC, int NeffP4SGC, double *cov, double *covSGC, double alpha_min, double alpha_max, double alpha_step, char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks, char *path_output, char *identifier, char *do_plot, char *do_power_spectrum, char *do_bispectrum,char *spacing_dataNGC,char *spacing_dataSGC, char *spacing_theory)
{
int modeP0,modeP2,modeP4;
long int l,j,i,i1,i2,imax,i_min,imax2;
int c1,c2;
int l1,l2,Ncov,NcovSGC;
double chi2_min,CHI2;
double alpha_input,alpha_input2;
imax=(int)((alpha_max-alpha_min)/alpha_step)+2;
if(imax<=0){printf("Wrong value of imax in function do_bao_analytic: %ld. Exiting now...\n",imax);exit(0);}
double *parameters1,*parameters2;
double *parameters1_pre,*parameters1_post;
double *parameters2_pre,*parameters2_post;
double *cov0,*cov0SGC,*cov2,*cov2SGC,*cov4,*cov4SGC;
int dimension;
double **parameter_output;
int Npoints,Ndof;
FILE *f;
char name_output[2000];
double prior;
double epsilon;
double sigma0,sigma_trial,step,direction,chi2_0,chi2_trial,convergence;
int iter,icheck;
int Nplan;
double params[5];
double kmax_Nlin,kmin_Nlin;
double kmax_data,kmin_data;
int processors;
double time_run, time_ini, time_fin;
int Nalphas,Nsigmas_free,Nsigmas_tot;
//int allsigmafixed;
int factor_sampling_mask=10;
double Sigma_nl_mean_P0,Sigma_nl_mean_P2,Sigma_nl_mean_P4;
time_ini=time(NULL);

modeP0=0;
modeP2=0;
modeP4=0;
Npoints=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;Npoints=Npoints+NeffP0;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;Npoints=Npoints+NeffP2;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;Npoints=Npoints+NeffP4;}

epsilon=1e-10;

if(Nchunks==2){
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){Npoints=Npoints+NeffP0SGC;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){Npoints=Npoints+NeffP2SGC;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){Npoints=Npoints+NeffP4SGC;}
}

//if( strcmp(Sigma_nl_type, "fixed") == 0 ){Ndof=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+1;}
//if( strcmp(Sigma_nl_type, "fixed") != 0 ){Ndof=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+2;}

//if(modeP0+modeP2+modeP4>1){Ndof=Ndof+1;}

//allsigmafixed=-1;
Nalphas=1;if(modeP0+modeP2+modeP4>1){Nalphas=2;}
Nsigmas_free=0;
if(strcmp(Sigma_independent, "yes") == 0 ){

    if(strcmp(Sigma_def_type, "para-perp") == 0)
    {
       if(Sigma_type[0]>0){Nsigmas_free=Nsigmas_free+1;}
       if(Sigma_type[1]>0){Nsigmas_free=Nsigmas_free+1;}
    }

    if(strcmp(Sigma_def_type, "effective") == 0)
    {
        if(modeP0==1 && Sigma_type[0]>0){Nsigmas_free=Nsigmas_free+1;}
        if(modeP2==1 && Sigma_type[1]>0){Nsigmas_free=Nsigmas_free+1;}
        if(modeP4==1 && Sigma_type[2]>0){Nsigmas_free=Nsigmas_free+1;}
    }

}
if(strcmp(Sigma_independent, "no") == 0 ){
//only one sigma (or none if fixed)
if(Sigma_type[0]>0){Nsigmas_free=1;}
if(Sigma_type[0]==0){Nsigmas_free=0;}
}

Ndof=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+Nalphas+Nsigmas_free;//number free  parameters



if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0){Nsigmas_tot=2;}
if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0){Nsigmas_tot=modeP0+modeP2+modeP4;}
if( strcmp(Sigma_independent, "no") == 0 ){Nsigmas_tot=1;}
dimension=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+Nalphas+Nsigmas_tot;//number total parameters


if(dimension<=0){printf("Wrong value of dimension in function do_bao_analytic: %d. Exiting now...\n",dimension);exit(0);}

parameter_output = (double**) calloc(1+dimension,sizeof(double*));

if(modeP0+modeP2+modeP4>1){imax2=imax*imax;}
if(modeP0+modeP2+modeP4==1){imax2=imax;}

for(j=0;j<1+dimension;j++)
{
   parameter_output[j] = (double*) calloc(imax2,sizeof(double));
}

 if(modeP0==1){  cov0= (double*) calloc( NeffP0*NeffP0, sizeof(double));}
 if(modeP2==1){  cov2= (double*) calloc( NeffP2*NeffP2, sizeof(double));}
 if(modeP4==1){  cov4= (double*) calloc( NeffP4*NeffP4, sizeof(double));}

if(strcmp(fit_BAO, "P0") == 0){Ncov=NeffP0;}
if(strcmp(fit_BAO, "P2") == 0){Ncov=NeffP2;}
if(strcmp(fit_BAO, "P4") == 0){Ncov=NeffP4;}
if(strcmp(fit_BAO, "P02") == 0){Ncov=NeffP0+NeffP2;}
if(strcmp(fit_BAO, "P04") == 0){Ncov=NeffP0+NeffP4;}
if(strcmp(fit_BAO, "P24") == 0){Ncov=NeffP2+NeffP4;}
if(strcmp(fit_BAO, "P024") == 0){Ncov=NeffP0+NeffP2+NeffP4;}

for(l1=0;l1<Ncov;l1++)
{
     for(l2=0;l2<Ncov;l2++)
     {
if(strcmp(fit_BAO, "P0") == 0){cov0[l1+NeffP0*l2]=cov[l1+Ncov*l2];}

if(strcmp(fit_BAO, "P2") == 0){cov2[l1+NeffP2*l2]=cov[l1+Ncov*l2];}

if(strcmp(fit_BAO, "P4") == 0){cov4[l1+NeffP4*l2]=cov[l1+Ncov*l2];}


if(strcmp(fit_BAO, "P02") == 0)
{if(l1<NeffP0 && l2<NeffP0)
{c1=l1;c2=l2;
cov0[c1+NeffP0*c2]=cov[l1+Ncov*l2];}
if(l1>NeffP0-1 && l2>NeffP0-1)
{c1=l1-NeffP0;c2=l2-NeffP0;
cov2[c1+NeffP2*c2]=cov[l1+Ncov*l2];}}

if(strcmp(fit_BAO, "P24") == 0)
{if(l1<NeffP2 && l2<NeffP2)
{c1=l1;c2=l2;
cov2[c1+NeffP2*c2]=cov[l1+Ncov*l2];}
if(l1>NeffP2-1 && l2>NeffP2-1)
{c1=l1-NeffP2;c2=l2-NeffP2;
cov4[c1+NeffP4*c2]=cov[l1+Ncov*l2];}}

if(strcmp(fit_BAO, "P04") == 0)
{if(l1<NeffP0 && l2<NeffP0)
{c1=l1;c2=l2;
cov0[c1+NeffP0*c2]=cov[l1+Ncov*l2];}
if(l1>NeffP0-1 && l2>NeffP0-1)
{c1=l1-NeffP0;c2=l2-NeffP0;
cov4[c1+NeffP4*c2]=cov[l1+Ncov*l2];}}

if(strcmp(fit_BAO, "P024") == 0)
{if(l1<NeffP0 && l2<NeffP0)
{c1=l1;c2=l2;
cov0[c1+NeffP0*c2]=cov[l1+Ncov*l2];}
if(l1<NeffP2+NeffP0 && l2<NeffP2+NeffP0 && l1>NeffP0-1 && l2>NeffP0-1)
{c1=l1-NeffP0;c2=l2-NeffP0;
cov2[c1+NeffP2*c2]=cov[l1+Ncov*l2];}

if(l1>NeffP0+NeffP2-1 && l2>NeffP0+NeffP2-1)
{c1=l1-NeffP0-NeffP2;c2=l2-NeffP0-NeffP2;
cov4[c1+NeffP4*c2]=cov[l1+Ncov*l2];}}

}
}

if(Nchunks==2){
 if(modeP0==1){  cov0SGC= (double*) calloc( NeffP0SGC*NeffP0SGC, sizeof(double));}
 if(modeP2==1){  cov2SGC= (double*) calloc( NeffP2SGC*NeffP2SGC, sizeof(double));}
 if(modeP4==1){  cov4SGC= (double*) calloc( NeffP4SGC*NeffP4SGC, sizeof(double));}

if(strcmp(fit_BAO, "P0") == 0){NcovSGC=NeffP0SGC;}
if(strcmp(fit_BAO, "P2") == 0){NcovSGC=NeffP2SGC;}
if(strcmp(fit_BAO, "P4") == 0){NcovSGC=NeffP4SGC;}
if(strcmp(fit_BAO, "P02") == 0){NcovSGC=NeffP0SGC+NeffP2SGC;}
if(strcmp(fit_BAO, "P04") == 0){NcovSGC=NeffP0SGC+NeffP4SGC;}
if(strcmp(fit_BAO, "P24") == 0){NcovSGC=NeffP2SGC+NeffP4SGC;}
if(strcmp(fit_BAO, "P024") == 0){NcovSGC=NeffP0SGC+NeffP2SGC+NeffP4SGC;}


for(l1=0;l1<NcovSGC;l1++)
{
     for(l2=0;l2<NcovSGC;l2++)
     {
if(strcmp(fit_BAO, "P0") == 0){cov0SGC[l1+NeffP0SGC*l2]=covSGC[l1+NcovSGC*l2];}

if(strcmp(fit_BAO, "P2") == 0){cov2SGC[l1+NeffP2SGC*l2]=covSGC[l1+NcovSGC*l2];}

if(strcmp(fit_BAO, "P4") == 0){cov4SGC[l1+NeffP4SGC*l2]=covSGC[l1+NcovSGC*l2];}

if(strcmp(fit_BAO, "P02") == 0)
{if(l1<NeffP0SGC && l2<NeffP0SGC)
{c1=l1;c2=l2;
cov0SGC[c1+NeffP0SGC*c2]=covSGC[l1+NcovSGC*l2];}
if(l1>NeffP0SGC-1 && l2>NeffP0SGC-1)
{c1=l1-NeffP0SGC;c2=l2-NeffP0SGC;
cov2SGC[c1+NeffP2SGC*c2]=covSGC[l1+NcovSGC*l2];}}

if(strcmp(fit_BAO, "P24") == 0)
{if(l1<NeffP2SGC && l2<NeffP2SGC)
{c1=l1;c2=l2;
cov2SGC[c1+NeffP2SGC*c2]=covSGC[l1+NcovSGC*l2];}
if(l1>NeffP2SGC-1 && l2>NeffP2SGC-1)
{c1=l1-NeffP2SGC;c2=l2-NeffP2SGC;
cov4SGC[c1+NeffP4SGC*c2]=covSGC[l1+NcovSGC*l2];}}

if(strcmp(fit_BAO, "P04") == 0)
{if(l1<NeffP0SGC && l2<NeffP0SGC)
{c1=l1;c2=l2;
cov0SGC[c1+NeffP0SGC*c2]=covSGC[l1+NcovSGC*l2];}
if(l1>NeffP0SGC-1 && l2>NeffP0SGC-1)
{c1=l1-NeffP0SGC;c2=l2-NeffP0SGC;
cov4SGC[c1+NeffP4SGC*c2]=covSGC[l1+NcovSGC*l2];}}

if(strcmp(fit_BAO, "P024") == 0)
{if(l1<NeffP0SGC && l2<NeffP0SGC)
{c1=l1;c2=l2;
cov0SGC[c1+NeffP0SGC*c2]=covSGC[l1+NcovSGC*l2];}
if(l1<NeffP2SGC+NeffP0SGC && l2<NeffP2SGC+NeffP0SGC && l1>NeffP0SGC-1 && l2>NeffP0SGC-1)
{c1=l1-NeffP0SGC;c2=l2-NeffP0SGC;
cov2SGC[c1+NeffP2SGC*c2]=covSGC[l1+NcovSGC*l2];}

if(l1>NeffP0SGC+NeffP2SGC-1 && l2>NeffP0SGC+NeffP2SGC-1)
{c1=l1-NeffP0SGC-NeffP2SGC;c2=l2-NeffP0SGC-NeffP2SGC;
cov4SGC[c1+NeffP4SGC*c2]=covSGC[l1+NcovSGC*l2];}}

}
}



}


// FFTW plans ARE NOT thread-safe and MUST BE created BEFORE the multi-thread starts
 
kmin_Nlin=k_Plin[0];
kmax_Nlin=k_Plin[Nlin-1];

kmax_data=-1;
kmin_data=999999;
if(modeP4==1)
{
if(kmin_data>k4[0]){kmin_data=k4[0];}
if(kmax_data<k4[NeffP4-1]){kmax_data=k4[NeffP4-1];}
if(Nchunks==2){
if(kmin_data>k4SGC[0]){kmin_data=k4SGC[0];}
if(kmax_data<k4SGC[NeffP4SGC-1]){kmax_data=k4SGC[NeffP4SGC-1];}}
}
if(modeP2==1)
{
if(kmin_data>k2[0]){kmin_data=k2[0];}
if(kmax_data<k2[NeffP2-1]){kmax_data=k2[NeffP2-1];}
if(Nchunks==2){
if(kmin_data>k2SGC[0]){kmin_data=k2SGC[0];}
if(kmax_data<k2SGC[NeffP2SGC-1]){kmax_data=k2SGC[NeffP2SGC-1];}}
}
if(modeP0==1)
{
if(kmin_data>k0[0]){kmin_data=k0[0];}
if(kmax_data<k0[NeffP0-1]){kmax_data=k0[NeffP0-1];}
if(Nchunks==2){
if(kmin_data>k0SGC[0]){kmin_data=k0SGC[0];}
if(kmax_data<k0SGC[NeffP0SGC-1]){kmax_data=k0SGC[NeffP0SGC-1];}}
}
if(kmin_data==999999 || kmax_data==-1){printf("Error with values of kmin,kmax for data. Exiting now...\n");exit(0);}

set_mask_params(params,kmin_Nlin,kmax_Nlin,Nlin,kmin_data,kmax_data);
Nplan=(int)(params[2]);

fftw_complex *a_pointer;
fftw_complex *b_pointer;

    fftw_plan plan1 = fftw_plan_dft_1d(Nplan,  a_pointer,  b_pointer,  -1, FFTW_ESTIMATE);//forward plan
    fftw_plan plan2 = fftw_plan_dft_1d(Nplan,  b_pointer,  b_pointer, +1, FFTW_ESTIMATE);//reverse plan
   

#pragma omp parallel for  private(Sigma_nl_mean_P0,Sigma_nl_mean_P2,Sigma_nl_mean_P4,icheck,l,j,i,i1,i2,alpha_input,alpha_input2,parameters1,parameters1_pre,parameters1_post,parameters2,parameters2_pre,parameters2_post,CHI2,prior,sigma0,step,direction,iter,sigma_trial,chi2_0,chi2_trial,convergence) shared(imax,imax2,parameter_output,dimension,k_Plin,Plin,Nlin,k_Olin,Olin,NOlin,pos,W0,W2,W4,W6,W8,Nmask,path_to_mask1,posSGC,W0SGC,W2SGC,W4SGC,W6SGC,W8SGC,NmaskSGC,path_to_mask2,k0,P0,k2,P2,k4,P4,k0SGC,P0SGC,k2SGC,P2SGC,k4SGC,P4SGC,NeffP0,NeffP0SGC,NeffP2,NeffP2SGC,NeffP4,NeffP4SGC,cov,covSGC,alpha_min,alpha_max,alpha_step,Npolynomial,Nchunks,path_output,identifier,do_plot,epsilon,plan1,plan2,type_of_analysis,fit_BAO,modeP0,modeP2,modeP4,cov0,cov2,cov4,cov0SGC,cov2SGC,cov4SGC, do_power_spectrum, do_bispectrum,type_BAO_fit,processors, Sigma_nl_mean, Sigma_nl_stddev, Sigma_type, Sigma_def_type, Sigma_independent, ffactor,Nalphas,Nsigmas_tot,Nsigmas_free,factor_sampling_mask)
for(i=0;i<imax2;i++)
{

parameters1 =  (double*) calloc( Npolynomial+1, sizeof(double));
parameters2 = (double*) calloc( dimension, sizeof(double));

if(i==0){processors=omp_get_num_threads();}

if(modeP0+modeP2+modeP4==1){alpha_input=alpha_min+i*alpha_step;}
if(modeP0+modeP2+modeP4>1){
i1=(long int)(i*1./imax*1.);
i2=i-imax*i1;

//alpha_input=0.996229;
//alpha_input2=1.009689;
alpha_input=alpha_min+i1*alpha_step;
alpha_input2=alpha_min+i2*alpha_step;

}

parameters2[0]=alpha_input;
if(modeP0+modeP2+modeP4>1){parameters2[1]=alpha_input2;}


    if( Nsigmas_free == 0 )//(any) sigma cannot be a free variable in the analytic analysis
    {

if(modeP0+modeP2+modeP4==1){

  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
  parameters2[1]=Sigma_nl_mean[0];
  parameters2[2]=Sigma_nl_mean[1];
Sigma_nl_mean_P0=pow(Sigma_nl_mean[0],2./6.)*pow(Sigma_nl_mean[1],4./6.);
Sigma_nl_mean_P2=pow(Sigma_nl_mean[0],6./10.)*pow(Sigma_nl_mean[1],4./10.);
Sigma_nl_mean_P4=pow(Sigma_nl_mean[0],10./14.)*pow(Sigma_nl_mean[1],4./14.);
  }

  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
  parameters2[1]=Sigma_nl_mean[0];
Sigma_nl_mean_P0=pow(Sigma_nl_mean[0],2./6.)*pow(Sigma_nl_mean[0]/(1.+ffactor),4./6.);
Sigma_nl_mean_P2=pow(Sigma_nl_mean[0],6./10.)*pow(Sigma_nl_mean[0]/(1.+ffactor),4./10.);
Sigma_nl_mean_P4=pow(Sigma_nl_mean[0],10./14.)*pow(Sigma_nl_mean[0]/(1.+ffactor),4./14.);
  }


  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
    if(modeP0==1){parameters2[1]=Sigma_nl_mean[0];Sigma_nl_mean_P0=Sigma_nl_mean[0];}
    if(modeP2==1){parameters2[1]=Sigma_nl_mean[1];Sigma_nl_mean_P2=Sigma_nl_mean[1];}
    if(modeP4==1){parameters2[1]=Sigma_nl_mean[2];Sigma_nl_mean_P4=Sigma_nl_mean[2];}

  }


  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
    if(modeP0==1){parameters2[1]=Sigma_nl_mean[0];Sigma_nl_mean_P0=Sigma_nl_mean[0];}
    if(modeP2==1){parameters2[1]=Sigma_nl_mean[0];Sigma_nl_mean_P2=Sigma_nl_mean[0];}
    if(modeP4==1){parameters2[1]=Sigma_nl_mean[0];Sigma_nl_mean_P4=Sigma_nl_mean[0];}

  }



  }
if(modeP0+modeP2+modeP4>1){

  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=Sigma_nl_mean[1];

  Sigma_nl_mean_P0=pow(Sigma_nl_mean[0],2./6.)*pow(Sigma_nl_mean[1],4./6.);
  Sigma_nl_mean_P2=pow(Sigma_nl_mean[0],6./10.)*pow(Sigma_nl_mean[1],4./10.);
  Sigma_nl_mean_P4=pow(Sigma_nl_mean[0],10./14.)*pow(Sigma_nl_mean[1],4./14.);

  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];

  Sigma_nl_mean_P0=pow(Sigma_nl_mean[0],2./6.)*pow(Sigma_nl_mean[0]/(1.+ffactor),4./6.);
  Sigma_nl_mean_P2=pow(Sigma_nl_mean[0],6./10.)*pow(Sigma_nl_mean[0]/(1.+ffactor),4./10.);
  Sigma_nl_mean_P4=pow(Sigma_nl_mean[0],10./14.)*pow(Sigma_nl_mean[0]/(1.+ffactor),4./14.);

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
    if(modeP0==1 && modeP2==1 && modeP4==1){parameters2[2]=Sigma_nl_mean[0];Sigma_nl_mean_P0=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];Sigma_nl_mean_P2=Sigma_nl_mean[1];parameters2[4]=Sigma_nl_mean[2];Sigma_nl_mean_P4=Sigma_nl_mean[2];}
    if(modeP0==1 && modeP2==1 && modeP4==0){parameters2[2]=Sigma_nl_mean[0];Sigma_nl_mean_P0=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];Sigma_nl_mean_P2=Sigma_nl_mean[1];}
    if(modeP0==1 && modeP2==0 && modeP4==1){parameters2[2]=Sigma_nl_mean[0];Sigma_nl_mean_P0=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[2];Sigma_nl_mean_P4=Sigma_nl_mean[2];}
    if(modeP0==0 && modeP2==1 && modeP4==1){parameters2[2]=Sigma_nl_mean[1];Sigma_nl_mean_P2=Sigma_nl_mean[1];parameters2[3]=Sigma_nl_mean[2];Sigma_nl_mean_P4=Sigma_nl_mean[2];}

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
     parameters2[2]=Sigma_nl_mean[0];
     Sigma_nl_mean_P0=Sigma_nl_mean[0];
     Sigma_nl_mean_P2=Sigma_nl_mean[0];
     Sigma_nl_mean_P4=Sigma_nl_mean[0];
  }



}

      
if(modeP0==1){
if(NeffP0<=0){printf("Error NeffP0=0. Exiting now...\n");exit(0);}

          parameter_bestfit(0,modeP0,modeP2,modeP4,alpha_input,alpha_input2, Sigma_nl_mean_P0, k0,P0,cov0, NeffP0, parameters1, Npolynomial, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1, k_Plin, Plin, Nlin, k_Olin, Olin, NOlin,plan1,plan2);

if(modeP0+modeP2+modeP4==1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j]=parameters1[j];}}
if(modeP0+modeP2+modeP4>1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j]=parameters1[j];}}
}

if(modeP2==1){
if(NeffP2<=0){printf("Error NeffP2=0. Exiting now...\n");exit(0);}
          parameter_bestfit(2,modeP0,modeP2,modeP4,alpha_input,alpha_input2, Sigma_nl_mean_P2, k2,P2, cov2,NeffP2, parameters1, Npolynomial, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1, k_Plin, Plin, Nlin, k_Olin, Olin, NOlin,plan1,plan2);


if(modeP0+modeP2+modeP4==1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j]=parameters1[j];}}
if(modeP0+modeP2+modeP4>1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)*modeP0]=parameters1[j];}}

}

if(modeP4==1){
if(NeffP4<=0){printf("Error NeffP4=0. Exiting now...\n");exit(0);}
          parameter_bestfit(4,modeP0,modeP2,modeP4,alpha_input,alpha_input2, Sigma_nl_mean_P4,k4,P4,cov4,NeffP4, parameters1, Npolynomial, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1, k_Plin, Plin, Nlin, k_Olin, Olin, NOlin,plan1,plan2);



if(modeP0+modeP2+modeP4==1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j]=parameters1[j];}}
if(modeP0+modeP2+modeP4>1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)*(modeP2+modeP0)]=parameters1[j];}}
}

       if(Nchunks==2){

if(modeP0==1){
if(NeffP0SGC<=0){printf("Error NeffP0SGC=0. Exiting now...\n");exit(0);}
          parameter_bestfit(0,modeP0,modeP2,modeP4,alpha_input,alpha_input2, Sigma_nl_mean_P0, k0SGC,P0SGC,cov0SGC, NeffP0SGC, parameters1, Npolynomial, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2, k_Plin, Plin, Nlin, k_Olin, Olin, NOlin,plan1,plan2);

if(modeP0+modeP2+modeP4==1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)]=parameters1[j];}}
if(modeP0+modeP2+modeP4>1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)*(modeP4+modeP2+modeP0)]=parameters1[j];}}
}

if(modeP2==1){
if(NeffP2SGC<=0){printf("Error NeffP2SGC=0. Exiting now...\n");exit(0);}
          parameter_bestfit(2,modeP0,modeP2,modeP4,alpha_input,alpha_input2, Sigma_nl_mean_P2, k2SGC, P2SGC, cov2SGC, NeffP2SGC, parameters1, Npolynomial, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2, k_Plin, Plin, Nlin, k_Olin, Olin, NOlin,plan1,plan2);

if(modeP0+modeP2+modeP4==1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)]=parameters1[j];}}
if(modeP0+modeP2+modeP4>1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)*(modeP4+modeP2+2*modeP0)]=parameters1[j];}}
}

if(modeP4==1){
if(NeffP4SGC<=0){printf("Error NeffP4SGC=0. Exiting now...\n");exit(0);}
          parameter_bestfit(4,modeP0,modeP2,modeP4,alpha_input,alpha_input2, Sigma_nl_mean_P4,k4SGC,P4SGC,cov4SGC, NeffP4SGC, parameters1, Npolynomial, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2, k_Plin, Plin, Nlin, k_Olin, Olin, NOlin,plan1,plan2);

if(modeP0+modeP2+modeP4==1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)]=parameters1[j];}}
if(modeP0+modeP2+modeP4>1){for(j=0;j<Npolynomial+1;j++){parameters2[Nalphas+Nsigmas_tot+j+(Npolynomial+1)*(2*modeP0+2*modeP2+modeP4)]=parameters1[j];}}
}
                    


                    }


 
                           CHI2=chi2_bao(type_BAO_fit,type_of_analysis,fit_BAO,parameters2,k_Plin,Plin,Nlin,k_Olin,Olin,NOlin,pos,W0,W2,W4,W6,W8,Nmask,path_to_mask1,spacing_maskNGC,posSGC,W0SGC,W2SGC,W4SGC,W6SGC,W8SGC,NmaskSGC,path_to_mask2,spacing_maskSGC,k0,P0,NeffP0,k2,P2,NeffP2,k4,P4,NeffP4,k0SGC,P0SGC,NeffP0SGC,k2SGC,P2SGC,NeffP2SGC,k4SGC,P4SGC,NeffP4SGC,cov,covSGC,Sigma_def_type, Sigma_independent, ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev,Npolynomial,Nchunks,plan1,plan2, do_power_spectrum, do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free,0,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory);
    
                    for(l=0;l<dimension;l++){parameter_output[l][i]=parameters2[l];}
                        parameter_output[dimension][i]=CHI2;

//printf("ch2 analytic=%lf, %lf %lf\n",CHI2,Sigma_nl_mean_P0,Sigma_nl_mean_P2);
//exit(0);
    }
   else
   {
    printf("Error. You shouldn't be here.Exiting.../n");
    exit(0);              //TBD
   }

free(parameters1);
free(parameters2);

}


sprintf(name_output,"%s/likelihood_%s.txt",path_output,identifier);
f=fopen(name_output,"w");
fprintf(f,"#");
if(modeP0+modeP2+modeP4==1){
fprintf(f,"B\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A%ld\t",j);}
}
else{

if(modeP0==1){
fprintf(f,"B0\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A0_%ld\t",j);}
}

if(modeP2==1){
fprintf(f,"B2\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A2_%ld\t",j);}
}

if(modeP4==1){
fprintf(f,"B4\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A4_%ld\t",j);}
}

}


if(Nchunks==2){

if(modeP0+modeP2+modeP4==1){
fprintf(f,"BSGC\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A%ldSGC\t",j);}
}
else
{
if(modeP0==1){
fprintf(f,"B0SGC\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A0_%ldSGC\t",j);}
}

if(modeP2==1){
fprintf(f,"B2SGC\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A2_%ldSGC\t",j);}
}

if(modeP4==1){
fprintf(f,"B4SGC\t");
for(j=1;j<=Npolynomial;j++){fprintf(f,"A4_%ldSGC\t",j);}
}


}

}

if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f,"Sigmanl_eff \t alpha0 \t chi2\n");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f,"Sigmanl_eff \t alpha2 \t chi2\n");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f,"Sigmanl_eff \t alpha4 \t chi2\n");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){fprintf(f,"Sigmanl_eff \t alpha_para \t alpha_perp \t chi2\n");}
}
if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f,"Sigmanl0_eff \t alpha0 \t chi2\n");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f,"Sigmanl2_eff \t alpha2 \t chi2\n");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f,"Sigmanl4_eff \t alpha4 \t chi2\n");}
if( strcmp(fit_BAO, "P02") == 0 ){fprintf(f,"Sigmanl0_eff \t Sigmanl2_eff \t alpha_para \t alpha_perp \t chi2\n");}
if( strcmp(fit_BAO, "P04") == 0 ){fprintf(f,"Sigmanl0_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\n");}
if( strcmp(fit_BAO, "P24") == 0 ){fprintf(f,"Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\n");}
if( strcmp(fit_BAO, "P024") == 0 ){fprintf(f,"Sigmanl0_eff \t Sigmanl2_eff \t Sigmanl4_eff \t alpha_para \t alpha_perp \t chi2\n");}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f,"Sigmanl_para \t alpha0 \t chi2\n");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f,"Sigmanl_para \t alpha2 \t chi2\n");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f,"Sigmanl_para \t alpha4 \t chi2\n");}
if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){fprintf(f,"Sigmanl_para \t alpha_para \t alpha_perp \t chi2\n");}
}
if( strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 ){fprintf(f,"Sigmanl_para \t Sigma_perp \t alpha0 \t chi2\n");}
if( strcmp(fit_BAO, "P2") == 0 ){fprintf(f,"Sigmanl_para \t Sigma_perp \t alpha2 \t chi2\n");}
if( strcmp(fit_BAO, "P4") == 0 ){fprintf(f,"Sigmanl_para \t Sigma_perp \t alpha4 \t chi2\n");}
if( strcmp(fit_BAO, "P02") == 0 ||  strcmp(fit_BAO, "P04") == 0 ||  strcmp(fit_BAO, "P24") == 0 ||  strcmp(fit_BAO, "P024") == 0){fprintf(f,"Sigmanl_para \t Sigmanl_perp \t alpha_para \t alpha_perp \t chi2\n");}
}

Ndof=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+1;
if(modeP0+modeP2+modeP4>1){Ndof++;}
for(i=0;i<Nsigmas_tot;i++){
if( Sigma_type[i] >0){Ndof++;}
}

fprintf(f,"#Npoints-dof=%d-%d\n",Npoints,Ndof);


for(i=0;i<imax2;i++)
{

//if(modeP0+modeP2+modeP4==1){for(j=2;j<=dimension-1;j++){fprintf(f,"%e\t",parameter_output[j][i]);}}
//if(modeP0+modeP2+modeP4>1){for(j=3;j<=dimension-1;j++){fprintf(f,"%e\t",parameter_output[j][i]);}}

for(j=Nalphas+Nsigmas_tot;j<=dimension-1;j++){fprintf(f,"%e\t",parameter_output[j][i]);}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "no") == 0 ){
fprintf(f,"%e\t %e\t %e\n",parameter_output[1][i],parameter_output[0][i], parameter_output[dimension][i]);
}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
fprintf(f,"%e\t %e\t %e\n",parameter_output[1][i],parameter_output[0][i], parameter_output[dimension][i]);
}

if(modeP0+modeP2+modeP4==1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){
fprintf(f,"%e\t %e\t %e\t %e\n",parameter_output[1][i], parameter_output[2][i], parameter_output[0][i], parameter_output[dimension][i]);
}



if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "no") == 0 ){
fprintf(f,"%e\t %e\t %e\t %e\n",parameter_output[2][i],parameter_output[0][i],parameter_output[1][i], parameter_output[dimension][i]);
}
if(modeP0+modeP2+modeP4>1 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "para-perp") == 0 ){
fprintf(f,"%e\t %e\t %e\t %e\t %e\n",parameter_output[2][i], parameter_output[3][i], parameter_output[0][i],parameter_output[1][i], parameter_output[dimension][i]);
}

if(modeP0+modeP2+modeP4==2 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
fprintf(f,"%e\t %e\t %e\t %e\t %e\n",parameter_output[2][i], parameter_output[3][i], parameter_output[0][i],parameter_output[1][i], parameter_output[dimension][i]);
}

if(modeP0+modeP2+modeP4==3 && strcmp(Sigma_independent, "yes") == 0 && strcmp(Sigma_def_type, "effective") == 0 ){
fprintf(f,"%e\t %e\t %e\t %e\t %e\t %e\n",parameter_output[2][i], parameter_output[3][i], parameter_output[4][i], parameter_output[0][i], parameter_output[1][i], parameter_output[dimension][i]);
}



if(i==0){chi2_min=parameter_output[dimension][i];i_min=i;}
if(i>0 && chi2_min>parameter_output[dimension][i]){chi2_min=parameter_output[dimension][i];i_min=i;}

}
fclose(f);

    if( strcmp(do_plot, "yes") == 0 )
    {

parameters2 =  (double*) calloc( dimension, sizeof(double));


for(j=0;j<dimension;j++){parameters2[j]=parameter_output[j][i_min];/*printf("%d %lf\n",j,parameters2[j]);*/}


make_a_bao_plot(type_BAO_fit,type_of_analysis,fit_BAO,parameters2,parameter_output[dimension][i_min], k0,P0,errP0,NeffP0,k0SGC,P0SGC,errP0SGC,NeffP0SGC,k2,P2,errP2,NeffP2,k2SGC,P2SGC,errP2SGC,NeffP2SGC,k4,P4,errP4,NeffP4,k4SGC,P4SGC,errP4SGC,NeffP4SGC, k_Plin, Plin, Nlin, k_Olin, Olin, NOlin, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, Sigma_def_type, Sigma_independent, ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks,Npoints,Ndof, path_output, identifier,plan1,plan2,Nalphas,Nsigmas_tot, Nsigmas_free,0,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory);

free(parameters2);
    }

time_fin=time(NULL);
time_run=(time_fin-time_ini)/(60.*60.);//in hours

do_log_file2(path_output,identifier,name_output,fit_BAO,Nchunks,Npolynomial,time_run,processors,alpha_min, alpha_max, alpha_step,Nalphas,Nsigmas_tot);

freeTokens(parameter_output,dimension+1);
fftw_destroy_plan(plan1);
fftw_destroy_plan(plan2);


}


