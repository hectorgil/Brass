#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>
#include "functions.h"
#include "structures.h"
#include "integrals_rsd.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "fftlog.h"
#include "mcmc_bao.h"
#include "cubature.h"

#define Pi (4.*atan(1.))

void do_Ptheo_RSD(char *ptmodel,char *rsdmodel_ps,char *fogmodel_ps, char *type_of_analysis, char *RSD_fit, char *fit_BAO,int modeP0,int modeP2,int modeP4,double **Theory,int N_Plin,  double k_theo[], double k_theo0[], double k_theo2[], double k_theo4[], double P_theo0[], double P_theo2[], double P_theo4[], double Pnoise,int NeffP0, int NeffP2, int NeffP4,int factor_for_sampling, double *parameters1,double *pos, double *W0,double *W2,double *W4,double *W6,double *W8,int Nmask,char *spacing_mask, double *k0, double *k2, double *k4, char *path_to_mask1,fftw_plan plan1,fftw_plan plan2,double k0min, double k0max, double k2min, double k2max, double k4min, double k4max, char *spacing_data,char *spacing_theory)
{

int i,j;
double b1,b2,Anoise,sigma8_scaling,f;
double bs2,b3nl,sigmaP,apara,aperp;
double ptheo,ptheo2,ptheo4,junk;
double XMIN1[1]={-1};
double XMAX1[1]={+1};
double precision=1e-1;
int Nelements;
double kmax,kmin;
int NeffP;
int NeffP_max;
if(k0max>=k2max && modeP0==1 && modeP2==1){kmax=k0max;}
if(k0max>=k4max && modeP0==1 && modeP4==1){kmax=k0max;}

if(k2max>=k0max && modeP0==1 && modeP2==1){kmax=k2max;}
if(k2max>=k4max && modeP4==1 && modeP2==1){kmax=k2max;}

if(k4max>=k0max && modeP0==1 && modeP4==1){kmax=k4max;}
if(k4max>=k2max && modeP4==1 && modeP2==1){kmax=k4max;}

if(k0min<=k2min && modeP0==1 && modeP2==1){kmin=k0min;}
if(k0min<=k4min && modeP0==1 && modeP4==1){kmin=k0min;}

if(k2min<=k0min && modeP0==1 && modeP2==1){kmin=k2min;}
if(k2min<=k4min && modeP4==1 && modeP2==1){kmin=k2min;}

if(k4min<=k0min && modeP0==1 && modeP4==1){kmin=k4min;}
if(k4min<=k2min && modeP4==1 && modeP2==1){kmin=k4min;}

if(modeP0==1 && modeP2==0 && modeP4==0){kmax=k0max;kmin=k0min;}
if(modeP2==1 && modeP0==0 && modeP4==0){kmax=k2max;kmin=k2min;}
if(modeP4==1 && modeP2==0 && modeP0==0){kmax=k4max;kmin=k4min;}

NeffP=0;
if(strcmp(spacing_data,"linear") == 0){

if(kmax==k0max && kmin==k0min){NeffP=(NeffP0+25-2+(int)(k0min/(k0max-k0min)*(NeffP0-1.)))*factor_for_sampling;NeffP_max=NeffP0;}
if(kmax==k2max && kmin==k2min){NeffP=(NeffP2+25-2+(int)(k2min/(k2max-k2min)*(NeffP2-1.)))*factor_for_sampling;NeffP_max=NeffP2;}
if(kmax==k4max && kmin==k4min){NeffP=(NeffP4+25-2+(int)(k4min/(k4max-k4min)*(NeffP4-1.)))*factor_for_sampling;NeffP_max=NeffP4;}
}
if(strcmp(spacing_data,"log") == 0){
if(kmax==k0max && kmin==k0min){NeffP=(NeffP0+25-2+(int)(log(k0min)/(log(k0max)-log(k0min))*(NeffP0-1.)))*factor_for_sampling;NeffP_max=NeffP0;}
if(kmax==k2max && kmin==k2min){NeffP=(NeffP2+25-2+(int)(log(k2min)/(log(k2max)-log(k2min))*(NeffP2-1.)))*factor_for_sampling;NeffP_max=NeffP2;}
if(kmax==k4max && kmin==k4min){NeffP=(NeffP4+25-2+(int)(log(k4min)/(log(k4max)-log(k4min))*(NeffP4-1.)))*factor_for_sampling;NeffP_max=NeffP4;}
}

if(strcmp(spacing_data,"log10") == 0){

if(kmax==k0max && kmin==k0min){NeffP=(NeffP0+25-2+(int)(log10(k0min)/(log10(k0max)-log10(k0min))*(NeffP0-1.)))*factor_for_sampling;NeffP_max=NeffP0;}
if(kmax==k2max && kmin==k2min){NeffP=(NeffP2+25-2+(int)(log10(k2min)/(log10(k2max)-log10(k2min))*(NeffP2-1.)))*factor_for_sampling;NeffP_max=NeffP2;}
if(kmax==k4max && kmin==k4min){NeffP=(NeffP4+25-2+(int)(log10(k4min)/(log10(k4max)-log10(k4min))*(NeffP4-1.)))*factor_for_sampling;NeffP_max=NeffP4;}

}

if(strcmp(spacing_data,"irregular") == 0){

NeffP=N_Plin;

}

if(NeffP==0 && strcmp(path_to_mask1, "none") != 0 ){printf("Warning, crossing intervals among P0 (%lf<k<%lf), P2 (%lf<k<%lf), P4 (%lf<k<%lf) make impossible to determine NeffP within mask application. Fix this\n",k0min,k0max,k2min,k2max,k4min,k4max);exit(0);}

NeffP0=NeffP0*modeP0;
NeffP2=NeffP2*modeP2;
NeffP4=NeffP4*modeP4;


b1=parameters1[4];
b2=parameters1[5];
Anoise=parameters1[6];
sigmaP=parameters1[9];
sigma8_scaling=pow(parameters1[3]/Theory[0][41],2);
f=parameters1[2];
apara=parameters1[0];
aperp=parameters1[1];
bs2=parameters1[7];//-4./7.*(b1-1)
b3nl=parameters1[8];//32./315.*(b1-1)

        f_params *function_parameters;

        function_parameters = (f_params *) malloc(sizeof(f_params));

(*function_parameters).type_fog=fogmodel_ps;
(*function_parameters).type_ptmodel=ptmodel;
(*function_parameters).type_rsdmodel=rsdmodel_ps;

(*function_parameters).N=N_Plin;
(*function_parameters).sigmaP=sigmaP;
(*function_parameters).b1=b1;
(*function_parameters).b2=b2;
(*function_parameters).bs2=bs2;
(*function_parameters).b3nl=b3nl;
(*function_parameters).A=Anoise;
(*function_parameters).Pnoise=Pnoise;
(*function_parameters).f=f;
(*function_parameters).sigma8=sigma8_scaling;
(*function_parameters).a_parallel=apara;
(*function_parameters).a_perpendicular=aperp;
(*function_parameters).theory=Theory;
(*function_parameters).spacing=spacing_theory;
if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP0;}

//for(i=0;i<NeffP0;i++)printf("data: %d, %lf\n",i,k0[i]);
//exit(0);
if(modeP0==1)
{
    for(j=0;j<Nelements;j++)
    {
        if(strcmp(path_to_mask1, "none") != 0){

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"irregular") == 0){(*function_parameters).kinput=Theory[j][0];k_theo[j]=Theory[j][0];}


}
        else{k_theo0[j]=k0[j];(*function_parameters).kinput=k_theo0[j];}

         (*function_parameters).mode=0;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);
         P_theo0[j]=ptheo;
    }
}


if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP2;}

if(modeP2==1)
{
    for(j=0;j<Nelements;j++)
    {
        if(strcmp(path_to_mask1, "none") != 0){

//k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"irregular") == 0){(*function_parameters).kinput=Theory[j][0];k_theo[j]=Theory[j][0];}


}
        else{k_theo2[j]=k2[j];(*function_parameters).kinput=k_theo2[j];}

         (*function_parameters).mode=2;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo2,&junk);
         P_theo2[j]=ptheo2;
    }
}

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP4;}

if(modeP4==1)
{
    for(j=0;j<Nelements;j++)
    {
        if(strcmp(path_to_mask1, "none") != 0){

//k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

if(strcmp(spacing_data,"irregular") == 0){(*function_parameters).kinput=Theory[j][0];k_theo[j]=Theory[j][0];}

}

        else{k_theo4[j]=k4[j];(*function_parameters).kinput=k_theo4[j];}

         (*function_parameters).mode=4;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo4,&junk);
         P_theo4[j]=ptheo4;
    }
}



free(function_parameters);
//for(i=0;i<Nelements;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i],P_theo2[i]);}
//exit(0);

//aply mask
if(strcmp(path_to_mask1, "none") != 0)
{
//apply_mask(type_of_analysis, modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,N_Plin, pos, W0, W2, W4, W6,W8,Nmask, plan1, plan2, k_theo[0], k_theo[N_Plin-1], k0min, k0max, k2min, k2max, k4min, k4max);

//apply_mask(type_of_analysis, modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,NeffP, pos, W0, W2, W4, W6,W8,Nmask,spacing_data, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max);

apply_mask(type_of_analysis, modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,NeffP, pos, W0, W2, W4, W6,W8,Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,0,0);

}

//for(i=0;i<Nelements;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i],P_theo2[i]);}


//exit(0);
}

double chi2_bao_rsd(char *type_BAO_fit, char *type_of_analysis,char *fit_BAO,char *fit_RSD, double *parameters2_bao,double  *parameters2_rsd, double *k_Plin, double *Plin,int N_Plin,double *k_Olin,double *Olin,int N_Olin,double **Theory,int Ntheory, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0bao, double *k0rsd, double *P0bao, double *P0rsd, double Pnoise,int NeffP0bao,int NeffP0rsd, double *k2bao, double *k2rsd, double *P2bao, double *P2rsd,int NeffP2bao,int NeffP2rsd, double *k4bao, double *k4rsd, double *P4bao, double *P4rsd, int NeffP4bao, int NeffP4rsd, double *k0baoSGC, double *k0rsdSGC, double *P0baoSGC, double *P0rsdSGC,double PnoiseSGC,int NeffP0baoSGC,int NeffP0rsdSGC, double *k2baoSGC,double *k2rsdSGC, double *P2baoSGC,double *P2rsdSGC,int NeffP2baoSGC,int NeffP2rsdSGC, double *k4baoSGC, double *k4rsdSGC, double *P4baoSGC,double *P4rsdSGC, int NeffP4baoSGC, int NeffP4rsdSGC,double *cov, double *covSGC, char *Sigma_def_type, char *Sigma_independent,  double ffactor,double *Sigma_type,  double *Sigma_nl_mean,  double *Sigma_nl_stddev, int Npolynomial, char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free, char *fog_bs, int Nchunks, fftw_plan plan1bao, fftw_plan plan2bao, fftw_plan plan1rsd, fftw_plan plan2rsd ,char *do_power_spectrum, char *do_bispectrum, int Nalphas, int Nsigmas_tot,int Nsigmas_free, double Sigma_smooth,int factor_sampling_mask_in,char *spacing_dataNGC_bao,char *spacing_dataNGC_rsd,char *spacing_dataSGC_bao,char *spacing_dataSGC_rsd, char *spacing_theory_bao, char *spacing_theory_rsd)
{
//double chi2_rsd,chi2_bao;
//double chi2_rsdSGC,chi2_baoSGC;
double ch2,ch2SGC,prior_chi2;
double *difference;
double *parameters1;
int i,j,i1;
double ptheo,pobs;
double *k_theobao,*k_theo0bao,*k_theo2bao,*k_theo4bao;
double *P_theo0bao,*P_theo2bao,*P_theo4bao;
double *k_theorsd,*k_theo0rsd,*k_theo2rsd,*k_theo4rsd;
double *P_theo0rsd,*P_theo2rsd,*P_theo4rsd;
int modeP0bao,modeP2bao,modeP4bao;
int modeP0rsd,modeP2rsd,modeP4rsd; 
int Ncov,Ncovbao,Ncovrsd;
int points;
int offsetbao;
int offsetbao_ini;
int offsetrsd;
int offsetrsd_ini;
double Sigmanl0,Sigmanl2,Sigmanl4;
int Nsigmas_for_param1;
int Nalphas_for_param1;
int Neffmax;
int factor_sampling_mask;
int interpolation_order,shiftN;
double w1,w2,w0;
int Ninterpol;
int dimension;

interpolation_order=1;

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}

factor_sampling_mask=factor_sampling_mask_in;

modeP0bao=0;
modeP2bao=0;
modeP4bao=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0bao=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2bao=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4bao=1;}

Ncovbao=NeffP0bao*modeP0bao+NeffP2bao*modeP2bao+NeffP4bao*modeP4bao;
points=Ncovbao;

modeP0rsd=0;
modeP2rsd=0;
modeP4rsd=0;
if(strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP0rsd=1;}
if(strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP2rsd=1;}
if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP4rsd=1;}
Ncovrsd=NeffP0rsd*modeP0rsd+NeffP2rsd*modeP2rsd+NeffP4rsd*modeP4rsd;
points=points+Ncovrsd;

Ncov=Ncovbao+Ncovrsd;

if(strcmp(path_to_mask1, "none") == 0 )//no mask
{
if(NeffP0bao*modeP0bao>0){k_theo0bao = (double*) calloc( NeffP0bao, sizeof(double));}
if(NeffP2bao*modeP2bao>0){k_theo2bao = (double*) calloc( NeffP2bao, sizeof(double));}
if(NeffP4bao*modeP4bao>0){k_theo4bao = (double*) calloc( NeffP4bao, sizeof(double));}
if(NeffP0bao*modeP0bao>0){P_theo0bao = (double*) calloc( NeffP0bao, sizeof(double));}
if(NeffP2bao*modeP2bao>0){P_theo2bao = (double*) calloc( NeffP2bao, sizeof(double));}
if(NeffP4bao*modeP4bao>0){P_theo4bao = (double*) calloc( NeffP4bao, sizeof(double));}
}
else
{

if( strcmp(spacing_dataNGC_bao,"linear") == 0  ){

if(NeffP0bao*modeP0bao>=NeffP2bao*modeP2bao){Neffmax=NeffP0bao-2+25+(int)(k0bao[0]/(k0bao[NeffP0bao-1]-k0bao[0])*(NeffP0bao-1));}
if(NeffP0bao*modeP0bao>=NeffP4bao*modeP4bao){Neffmax=NeffP0bao-2+25+(int)(k0bao[0]/(k0bao[NeffP0bao-1]-k0bao[0])*(NeffP0bao-1));}

if(NeffP2bao*modeP2bao>=NeffP0bao*modeP0bao){Neffmax=NeffP2bao-2+25+(int)(k2bao[0]/(k2bao[NeffP2bao-1]-k2bao[0])*(NeffP2bao-1));}
if(NeffP2bao*modeP2bao>=NeffP4bao*modeP4bao){Neffmax=NeffP2bao-2+25+(int)(k2bao[0]/(k2bao[NeffP2bao-1]-k2bao[0])*(NeffP2bao-1));}

if(NeffP4bao*modeP4bao>=NeffP0bao*modeP0bao){Neffmax=NeffP4bao-2+25+(int)(k4bao[0]/(k4bao[NeffP4bao-1]-k4bao[0])*(NeffP4bao-1));}
if(NeffP4bao*modeP4bao>=NeffP2bao*modeP2bao){Neffmax=NeffP4bao-2+25+(int)(k4bao[0]/(k4bao[NeffP4bao-1]-k4bao[0])*(NeffP4bao-1));}
}

if( strcmp(spacing_dataNGC_bao,"log") == 0  ){

if(NeffP0bao*modeP0bao>=NeffP2bao*modeP2bao){Neffmax=NeffP0bao-2+25+(int)(log(k0bao[0])/(log(k0bao[NeffP0bao-1])-log(k0bao[0]))*(NeffP0bao-1));}
if(NeffP0bao*modeP0bao>=NeffP4bao*modeP4bao){Neffmax=NeffP0bao-2+25+(int)(log(k0bao[0])/(log(k0bao[NeffP0bao-1])-log(k0bao[0]))*(NeffP0bao-1));}

if(NeffP2bao*modeP2bao>=NeffP0bao*modeP0bao){Neffmax=NeffP2bao-2+25+(int)(log(k2bao[0])/(log(k2bao[NeffP2bao-1])-log(k2bao[0]))*(NeffP2bao-1));}
if(NeffP2bao*modeP2bao>=NeffP4bao*modeP4bao){Neffmax=NeffP2bao-2+25+(int)(log(k2bao[0])/(log(k2bao[NeffP2bao-1])-log(k2bao[0]))*(NeffP2bao-1));}

if(NeffP4bao*modeP4bao>=NeffP0bao*modeP0bao){Neffmax=NeffP4bao-2+25+(int)(log(k4bao[0])/(log(k4bao[NeffP4bao-1])-log(k4bao[0]))*(NeffP4bao-1));}
if(NeffP4bao*modeP4bao>=NeffP2bao*modeP2bao){Neffmax=NeffP4bao-2+25+(int)(log(k4bao[0])/(log(k4bao[NeffP4bao-1])-log(k4bao[0]))*(NeffP4bao-1));}

}

if( strcmp(spacing_dataNGC_bao,"log10") == 0  ){

if(NeffP0bao*modeP0bao>=NeffP2bao*modeP2bao){Neffmax=NeffP0bao-2+25+(int)(log10(k0bao[0])/(log10(k0bao[NeffP0bao-1])-log10(k0bao[0]))*(NeffP0bao-1));}
if(NeffP0bao*modeP0bao>=NeffP4bao*modeP4bao){Neffmax=NeffP0bao-2+25+(int)(log10(k0bao[0])/(log10(k0bao[NeffP0bao-1])-log10(k0bao[0]))*(NeffP0bao-1));}

if(NeffP2bao*modeP2bao>=NeffP0bao*modeP0bao){Neffmax=NeffP2bao-2+25+(int)(log10(k2bao[0])/(log10(k2bao[NeffP2bao-1])-log10(k2bao[0]))*(NeffP2bao-1));}
if(NeffP2bao*modeP2bao>=NeffP4bao*modeP4bao){Neffmax=NeffP2bao-2+25+(int)(log10(k2bao[0])/(log10(k2bao[NeffP2bao-1])-log10(k2bao[0]))*(NeffP2bao-1));}

if(NeffP4bao*modeP4bao>=NeffP0bao*modeP0bao){Neffmax=NeffP4bao-2+25+(int)(log10(k4bao[0])/(log10(k4bao[NeffP4bao-1])-log10(k4bao[0]))*(NeffP4bao-1));}
if(NeffP4bao*modeP4bao>=NeffP2bao*modeP2bao){Neffmax=NeffP4bao-2+25+(int)(log10(k4bao[0])/(log10(k4bao[NeffP4bao-1])-log10(k4bao[0]))*(NeffP4bao-1));}

}
if( strcmp(spacing_dataNGC_bao,"irregular") == 0  ){
Neffmax=N_Plin;
factor_sampling_mask=1;
}


k_theobao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0bao==1){P_theo0bao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2bao==1){P_theo2bao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4bao==1){P_theo4bao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
}

difference = (double*) calloc( Ncov, sizeof(double));


if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){

if(modeP0bao+modeP2bao+modeP4bao==1){Nalphas_for_param1=1;}
else{Nalphas_for_param1=2;}

Nsigmas_for_param1=modeP0bao+modeP2bao+modeP4bao;

parameters1 =  (double*) calloc( (modeP0bao+modeP2bao+modeP4bao)*(Npolynomial+1)+Nalphas_for_param1+Nsigmas_for_param1, sizeof(double));

offsetbao=Nalphas+Nsigmas_tot;

}

if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){

Nalphas_for_param1=2;
Nsigmas_for_param1=2;
parameters1 =  (double*) calloc( (modeP0bao+modeP2bao+modeP4bao)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1, sizeof(double));

offsetbao=Nalphas+Nsigmas_tot+1;

}

if(strcmp(Sigma_def_type, "effective") == 0)
{
if(modeP0bao+modeP2bao+modeP4bao==1){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];offsetbao_ini=2;}

if(modeP0bao+modeP2bao+modeP4bao==2 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[2];offsetbao_ini=4;}
if(modeP0bao+modeP2bao+modeP4bao==3 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[2];parameters1[4]=parameters2_bao[2];offsetbao_ini=5;}
if(modeP0bao+modeP2bao+modeP4bao==2 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[3];offsetbao_ini=4;}
if(modeP0bao+modeP2bao+modeP4bao==3 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[3];parameters1[4]=parameters2_bao[4];offsetbao_ini=5;}


}
else//para-perp
{

if(modeP0bao+modeP2bao+modeP4bao==1){

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[1],2./6.)*pow(parameters2_bao[2],4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[1],6./10.)*pow(parameters2_bao[2],4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[1],10./14.)*pow(parameters2_bao[2],4./14.);}
}
else
{
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[1],2./6.)*pow(parameters2_bao[1]/(1.+ffactor),4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[1],6./10.)*pow(parameters2_bao[1]/(1.+ffactor),4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[1],10./14.)*pow(parameters2_bao[1]/(1.+ffactor),4./14.);}
}

     parameters1[0]=parameters2_bao[0];
if(modeP0bao==1){parameters1[1]=Sigmanl0;}
if(modeP2bao==1){parameters1[1]=Sigmanl2;}
if(modeP4bao==1){parameters1[1]=Sigmanl4;}
offsetbao_ini=2;

}else{

     parameters1[0]=parameters2_bao[0];
     parameters1[1]=parameters2_bao[1];

if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){

if(strcmp(Sigma_independent, "yes") == 0)
{ 
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[2],2./6.)*pow(parameters2_bao[3],4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[2],6./10.)*pow(parameters2_bao[3],4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[2],10./14.)*pow(parameters2_bao[3],4./14.);}
}
else
{
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[2],2./6.)*pow(parameters2_bao[2]/(1.+ffactor),4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[2],6./10.)*pow(parameters2_bao[2]/(1.+ffactor),4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[2],10./14.)*pow(parameters2_bao[2]/(1.+ffactor),4./14.);}
}


if(modeP0bao==0){parameters1[2]=Sigmanl2;parameters1[3]=Sigmanl4;offsetbao_ini=4;}
if(modeP2bao==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl4;offsetbao_ini=4;}
if(modeP4bao==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;offsetbao_ini=4;}
if(modeP0bao+modeP2bao+modeP4bao==3){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;parameters1[4]=Sigmanl4;offsetbao_ini=5;}

}
if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){

parameters1[2]=parameters2_bao[2];
if(strcmp(Sigma_independent, "yes") == 0)
{
parameters1[3]=parameters2_bao[3];
parameters1[4]=parameters2_bao[4];
}
else
{
parameters1[3]=parameters2_bao[2]/(1+ffactor);
parameters1[4]=parameters2_bao[3];
}
offsetbao_ini=5;

}
}
}


if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){for(i=offsetbao_ini;i<(Npolynomial+1)*(modeP0bao+modeP2bao+modeP4bao)+offsetbao_ini;i++){parameters1[i]=parameters2_bao[i-offsetbao_ini+offsetbao];}}
if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){for(i=offsetbao_ini;i<1+(Npolynomial)*(modeP0bao+modeP2bao+modeP4bao)+offsetbao_ini;i++){parameters1[i]=parameters2_bao[i-offsetbao_ini+offsetbao];}}



if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){


do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0bao,NeffP2bao,NeffP4bao,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin,pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0bao,k2bao,k4bao, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0bao[0],k0bao[NeffP0bao-1],k2bao[0],k2bao[NeffP2bao-1], k4bao[0],k4bao[NeffP4bao-1], 1,spacing_dataNGC_bao,spacing_theory_bao,Sigma_smooth);
}

if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){

//for(i=0;i<(modeP0bao+modeP2bao+modeP4bao)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1;i++){printf("BAO-NGC: %d %lf\n",i,parameters1[i]);}

do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0bao,NeffP2bao,NeffP4bao,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, Sigma_smooth, pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0bao,k2bao,k4bao, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0bao[0],k0bao[NeffP0bao-1],k2bao[0],k2bao[NeffP2bao-1], k4bao[0],k4bao[NeffP4bao-1], 1,spacing_dataNGC_bao,spacing_theory_bao);
//if(strcmp(path_to_mask1, "none") != 0 ){printf("%lf %lf %lf\n",k_theobao[300],P_theo0bao[300],P_theo2bao[300]);}
//if(strcmp(path_to_mask1, "none") == 0 ){printf("%lf %lf\n",P_theo0bao[10],P_theo2bao[10]);}

}

free(parameters1);

if(strcmp(path_to_mask1, "none") == 0 )
{
i=-1;

        if(modeP0bao==1){
        for(j=0;j<NeffP0bao;j++){
        ptheo=P_theo0bao[j];
        pobs=P0bao[j];
        i++;
        difference[i]=ptheo-pobs;//printf("P0 %lf %lf %lf %lf\n",k0[j],k_theo0[j],pobs,ptheo);
        }
        }
       if(modeP2bao==1){
        for(j=0;j<NeffP2bao;j++){
        ptheo=P_theo2bao[j];
        pobs=P2bao[j];
        i++;
        difference[i]=ptheo-pobs;//printf("P2 %lf %lf %lf %lf\n",k2[j],k_theo2[j],pobs,ptheo);
        }
        }

        if(modeP4bao==1){
        for(j=0;j<NeffP4bao;j++){
        ptheo=P_theo4bao[j];
        pobs=P4bao[j];
        i++;
        difference[i]=ptheo-pobs;//printf("P4 %lf %lf %lf %lf\n",k4[j],k_theo4[j],pobs,ptheo);
        }
        }
}
else{

j=0;
    for(i=0;i<Ncovbao;i++)
    {

        if(modeP4bao==1 && modeP0bao==0 && modeP2bao==0){
Ninterpol=determine_N_singlearray(k_theobao,k4bao[j],Neffmax*factor_sampling_mask,spacing_dataNGC_bao);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4bao[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theobao,k4bao[j],Ninterpol,spacing_dataNGC_bao);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theobao,k4bao[j],Ninterpol,spacing_dataNGC_bao);
w1=determine_w1_2ndorder_singlearray(k_theobao,k4bao[j],Ninterpol,spacing_dataNGC_bao);
w2=determine_w2_2ndorder_singlearray(k_theobao,k4bao[j],Ninterpol,spacing_dataNGC_bao);
}
ptheo=P_interpol_fast(k4bao[j],P_theo4bao,Neffmax*factor_sampling_mask,spacing_dataNGC_bao,interpolation_order,Ninterpol,w0,w1,w2);
}



        pobs=P4bao[j];
        j++;
        if(j==NeffP4bao){j=0;modeP4bao=0;}
        }

        if(modeP2bao==1 && modeP0bao==0){
Ninterpol=determine_N_singlearray(k_theobao,k2bao[j],Neffmax*factor_sampling_mask,spacing_dataNGC_bao);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2bao[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theobao,k2bao[j],Ninterpol,spacing_dataNGC_bao);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theobao,k2bao[j],Ninterpol,spacing_dataNGC_bao);
w1=determine_w1_2ndorder_singlearray(k_theobao,k2bao[j],Ninterpol,spacing_dataNGC_bao);
w2=determine_w2_2ndorder_singlearray(k_theobao,k2bao[j],Ninterpol,spacing_dataNGC_bao);
}
ptheo=P_interpol_fast(k2bao[j],P_theo2bao,Neffmax*factor_sampling_mask,spacing_dataNGC_bao,interpolation_order,Ninterpol,w0,w1,w2);
}



        pobs=P2bao[j];
        j++;
        if(j==NeffP2bao){j=0;modeP2bao=0;}
        }

        if(modeP0bao==1){
Ninterpol=determine_N_singlearray(k_theobao,k0bao[j],Neffmax*factor_sampling_mask,spacing_dataNGC_bao);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0bao[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theobao,k0bao[j],Ninterpol,spacing_dataNGC_bao);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theobao,k0bao[j],Ninterpol,spacing_dataNGC_bao);
w1=determine_w1_2ndorder_singlearray(k_theobao,k0bao[j],Ninterpol,spacing_dataNGC_bao);
w2=determine_w2_2ndorder_singlearray(k_theobao,k0bao[j],Ninterpol,spacing_dataNGC_bao);
}
ptheo=P_interpol_fast(k0bao[j],P_theo0bao,Neffmax*factor_sampling_mask,spacing_dataNGC_bao,interpolation_order,Ninterpol,w0,w1,w2);
}


        pobs=P0bao[j];
        j++;
        if(j==NeffP0bao){j=0;modeP0bao=0;}
        }

        difference[i]=ptheo-pobs;
   }

}

//RSD aqui
factor_sampling_mask=factor_sampling_mask_in;

offsetrsd=3;
if(strcmp(local_b2s2, "no") == 0){offsetrsd++;}//b2s2
if(strcmp(local_b3nl, "no") == 0){offsetrsd++;}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){offsetrsd++;}//fog_ps
if(strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 ){offsetrsd++;}//fog_bs

dimension=11;//This is a fixed quantity.


if(strcmp(path_to_mask1, "none") == 0)
{
if(NeffP0rsd*modeP0rsd>0){k_theo0rsd = (double*) calloc( NeffP0rsd, sizeof(double));}
if(NeffP2rsd*modeP2rsd>0){k_theo2rsd = (double*) calloc( NeffP2rsd, sizeof(double));}
if(NeffP4rsd*modeP4rsd>0){k_theo4rsd = (double*) calloc( NeffP4rsd, sizeof(double));}
if(NeffP0rsd*modeP0rsd>0){P_theo0rsd = (double*) calloc( NeffP0rsd, sizeof(double));}
if(NeffP2rsd*modeP2rsd>0){P_theo2rsd = (double*) calloc( NeffP2rsd, sizeof(double));}
if(NeffP4rsd*modeP4rsd>0){P_theo4rsd = (double*) calloc( NeffP4rsd, sizeof(double));}
}
else{
if( strcmp(spacing_dataNGC_rsd,"linear") == 0  ){

if(NeffP0rsd*modeP0rsd>=NeffP2rsd*modeP2rsd){Neffmax=NeffP0rsd-2+25+(int)(k0rsd[0]/(k0rsd[NeffP0rsd-1]-k0rsd[0])*(NeffP0rsd-1));}
if(NeffP0rsd*modeP0rsd>=NeffP4rsd*modeP4rsd){Neffmax=NeffP0rsd-2+25+(int)(k0rsd[0]/(k0rsd[NeffP0rsd-1]-k0rsd[0])*(NeffP0rsd-1));}

if(NeffP2rsd*modeP2rsd>=NeffP0rsd*modeP0rsd){Neffmax=NeffP2rsd-2+25+(int)(k2rsd[0]/(k2rsd[NeffP2rsd-1]-k2rsd[0])*(NeffP2rsd-1));}
if(NeffP2rsd*modeP2rsd>=NeffP4rsd*modeP4rsd){Neffmax=NeffP2rsd-2+25+(int)(k2rsd[0]/(k2rsd[NeffP2rsd-1]-k2rsd[0])*(NeffP2rsd-1));}

if(NeffP4rsd*modeP4rsd>=NeffP0rsd*modeP0rsd){Neffmax=NeffP4rsd-2+25+(int)(k4rsd[0]/(k4rsd[NeffP4rsd-1]-k4rsd[0])*(NeffP4rsd-1));}
if(NeffP4rsd*modeP4rsd>=NeffP2rsd*modeP2rsd){Neffmax=NeffP4rsd-2+25+(int)(k4rsd[0]/(k4rsd[NeffP4rsd-1]-k4rsd[0])*(NeffP4rsd-1));}

}

if( strcmp(spacing_dataNGC_rsd,"log") == 0  ){

if(NeffP0rsd*modeP0rsd>=NeffP2rsd*modeP2rsd){Neffmax=NeffP0rsd-2+25+(int)(log(k0rsd[0])/(log(k0rsd[NeffP0rsd-1])-log(k0rsd[0]))*(NeffP0rsd-1));}
if(NeffP0rsd*modeP0rsd>=NeffP4rsd*modeP4rsd){Neffmax=NeffP0rsd-2+25+(int)(log(k0rsd[0])/(log(k0rsd[NeffP0rsd-1])-log(k0rsd[0]))*(NeffP0rsd-1));}

if(NeffP2rsd*modeP2rsd>=NeffP0rsd*modeP0rsd){Neffmax=NeffP2rsd-2+25+(int)(log(k2rsd[0])/(log(k2rsd[NeffP2rsd-1])-log(k2rsd[0]))*(NeffP2rsd-1));}
if(NeffP2rsd*modeP2rsd>=NeffP4rsd*modeP4rsd){Neffmax=NeffP2rsd-2+25+(int)(log(k2rsd[0])/(log(k2rsd[NeffP2rsd-1])-log(k2rsd[0]))*(NeffP2rsd-1));}

if(NeffP4rsd*modeP4rsd>=NeffP0rsd*modeP0rsd){Neffmax=NeffP4rsd-2+25+(int)(log(k4rsd[0])/(log(k4rsd[NeffP4rsd-1])-log(k4rsd[0]))*(NeffP4rsd-1));}
if(NeffP4rsd*modeP4rsd>=NeffP2rsd*modeP2rsd){Neffmax=NeffP4rsd-2+25+(int)(log(k4rsd[0])/(log(k4rsd[NeffP4rsd-1])-log(k4rsd[0]))*(NeffP4rsd-1));}

}

if( strcmp(spacing_dataNGC_rsd,"log10") == 0  ){

if(NeffP0rsd*modeP0rsd>=NeffP2rsd*modeP2rsd){Neffmax=NeffP0rsd-2+25+(int)(log10(k0rsd[0])/(log10(k0rsd[NeffP0rsd-1])-log10(k0rsd[0]))*(NeffP0rsd-1));}
if(NeffP0rsd*modeP0rsd>=NeffP4rsd*modeP4rsd){Neffmax=NeffP0rsd-2+25+(int)(log10(k0rsd[0])/(log10(k0rsd[NeffP0rsd-1])-log10(k0rsd[0]))*(NeffP0rsd-1));}

if(NeffP2rsd*modeP2rsd>=NeffP0rsd*modeP0rsd){Neffmax=NeffP2rsd-2+25+(int)(log10(k2rsd[0])/(log10(k2rsd[NeffP2rsd-1])-log10(k2rsd[0]))*(NeffP2rsd-1));}
if(NeffP2rsd*modeP2rsd>=NeffP4rsd*modeP4rsd){Neffmax=NeffP2rsd-2+25+(int)(log10(k2rsd[0])/(log10(k2rsd[NeffP2rsd-1])-log10(k2rsd[0]))*(NeffP2rsd-1));}

if(NeffP4rsd*modeP4rsd>=NeffP0rsd*modeP0rsd){Neffmax=NeffP4rsd-2+25+(int)(log10(k4rsd[0])/(log10(k4rsd[NeffP4rsd-1])-log10(k4rsd[0]))*(NeffP4rsd-1));}
if(NeffP4rsd*modeP4rsd>=NeffP2rsd*modeP2rsd){Neffmax=NeffP4rsd-2+25+(int)(log10(k4rsd[0])/(log10(k4rsd[NeffP4rsd-1])-log10(k4rsd[0]))*(NeffP4rsd-1));}

}
if( strcmp(spacing_dataNGC_rsd,"irregular") == 0  ){
Neffmax=Ntheory;
factor_sampling_mask=1;
}


k_theorsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0rsd==1){P_theo0rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2rsd==1){P_theo2rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4rsd==1){P_theo4rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}


}

//difference = (double*) calloc( Ncov, sizeof(double));

parameters1 = (double*) calloc( dimension, sizeof(double));

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<3){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}
}

if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==3 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}
if(i>=4 && i<=6){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==7 && strcmp(local_b2s2, "yes") == 0){
parameters1[i]=-4./7.*(parameters1[i-3]-1);
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==8 && strcmp(local_b3nl, "yes") == 0){
parameters1[i]=32./315.*(parameters1[i-4]-1);
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==9 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=0;
}

if(i==9 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){
parameters1[i]=parameters2_rsd[i1];
i1++;
}
if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){
parameters1[i]=parameters1[i-1];
}

if(i==10 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){
parameters1[i]=0;
}


}

//for(i=0;i<dimension;i++){printf("RSD-NGC: %d %lf\n",i,parameters1[i]);}

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_RSD,modeP0rsd,modeP2rsd,modeP4rsd,Theory,Ntheory, k_theorsd,k_theo0rsd,k_theo2rsd,k_theo4rsd, P_theo0rsd, P_theo2rsd, P_theo4rsd,Pnoise, NeffP0rsd,NeffP2rsd,NeffP4rsd,factor_sampling_mask, parameters1,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC, k0rsd,k2rsd,k4rsd,path_to_mask1,plan1rsd, plan2rsd, k0rsd[0], k0rsd[NeffP0rsd-1], k2rsd[0], k2rsd[NeffP2rsd-1], k4rsd[0], k4rsd[NeffP4rsd-1],spacing_dataNGC_rsd,spacing_theory_rsd);
//if(strcmp(path_to_mask1, "none") != 0 ){printf("%lf %lf %lf %lf\n",k_theorsd[300],P_theo0rsd[300],P_theo2rsd[300],P_theo4rsd[300]);}
//if(strcmp(path_to_mask1, "none") == 0 ){printf("%lf %lf %lf\n",P_theo0rsd[10],P_theo2rsd[10],P_theo4rsd[10]);}

free(parameters1);

if(strcmp(path_to_mask1, "none") == 0 )
{
i=-1;

   if(modeP0rsd==1){
        for(j=0;j<NeffP0rsd;j++){
        ptheo=P_theo0rsd[j]-Pnoise;
        pobs=P0rsd[j];
        i++;
        difference[i+Ncovbao]=ptheo-pobs;//printf("%lf %lf %lf\n",k_theo0[j],P0[j],P_theo0[j]);
        }
        }

        if(modeP2rsd==1){
        for(j=0;j<NeffP2rsd;j++){
        ptheo=P_theo2rsd[j];
        pobs=P2rsd[j];
        i++;
        difference[i+Ncovbao]=ptheo-pobs;
        }
        }

        if(modeP4rsd==1){
        for(j=0;j<NeffP4rsd;j++){
        ptheo=P_theo4rsd[j];
        pobs=P4rsd[j];
        i++;
        difference[i+Ncovbao]=ptheo-pobs;
        }
        }


}
else{
j=0;
    for(i=0;i<Ncovrsd;i++)
    {

        if(modeP4rsd==1 && modeP0rsd==0 && modeP2rsd==0){
Ninterpol=determine_N_singlearray(k_theorsd,k4rsd[j],Neffmax*factor_sampling_mask,spacing_dataNGC_rsd);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4rsd[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theorsd,k4rsd[j],Ninterpol,spacing_dataNGC_rsd);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theorsd,k4rsd[j],Ninterpol,spacing_dataNGC_rsd);
w1=determine_w1_2ndorder_singlearray(k_theorsd,k4rsd[j],Ninterpol,spacing_dataNGC_rsd);
w2=determine_w2_2ndorder_singlearray(k_theorsd,k4rsd[j],Ninterpol,spacing_dataNGC_rsd);
}

ptheo=P_interpol_fast(k4rsd[j],P_theo4rsd,Neffmax*factor_sampling_mask,spacing_dataNGC_rsd,interpolation_order,Ninterpol,w0,w1,w2);

}


        pobs=P4rsd[j];
        j++;
        if(j==NeffP4rsd){j=0;modeP4rsd=0;}
        }

        if(modeP2rsd==1 && modeP0rsd==0){
Ninterpol=determine_N_singlearray(k_theorsd,k2rsd[j],Neffmax*factor_sampling_mask,spacing_dataNGC_rsd);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2rsd[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theorsd,k2rsd[j],Ninterpol,spacing_dataNGC_rsd);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theorsd,k2rsd[j],Ninterpol,spacing_dataNGC_rsd);
w1=determine_w1_2ndorder_singlearray(k_theorsd,k2rsd[j],Ninterpol,spacing_dataNGC_rsd);
w2=determine_w2_2ndorder_singlearray(k_theorsd,k2rsd[j],Ninterpol,spacing_dataNGC_rsd);
}
ptheo=P_interpol_fast(k2rsd[j],P_theo2rsd,Neffmax*factor_sampling_mask,spacing_dataNGC_rsd,interpolation_order,Ninterpol,w0,w1,w2);
}
        pobs=P2rsd[j];
        j++;
        if(j==NeffP2rsd){j=0;modeP2rsd=0;}
        }

        if(modeP0rsd==1){
Ninterpol=determine_N_singlearray(k_theorsd,k0rsd[j],Neffmax*factor_sampling_mask,spacing_dataNGC_rsd);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0rsd[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theorsd,k0rsd[j],Ninterpol,spacing_dataNGC_rsd);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theorsd,k0rsd[j],Ninterpol,spacing_dataNGC_rsd);
w1=determine_w1_2ndorder_singlearray(k_theorsd,k0rsd[j],Ninterpol,spacing_dataNGC_rsd);
w2=determine_w2_2ndorder_singlearray(k_theorsd,k0rsd[j],Ninterpol,spacing_dataNGC_rsd);
}
ptheo=P_interpol_fast(k0rsd[j],P_theo0rsd,Neffmax*factor_sampling_mask,spacing_dataNGC_rsd,interpolation_order,Ninterpol,w0,w1,w2)-Pnoise;
}


        pobs=P0rsd[j];
        j++;
        if(j==NeffP0rsd){j=0;modeP0rsd=0;}
        }

        difference[i+Ncovbao]=ptheo-pobs;
   }

}




ch2=0;
//chi2_rsd=0;
//chi2_bao=0;
for(i=0;i<Ncov;i++)
{

   for(j=0;j<Ncov;j++)
      {

              ch2=ch2+difference[i]*1./cov[i+Ncov*j]*difference[j];

//if(i<Ncovbao && j<Ncovbao){chi2_bao=chi2_bao+difference[i]*1./cov[i+Ncov*j]*difference[j];}
//if(i>=Ncovbao && j>=Ncovbao){chi2_rsd=chi2_rsd+difference[i]*1./cov[i+Ncov*j]*difference[j];}

      }

}

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0bao=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2bao=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4bao=1;}

if(strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP0rsd=1;}
if(strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP2rsd=1;}
if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP4rsd=1;}

if(strcmp(path_to_mask1, "none") != 0 ){free(k_theobao);}
if(strcmp(path_to_mask1, "none") != 0 ){free(k_theorsd);}

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo0bao);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo0bao);}
}

if(strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0)
{
free(P_theo0rsd);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo0rsd);}
}

if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo2bao);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo2bao);}
}

if(strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0)
{
free(P_theo2rsd);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo2rsd);}
}


if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo4bao);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo4bao);}
}

if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0)
{
free(P_theo4rsd);
if(strcmp(path_to_mask1, "none") == 0 ){free(k_theo4rsd);}
}


free(difference);


if(Nchunks==2)
{
factor_sampling_mask=factor_sampling_mask_in;

modeP0bao=0;
modeP2bao=0;
modeP4bao=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0bao=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2bao=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4bao=1;}

Ncovbao=NeffP0baoSGC*modeP0bao+NeffP2baoSGC*modeP2bao+NeffP4baoSGC*modeP4bao;
points=points+Ncovbao;

modeP0rsd=0;
modeP2rsd=0;
modeP4rsd=0;
if(strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP0rsd=1;}
if(strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP2rsd=1;}
if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP4rsd=1;}

Ncovrsd=NeffP0rsdSGC*modeP0rsd+NeffP2rsdSGC*modeP2rsd+NeffP4rsdSGC*modeP4rsd;
points=points+Ncovrsd;

Ncov=Ncovbao+Ncovrsd;

if(strcmp(path_to_mask2, "none") == 0 )
{
if(NeffP0baoSGC*modeP0bao>0){k_theo0bao = (double*) calloc( NeffP0baoSGC, sizeof(double));}
if(NeffP2baoSGC*modeP2bao>0){k_theo2bao = (double*) calloc( NeffP2baoSGC, sizeof(double));}
if(NeffP4baoSGC*modeP4bao>0){k_theo4bao = (double*) calloc( NeffP4baoSGC, sizeof(double));}
if(NeffP0baoSGC*modeP0bao>0){P_theo0bao = (double*) calloc( NeffP0baoSGC, sizeof(double));}
if(NeffP2baoSGC*modeP2bao>0){P_theo2bao = (double*) calloc( NeffP2baoSGC, sizeof(double));}
if(NeffP4baoSGC*modeP4bao>0){P_theo4bao = (double*) calloc( NeffP4baoSGC, sizeof(double));}
}
else{
if( strcmp(spacing_dataSGC_bao,"linear") == 0  ){

if(NeffP0baoSGC*modeP0bao>=NeffP2baoSGC*modeP2bao){Neffmax=NeffP0baoSGC-2+25+(int)(k0baoSGC[0]/(k0baoSGC[NeffP0baoSGC-1]-k0baoSGC[0])*(NeffP0baoSGC-1));}
if(NeffP0baoSGC*modeP0bao>=NeffP4baoSGC*modeP4bao){Neffmax=NeffP0baoSGC-2+25+(int)(k0baoSGC[0]/(k0baoSGC[NeffP0baoSGC-1]-k0baoSGC[0])*(NeffP0baoSGC-1));}

if(NeffP2baoSGC*modeP2bao>=NeffP0baoSGC*modeP0bao){Neffmax=NeffP2baoSGC-2+25+(int)(k2baoSGC[0]/(k2baoSGC[NeffP2baoSGC-1]-k2baoSGC[0])*(NeffP2baoSGC-1));}
if(NeffP2baoSGC*modeP2bao>=NeffP4baoSGC*modeP4bao){Neffmax=NeffP2baoSGC-2+25+(int)(k2baoSGC[0]/(k2baoSGC[NeffP2baoSGC-1]-k2baoSGC[0])*(NeffP2baoSGC-1));}

if(NeffP4baoSGC*modeP4bao>=NeffP0baoSGC*modeP0bao){Neffmax=NeffP4baoSGC-2+25+(int)(k4baoSGC[0]/(k4baoSGC[NeffP4baoSGC-1]-k4baoSGC[0])*(NeffP4baoSGC-1));}
if(NeffP4baoSGC*modeP4bao>=NeffP2baoSGC*modeP2bao){Neffmax=NeffP4baoSGC-2+25+(int)(k4baoSGC[0]/(k4baoSGC[NeffP4baoSGC-1]-k4baoSGC[0])*(NeffP4baoSGC-1));}

}

if( strcmp(spacing_dataSGC_bao,"log") == 0  ){

if(NeffP0baoSGC*modeP0bao>=NeffP2baoSGC*modeP2bao){Neffmax=NeffP0baoSGC-2+25+(int)(log(k0baoSGC[0])/(log(k0baoSGC[NeffP0baoSGC-1])-log(k0baoSGC[0]))*(NeffP0baoSGC-1));}
if(NeffP0baoSGC*modeP0bao>=NeffP4baoSGC*modeP4bao){Neffmax=NeffP0baoSGC-2+25+(int)(log(k0baoSGC[0])/(log(k0baoSGC[NeffP0baoSGC-1])-log(k0baoSGC[0]))*(NeffP0baoSGC-1));}

if(NeffP2baoSGC*modeP2bao>=NeffP0baoSGC*modeP0bao){Neffmax=NeffP2baoSGC-2+25+(int)(log(k2baoSGC[0])/(log(k2baoSGC[NeffP2baoSGC-1])-log(k2baoSGC[0]))*(NeffP2baoSGC-1));}
if(NeffP2baoSGC*modeP2bao>=NeffP4baoSGC*modeP4bao){Neffmax=NeffP2baoSGC-2+25+(int)(log(k2baoSGC[0])/(log(k2baoSGC[NeffP2baoSGC-1])-log(k2baoSGC[0]))*(NeffP2baoSGC-1));}

if(NeffP4baoSGC*modeP4bao>=NeffP0baoSGC*modeP0bao){Neffmax=NeffP4baoSGC-2+25+(int)(log(k4baoSGC[0])/(log(k4baoSGC[NeffP4baoSGC-1])-log(k4baoSGC[0]))*(NeffP4baoSGC-1));}
if(NeffP4baoSGC*modeP4bao>=NeffP2baoSGC*modeP2bao){Neffmax=NeffP4baoSGC-2+25+(int)(log(k4baoSGC[0])/(log(k4baoSGC[NeffP4baoSGC-1])-log(k4baoSGC[0]))*(NeffP4baoSGC-1));}

}

if( strcmp(spacing_dataSGC_bao,"log10") == 0  ){

if(NeffP0baoSGC*modeP0bao>=NeffP2baoSGC*modeP2bao){Neffmax=NeffP0baoSGC-2+25+(int)(log10(k0baoSGC[0])/(log10(k0baoSGC[NeffP0baoSGC-1])-log10(k0baoSGC[0]))*(NeffP0baoSGC-1));}
if(NeffP0baoSGC*modeP0bao>=NeffP4baoSGC*modeP4bao){Neffmax=NeffP0baoSGC-2+25+(int)(log10(k0baoSGC[0])/(log10(k0baoSGC[NeffP0baoSGC-1])-log10(k0baoSGC[0]))*(NeffP0baoSGC-1));}

if(NeffP2baoSGC*modeP2bao>=NeffP0baoSGC*modeP0bao){Neffmax=NeffP2baoSGC-2+25+(int)(log10(k2baoSGC[0])/(log10(k2baoSGC[NeffP2baoSGC-1])-log10(k2baoSGC[0]))*(NeffP2baoSGC-1));}
if(NeffP2baoSGC*modeP2bao>=NeffP4baoSGC*modeP4bao){Neffmax=NeffP2baoSGC-2+25+(int)(log10(k2baoSGC[0])/(log10(k2baoSGC[NeffP2baoSGC-1])-log10(k2baoSGC[0]))*(NeffP2baoSGC-1));}

if(NeffP4baoSGC*modeP4bao>=NeffP0baoSGC*modeP0bao){Neffmax=NeffP4baoSGC-2+25+(int)(log10(k4baoSGC[0])/(log10(k4baoSGC[NeffP4baoSGC-1])-log10(k4baoSGC[0]))*(NeffP4baoSGC-1));}
if(NeffP4baoSGC*modeP4bao>=NeffP2baoSGC*modeP2bao){Neffmax=NeffP4baoSGC-2+25+(int)(log10(k4baoSGC[0])/(log10(k4baoSGC[NeffP4baoSGC-1])-log10(k4baoSGC[0]))*(NeffP4baoSGC-1));}

}
if( strcmp(spacing_dataSGC_bao,"irregular") == 0  ){
Neffmax=N_Plin;
factor_sampling_mask=1;
}


k_theobao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0bao==1){P_theo0bao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2bao==1){P_theo2bao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4bao==1){P_theo4bao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
}

difference = (double*) calloc( Ncov, sizeof(double));

if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){
if(modeP0bao+modeP2bao+modeP4bao==1){Nalphas_for_param1=1;}
else{Nalphas_for_param1=2;}
Nsigmas_for_param1=modeP0bao+modeP2bao+modeP4bao;
parameters1 =  (double*) calloc( (modeP0bao+modeP2bao+modeP4bao)*(Npolynomial+1)+Nalphas_for_param1+Nsigmas_for_param1, sizeof(double));
offsetbao=Nalphas+Nsigmas_tot;
}

if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){
Nalphas_for_param1=2;
Nsigmas_for_param1=2;;
parameters1 =  (double*) calloc( (modeP0bao+modeP2bao+modeP4bao)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1, sizeof(double));
offsetbao=Nalphas+Nsigmas_tot+1;
}


if(strcmp(Sigma_def_type, "effective") == 0)
{
if(modeP0bao+modeP2bao+modeP4bao==1){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];offsetbao_ini=2;}
if(modeP0bao+modeP2bao+modeP4bao==2 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[2];offsetbao_ini=4;}
if(modeP0bao+modeP2bao+modeP4bao==3 && strcmp(Sigma_independent, "no") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[2];parameters1[4]=parameters2_bao[2];offsetbao_ini=5;}
if(modeP0bao+modeP2bao+modeP4bao==2 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[3];offsetbao_ini=4;}
if(modeP0bao+modeP2bao+modeP4bao==3 && strcmp(Sigma_independent, "yes") == 0){parameters1[0]=parameters2_bao[0];parameters1[1]=parameters2_bao[1];parameters1[2]=parameters2_bao[2];parameters1[3]=parameters2_bao[3];parameters1[4]=parameters2_bao[4];offsetbao_ini=5;}


}
else
{


if(modeP0bao+modeP2bao+modeP4bao==1)
{

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[1],2./6.)*pow(parameters2_bao[2],4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[1],6./10.)*pow(parameters2_bao[2],4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[1],10./14.)*pow(parameters2_bao[2],4./14.);}
}
else
{
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[1],2./6.)*pow(parameters2_bao[1]/(1.+ffactor),4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[1],6./10.)*pow(parameters2_bao[1]/(1.+ffactor),4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[1],10./14.)*pow(parameters2_bao[1]/(1.+ffactor),4./14.);}
}

     parameters1[0]=parameters2_bao[0];
if(modeP0bao==1){parameters1[1]=Sigmanl0;}
if(modeP2bao==1){parameters1[1]=Sigmanl2;}
if(modeP4bao==1){parameters1[1]=Sigmanl4;}
offsetbao_ini=2;
}
else{
    parameters1[0]=parameters2_bao[0];
     parameters1[1]=parameters2_bao[1];

if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){

if(strcmp(Sigma_independent, "yes") == 0)
{
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[2],2./6.)*pow(parameters2_bao[3],4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[2],6./10.)*pow(parameters2_bao[3],4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[2],10./14.)*pow(parameters2_bao[3],4./14.);}
}else
{
if(modeP0bao==1){Sigmanl0=pow(parameters2_bao[2],2./6.)*pow(parameters2_bao[2]/(1.+ffactor),4./6.);}
if(modeP2bao==1){Sigmanl2=pow(parameters2_bao[2],6./10.)*pow(parameters2_bao[2]/(1.+ffactor),4./10.);}
if(modeP4bao==1){Sigmanl4=pow(parameters2_bao[2],10./14.)*pow(parameters2_bao[2]/(1.+ffactor),4./14.);}
}


if(modeP0bao==0){parameters1[2]=Sigmanl2;parameters1[3]=Sigmanl4;offsetbao_ini=4;}
if(modeP2bao==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl4;offsetbao_ini=4;}
if(modeP4bao==0){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;offsetbao_ini=4;}
if(modeP0bao+modeP2bao+modeP4bao==3){parameters1[2]=Sigmanl0;parameters1[3]=Sigmanl2;parameters1[4]=Sigmanl4;offsetbao_ini=5;}


}
if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){

parameters1[2]=parameters2_bao[2];
if(strcmp(Sigma_independent, "yes") == 0)
{
parameters1[3]=parameters2_bao[3];
parameters1[4]=parameters2_bao[4];
}
else
{
parameters1[3]=parameters2_bao[2]/(1+ffactor);
parameters1[4]=parameters2_bao[3];
}
offsetbao_ini=5;


}


}

}



if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){for(i=offsetbao_ini;i<(Npolynomial+1)*(modeP0bao+modeP2bao+modeP4bao)+offsetbao_ini;i++){parameters1[i]=parameters2_bao[i+(Npolynomial+1)*(modeP0bao+modeP2bao+modeP4bao)-offsetbao_ini+offsetbao];}}

if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){for(i=offsetbao_ini;i<(Npolynomial)*(modeP0bao+modeP2bao+modeP4bao)+1+offsetbao_ini;i++){parameters1[i]=parameters2_bao[i+1+(Npolynomial)*(modeP0bao+modeP2bao+modeP4bao)-offsetbao_ini+offsetbao];}}

if( strcmp(type_of_analysis, "FSBAOISO") == 0 ){
do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0baoSGC,NeffP2baoSGC,NeffP4baoSGC,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, spacing_maskSGC,path_to_mask2,k0baoSGC,k2baoSGC,k4baoSGC, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0baoSGC[0],k0baoSGC[NeffP0baoSGC-1],k2baoSGC[0],k2baoSGC[NeffP2baoSGC-1], k4baoSGC[0],k4baoSGC[NeffP4baoSGC-1], 1,spacing_dataSGC_bao,spacing_theory_bao,Sigma_smooth);
}

if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){

//for(i=0;i<(modeP0bao+modeP2bao+modeP4bao)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1;i++){printf("BAO-SGC: %d %lf\n",i,parameters1[i]);}


do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0baoSGC,NeffP2baoSGC,NeffP4baoSGC,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, Sigma_smooth,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, spacing_maskSGC, path_to_mask2,k0baoSGC,k2baoSGC,k4baoSGC, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0baoSGC[0],k0baoSGC[NeffP0baoSGC-1],k2baoSGC[0],k2baoSGC[NeffP2baoSGC-1], k4baoSGC[0],k4baoSGC[NeffP4baoSGC-1], 1,spacing_dataSGC_bao,spacing_theory_bao);
}

free(parameters1);

if(strcmp(path_to_mask2, "none") == 0 )
{
i=-1;
        if(modeP0bao==1){
        for(j=0;j<NeffP0baoSGC;j++){
        ptheo=P_theo0bao[j];
        pobs=P0baoSGC[j];
        i++;
        difference[i]=ptheo-pobs;
        }
        }

        if(modeP2bao==1){
        for(j=0;j<NeffP2baoSGC;j++){
        ptheo=P_theo2bao[j];
        pobs=P2baoSGC[j];
        i++;
        difference[i]=ptheo-pobs;
        }
        }

        if(modeP4bao==1){
        for(j=0;j<NeffP4baoSGC;j++){
        ptheo=P_theo4bao[j];
        pobs=P4baoSGC[j];
        i++;
        difference[i]=ptheo-pobs;
        }
        }
}
else{


j=0;
    for(i=0;i<Ncovbao;i++)
    {

        if(modeP4bao==1 && modeP0bao==0 && modeP2bao==0){
Ninterpol=determine_N_singlearray(k_theobao,k4baoSGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC_bao);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4baoSGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theobao,k4baoSGC[j],Ninterpol,spacing_dataSGC_bao);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theobao,k4baoSGC[j],Ninterpol,spacing_dataSGC_bao);
w1=determine_w1_2ndorder_singlearray(k_theobao,k4baoSGC[j],Ninterpol,spacing_dataSGC_bao);
w2=determine_w2_2ndorder_singlearray(k_theobao,k4baoSGC[j],Ninterpol,spacing_dataSGC_bao);
}
ptheo=P_interpol_fast(k4baoSGC[j],P_theo4bao,Neffmax*factor_sampling_mask,spacing_dataSGC_bao,interpolation_order,Ninterpol,w0,w1,w2);
}


        pobs=P4baoSGC[j];
        j++;
        if(j==NeffP4baoSGC){j=0;modeP4bao=0;}
        }

        if(modeP2bao==1 && modeP0bao==0){
Ninterpol=determine_N_singlearray(k_theobao,k2baoSGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC_bao);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2baoSGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theobao,k2baoSGC[j],Ninterpol,spacing_dataSGC_bao);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theobao,k2baoSGC[j],Ninterpol,spacing_dataSGC_bao);
w1=determine_w1_2ndorder_singlearray(k_theobao,k2baoSGC[j],Ninterpol,spacing_dataSGC_bao);
w2=determine_w2_2ndorder_singlearray(k_theobao,k2baoSGC[j],Ninterpol,spacing_dataSGC_bao);
}
ptheo=P_interpol_fast(k2baoSGC[j],P_theo2bao,Neffmax*factor_sampling_mask,spacing_dataSGC_bao,interpolation_order,Ninterpol,w0,w1,w2);
}

        pobs=P2baoSGC[j];//printf("%lf %lf %lf\n",k2SGC[j],P2SGC[j],ptheo);
        j++;
        if(j==NeffP2baoSGC){j=0;modeP2bao=0;}
        }
        if(modeP0bao==1){
Ninterpol=determine_N_singlearray(k_theobao,k0baoSGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC_bao);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0baoSGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theobao,k0baoSGC[j],Ninterpol,spacing_dataSGC_bao);

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theobao,k0baoSGC[j],Ninterpol,spacing_dataSGC_bao);
w1=determine_w1_2ndorder_singlearray(k_theobao,k0baoSGC[j],Ninterpol,spacing_dataSGC_bao);
w2=determine_w2_2ndorder_singlearray(k_theobao,k0baoSGC[j],Ninterpol,spacing_dataSGC_bao);
}
ptheo=P_interpol_fast(k0baoSGC[j],P_theo0bao,Neffmax*factor_sampling_mask,spacing_dataSGC_bao,interpolation_order,Ninterpol,w0,w1,w2);
}


        pobs=P0baoSGC[j];//printf("%lf %lf %lf\n",k0SGC[j],P0SGC[j],ptheo);
        j++;
        if(j==NeffP0baoSGC){j=0;modeP0bao=0;}
        }

        difference[i]=ptheo-pobs;
   }
}
//RSD aqui

if(strcmp(path_to_mask2, "none") == 0)
{
if(NeffP0rsdSGC*modeP0rsd>0){k_theo0rsd = (double*) calloc( NeffP0rsdSGC, sizeof(double));}
if(NeffP2rsdSGC*modeP2rsd>0){k_theo2rsd = (double*) calloc( NeffP2rsdSGC, sizeof(double));}
if(NeffP4rsdSGC*modeP4rsd>0){k_theo4rsd = (double*) calloc( NeffP4rsdSGC, sizeof(double));}
if(NeffP0rsdSGC*modeP0rsd>0){P_theo0rsd = (double*) calloc( NeffP0rsdSGC, sizeof(double));}
if(NeffP2rsdSGC*modeP2rsd>0){P_theo2rsd = (double*) calloc( NeffP2rsdSGC, sizeof(double));}
if(NeffP4rsdSGC*modeP4rsd>0){P_theo4rsd = (double*) calloc( NeffP4rsdSGC, sizeof(double));}
}
else{
if( strcmp(spacing_dataSGC_rsd,"linear") == 0  ){

if(NeffP0rsdSGC*modeP0rsd>=NeffP2rsdSGC*modeP2rsd){Neffmax=NeffP0rsdSGC-2+25+(int)(k0rsdSGC[0]/(k0rsdSGC[NeffP0rsdSGC-1]-k0rsdSGC[0])*(NeffP0rsdSGC-1));}
if(NeffP0rsdSGC*modeP0rsd>=NeffP4rsdSGC*modeP4rsd){Neffmax=NeffP0rsdSGC-2+25+(int)(k0rsdSGC[0]/(k0rsdSGC[NeffP0rsdSGC-1]-k0rsdSGC[0])*(NeffP0rsdSGC-1));}

if(NeffP2rsdSGC*modeP2rsd>=NeffP0rsdSGC*modeP0rsd){Neffmax=NeffP2rsdSGC-2+25+(int)(k2rsdSGC[0]/(k2rsdSGC[NeffP2rsdSGC-1]-k2rsdSGC[0])*(NeffP2rsdSGC-1));}
if(NeffP2rsdSGC*modeP2rsd>=NeffP4rsdSGC*modeP4rsd){Neffmax=NeffP2rsdSGC-2+25+(int)(k2rsdSGC[0]/(k2rsdSGC[NeffP2rsdSGC-1]-k2rsdSGC[0])*(NeffP2rsdSGC-1));}

if(NeffP4rsdSGC*modeP4rsd>=NeffP0rsdSGC*modeP0rsd){Neffmax=NeffP4rsdSGC-2+25+(int)(k4rsdSGC[0]/(k4rsdSGC[NeffP4rsdSGC-1]-k4rsdSGC[0])*(NeffP4rsdSGC-1));}
if(NeffP4rsdSGC*modeP4rsd>=NeffP2rsdSGC*modeP2rsd){Neffmax=NeffP4rsdSGC-2+25+(int)(k4rsdSGC[0]/(k4rsdSGC[NeffP4rsdSGC-1]-k4rsdSGC[0])*(NeffP4rsdSGC-1));}

}

if( strcmp(spacing_dataSGC_rsd,"log") == 0  ){

if(NeffP0rsdSGC*modeP0rsd>=NeffP2rsdSGC*modeP2rsd){Neffmax=NeffP0rsdSGC-2+25+(int)(log(k0rsdSGC[0])/(log(k0rsdSGC[NeffP0rsdSGC-1])-log(k0rsdSGC[0]))*(NeffP0rsdSGC-1));}
if(NeffP0rsdSGC*modeP0rsd>=NeffP4rsdSGC*modeP4rsd){Neffmax=NeffP0rsdSGC-2+25+(int)(log(k0rsdSGC[0])/(log(k0rsdSGC[NeffP0rsdSGC-1])-log(k0rsdSGC[0]))*(NeffP0rsdSGC-1));}

if(NeffP2rsdSGC*modeP2rsd>=NeffP0rsdSGC*modeP0rsd){Neffmax=NeffP2rsdSGC-2+25+(int)(log(k2rsdSGC[0])/(log(k2rsdSGC[NeffP2rsdSGC-1])-log(k2rsdSGC[0]))*(NeffP2rsdSGC-1));}
if(NeffP2rsdSGC*modeP2rsd>=NeffP4rsdSGC*modeP4rsd){Neffmax=NeffP2rsdSGC-2+25+(int)(log(k2rsdSGC[0])/(log(k2rsdSGC[NeffP2rsdSGC-1])-log(k2rsdSGC[0]))*(NeffP2rsdSGC-1));}

if(NeffP4rsdSGC*modeP4rsd>=NeffP0rsdSGC*modeP0rsd){Neffmax=NeffP4rsdSGC-2+25+(int)(log(k4rsdSGC[0])/(log(k4rsdSGC[NeffP4rsdSGC-1])-log(k4rsdSGC[0]))*(NeffP4rsdSGC-1));}
if(NeffP4rsdSGC*modeP4rsd>=NeffP2rsdSGC*modeP2rsd){Neffmax=NeffP4rsdSGC-2+25+(int)(log(k4rsdSGC[0])/(log(k4rsdSGC[NeffP4rsdSGC-1])-log(k4rsdSGC[0]))*(NeffP4rsdSGC-1));}

}

if( strcmp(spacing_dataSGC_rsd,"log10") == 0  ){

if(NeffP0rsdSGC*modeP0rsd>=NeffP2rsdSGC*modeP2rsd){Neffmax=NeffP0rsdSGC-2+25+(int)(log10(k0rsdSGC[0])/(log10(k0rsdSGC[NeffP0rsdSGC-1])-log10(k0rsdSGC[0]))*(NeffP0rsdSGC-1));}
if(NeffP0rsdSGC*modeP0rsd>=NeffP4rsdSGC*modeP4rsd){Neffmax=NeffP0rsdSGC-2+25+(int)(log10(k0rsdSGC[0])/(log10(k0rsdSGC[NeffP0rsdSGC-1])-log10(k0rsdSGC[0]))*(NeffP0rsdSGC-1));}

if(NeffP2rsdSGC*modeP2rsd>=NeffP0rsdSGC*modeP0rsd){Neffmax=NeffP2rsdSGC-2+25+(int)(log10(k2rsdSGC[0])/(log10(k2rsdSGC[NeffP2rsdSGC-1])-log10(k2rsdSGC[0]))*(NeffP2rsdSGC-1));}
if(NeffP2rsdSGC*modeP2rsd>=NeffP4rsdSGC*modeP4rsd){Neffmax=NeffP2rsdSGC-2+25+(int)(log10(k2rsdSGC[0])/(log10(k2rsdSGC[NeffP2rsdSGC-1])-log10(k2rsdSGC[0]))*(NeffP2rsdSGC-1));}

if(NeffP4rsdSGC*modeP4rsd>=NeffP0rsdSGC*modeP0rsd){Neffmax=NeffP4rsdSGC-2+25+(int)(log10(k4rsdSGC[0])/(log10(k4rsdSGC[NeffP4rsdSGC-1])-log10(k4rsdSGC[0]))*(NeffP4rsdSGC-1));}
if(NeffP4rsdSGC*modeP4rsd>=NeffP2rsdSGC*modeP2rsd){Neffmax=NeffP4rsdSGC-2+25+(int)(log10(k4rsdSGC[0])/(log10(k4rsdSGC[NeffP4rsdSGC-1])-log10(k4rsdSGC[0]))*(NeffP4rsdSGC-1));}

}
if( strcmp(spacing_dataSGC_rsd,"irregular") == 0  ){
Neffmax=Ntheory;
factor_sampling_mask=1;
}


k_theorsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0rsd==1){P_theo0rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2rsd==1){P_theo2rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4rsd==1){P_theo4rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
}


parameters1 =  (double*) calloc( dimension, sizeof(double));

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<3){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}
}

if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2_rsd[i1];
i1++;
}
if(i==3 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=4 && i<=6){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}

if(i==7 && strcmp(local_b2s2, "yes") == 0){
parameters1[i]=-4./7.*(parameters1[i-3]-1);
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}
if(i==8 && strcmp(local_b3nl, "yes") == 0){
parameters1[i]=32./315.*(parameters1[i-4]-1);
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}

if(i==9 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=0;
}

if(i==9 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){
parameters1[i]=parameters1[i-1];
}

if(i==10 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){
parameters1[i]=0;
}


}

//for(i=0;i<dimension;i++){printf("RSD-SGC: %d %lf\n",i,parameters1[i]);}

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps,type_of_analysis, RSD_fit, fit_RSD,modeP0rsd,modeP2rsd,modeP4rsd,Theory,Ntheory, k_theorsd,k_theo0rsd,k_theo2rsd,k_theo4rsd, P_theo0rsd, P_theo2rsd, P_theo4rsd,PnoiseSGC,NeffP0rsdSGC,NeffP2rsdSGC,NeffP4rsdSGC,factor_sampling_mask, parameters1,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC,k0rsdSGC,k2rsdSGC,k4rsdSGC, path_to_mask2,plan1rsd, plan2rsd, k0rsdSGC[0], k0rsdSGC[NeffP0rsdSGC-1], k2rsdSGC[0], k2rsdSGC[NeffP2rsdSGC-1], k4rsdSGC[0], k4rsdSGC[NeffP4rsdSGC-1],spacing_dataSGC_rsd,spacing_theory_rsd);

free(parameters1);

if(strcmp(path_to_mask2, "none") == 0 )
{
i=-1;

        if(modeP0rsd==1){
        for(j=0;j<NeffP0rsdSGC;j++){
        ptheo=P_theo0rsd[j]-PnoiseSGC;
        pobs=P0rsdSGC[j];
        i++;
        difference[i+Ncovbao]=ptheo-pobs;
        }
        }

        if(modeP2rsd==1){
        for(j=0;j<NeffP2rsdSGC;j++){
        ptheo=P_theo2rsd[j];
        pobs=P2rsdSGC[j];
        i++;
        difference[i+Ncovbao]=ptheo-pobs;
        }
        }

        if(modeP4rsd==1){
        for(j=0;j<NeffP4rsdSGC;j++){
        ptheo=P_theo4rsd[j];
        pobs=P4rsdSGC[j];
        i++;
        difference[i+Ncovbao]=ptheo-pobs;
        }
        }
}
else{


j=0;
    for(i=0;i<Ncovrsd;i++)
    {

        if(modeP4rsd==1 && modeP0rsd==0 && modeP2rsd==0){
Ninterpol=determine_N_singlearray(k_theorsd,k4rsdSGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC_rsd);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4rsdSGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theorsd,k4rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theorsd,k4rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
w1=determine_w1_2ndorder_singlearray(k_theorsd,k4rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
w2=determine_w2_2ndorder_singlearray(k_theorsd,k4rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
}
ptheo=P_interpol_fast(k4rsdSGC[j],P_theo4rsd,Neffmax*factor_sampling_mask,spacing_dataSGC_rsd,interpolation_order,Ninterpol,w0,w1,w2);

}


        pobs=P4rsdSGC[j];
        j++;
        if(j==NeffP4rsdSGC){j=0;modeP4rsd=0;}
        }

        if(modeP2rsd==1 && modeP0rsd==0){
Ninterpol=determine_N_singlearray(k_theorsd,k2rsdSGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC_rsd);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2rsdSGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theorsd,k2rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theorsd,k2rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
w1=determine_w1_2ndorder_singlearray(k_theorsd,k2rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
w2=determine_w2_2ndorder_singlearray(k_theorsd,k2rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
}
ptheo=P_interpol_fast(k2rsdSGC[j],P_theo2rsd,Neffmax*factor_sampling_mask,spacing_dataSGC_rsd,interpolation_order,Ninterpol,w0,w1,w2);
}



        pobs=P2rsdSGC[j];//printf("%lf %lf %lf\n",k2SGC[j],P2SGC[j],ptheo);
        j++;
        if(j==NeffP2rsdSGC){j=0;modeP2rsd=0;}
        }

        if(modeP0rsd==1){
Ninterpol=determine_N_singlearray(k_theorsd,k0rsdSGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC_rsd);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0rsdSGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theorsd,k0rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theorsd,k0rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
w1=determine_w1_2ndorder_singlearray(k_theorsd,k0rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
w2=determine_w2_2ndorder_singlearray(k_theorsd,k0rsdSGC[j],Ninterpol,spacing_dataSGC_rsd);
}
ptheo=P_interpol_fast(k0rsdSGC[j],P_theo0rsd,Neffmax*factor_sampling_mask,spacing_dataSGC_rsd,interpolation_order,Ninterpol,w0,w1,w2)-PnoiseSGC;
}


        pobs=P0rsdSGC[j];//printf("%lf %lf %lf\n",k0SGC[j],P0SGC[j],ptheo);
        j++;
        if(j==NeffP0rsdSGC){j=0;modeP0rsd=0;}
        }

        difference[i+Ncovbao]=ptheo-pobs;
   }


}



if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0bao=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2bao=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4bao=1;}

if(strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP0rsd=1;}
if(strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP2rsd=1;}
if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0 ){modeP4rsd=1;}


ch2SGC=0;
//chi2_baoSGC=0;
//chi2_rsdSGC=0;
for(i=0;i<Ncov;i++)
{

   for(j=0;j<Ncov;j++)
      {
                  ch2SGC=ch2SGC+difference[i]*1./covSGC[i+Ncov*j]*difference[j];

//if(i<Ncovbao && j<Ncovbao){chi2_baoSGC=chi2_baoSGC+difference[i]*1./covSGC[i+Ncov*j]*difference[j];}
//if(i>=Ncovbao && j>=Ncovbao){chi2_rsdSGC=chi2_rsdSGC+difference[i]*1./covSGC[i+Ncov*j]*difference[j];}


      }

}
if(strcmp(path_to_mask2, "none") != 0 ){free(k_theobao);}
if(strcmp(path_to_mask2, "none") != 0 ){free(k_theorsd);}

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo0bao);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo0bao);}
}
if(strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0)
{
free(P_theo0rsd);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo0rsd);}
}


if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo2bao);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo2bao);}
}
if(strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0)
{
free(P_theo2rsd);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo2rsd);}
}

if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo4bao);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo4bao);}
}
if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0)
{
free(P_theo4rsd);
if(strcmp(path_to_mask2, "none") == 0 ){free(k_theo4rsd);}
}


free(difference);


ch2=ch2+ch2SGC;
//chi2_bao=chi2_bao+chi2_baoSGC;
//chi2_rsd=chi2_rsd+chi2_rsdSGC;

}//Nchunks=2


for(i=0;i<Nsigmas_tot;i++)
{

if(Sigma_type[i]==1){prior_chi2=gauss(parameters2_bao[Nalphas+i],Sigma_nl_mean[i],Sigma_nl_stddev[i]);ch2=ch2+prior_chi2;}

}

//printf("%lf %lf %lf (%d)\n",ch2-ch2SGC-prior_chi2,ch2SGC,prior_chi2,Ncov);
//printf("%lf %lf, %lf %lf, %lf %lf\n",ch2,ch2SGC,chi2_bao,chi2_baoSGC,chi2_rsd,chi2_rsdSGC);
//printf("%lf %lf %lf\n",ch2,ch2-ch2SGC,ch2SGC);
//exit(0);
return ch2;
}

double chi2_rsd_mcmc(char *type_of_analysis,double *parameters2, double **Theory,int N_Plin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC,double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0, double *P0, double Pnoise, double *k2, double *P2, double *k4, double *P4, double *k0SGC, double *P0SGC, double PnoiseSGC, double *k2SGC, double *P2SGC, double *k4SGC, double *P4SGC, int NeffP0, int NeffP2, int NeffP4,int NeffP0SGC, int NeffP2SGC, int NeffP4SGC, double *cov, double *covSGC,  char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free, char *fog_bs, int Nchunks, fftw_plan plan1, fftw_plan plan2, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum, int factor_sampling_mask_in,char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory)
{
double ch2,ch2SGC,prior_chi2;
double *difference,*parameters1,*k_theo,*P_theo0,*P_theo2,*P_theo4;
double *k_theo0,*k_theo2,*k_theo4;
int i,j,Ndifference;
int i1;
double ptheo,pobs;
int dimension;
int Ncov;
int modeP0,modeP2,modeP4;
int Ntheo=N_Plin;
int offset;
int Neffmax;
int factor_sampling_mask;
int interpolation_order,shiftN;
double w1,w2,w0;
int Ninterpol;

interpolation_order=1;

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}



factor_sampling_mask=factor_sampling_mask_in;

offset=3;
if(strcmp(local_b2s2, "no") == 0){offset++;}//b2s2
if(strcmp(local_b3nl, "no") == 0){offset++;}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){offset++;}//fog_ps
if(strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 ){offset++;}//fog_bs

dimension=11;//This is a fixed quantity.

modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

Ncov=NeffP0*modeP0+NeffP2*modeP2+NeffP4*modeP4;

if(strcmp(path_to_mask1, "none") == 0)
{
if(NeffP0*modeP0>0){k_theo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2*modeP2>0){k_theo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4*modeP4>0){k_theo4 = (double*) calloc( NeffP4, sizeof(double));}
if(NeffP0*modeP0>0){P_theo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2*modeP2>0){P_theo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4*modeP4>0){P_theo4 = (double*) calloc( NeffP4, sizeof(double));}
}
else{
//k_theo = (double*) calloc( Ntheo, sizeof(double));
//if(modeP0==1){P_theo0 = (double*) calloc( Ntheo, sizeof(double));}
//if(modeP2==1){P_theo2 = (double*) calloc( Ntheo, sizeof(double));}
//if(modeP4==1){P_theo4 = (double*) calloc( Ntheo, sizeof(double));}

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


k_theo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0==1){P_theo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2==1){P_theo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4==1){P_theo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}


}

difference = (double*) calloc( Ncov, sizeof(double));

parameters1 = (double*) calloc( dimension, sizeof(double));

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
parameters1[i]=parameters2[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<3){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}
}

if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==3 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=4 && i<=6){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==7 && strcmp(local_b2s2, "yes") == 0){
parameters1[i]=-4./7.*(parameters1[i-3]-1);
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==8 && strcmp(local_b3nl, "yes") == 0){
parameters1[i]=32./315.*(parameters1[i-4]-1);
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==9 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=0;
}

if(i==9 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){
parameters1[i]=parameters1[i-1];
}

if(i==10 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){
parameters1[i]=0;
}


}


//for(i=0;i<dimension;i++){printf("NGC %d %lf\n",i,parameters1[i]);}

//call here the Pmodel 

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,Pnoise, NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC, k0,k2,k4,path_to_mask1,plan1, plan2, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1],spacing_dataNGC,spacing_theory);
//printf("%lf %lf %lf %lf\n",k_theo[300],P_theo0[300],P_theo2[300],P_theo4[300]);
//printf("%lf %lf %lf\n",P_theo0[10],P_theo2[10],P_theo4[10]);
free(parameters1);

if(strcmp(path_to_mask1, "none") == 0 )
{
i=-1;

   if(modeP0==1){
        for(j=0;j<NeffP0;j++){
        ptheo=P_theo0[j]-Pnoise;
        pobs=P0[j];
        i++;
        difference[i]=ptheo-pobs;//printf("%lf %lf %lf\n",k_theo0[j],P0[j],P_theo0[j]);
        }
        }

        if(modeP2==1){
        for(j=0;j<NeffP2;j++){
        ptheo=P_theo2[j];
        pobs=P2[j];
        i++;
        difference[i]=ptheo-pobs;
        }
        }

        if(modeP4==1){
        for(j=0;j<NeffP4;j++){
        ptheo=P_theo4[j];
        pobs=P4[j];
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

//        ptheo=P_interpol(k4[j],k_theo,P_theo4,Neffmax*factor_sampling_mask);
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
//        ptheo=P_interpol(k2[j],k_theo,P_theo2,Neffmax*factor_sampling_mask);
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
//        ptheo=P_interpol(k0[j],k_theo,P_theo0,Neffmax*factor_sampling_mask)-Pnoise;
Ninterpol=determine_N_singlearray(k_theo,k0[j],Neffmax*factor_sampling_mask,spacing_dataNGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso(P_theo0,Ninterpol,w1)-Pnoise;

}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k0[j],Ninterpol,spacing_dataNGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo0,Ninterpol,w0,w1,w2)-Pnoise;
}
ptheo=P_interpol_fast(k0[j],P_theo0,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2)-Pnoise;
}


        pobs=P0[j];
        j++;
        if(j==NeffP0){j=0;modeP0=0;}
        }

        difference[i]=ptheo-pobs;
   }

}
//exit(0);
/*
j=0;
    for(i=0;i<Ncov;i++)
    {

        if(modeP4==1 && modeP0==0 && modeP2==0){
        ptheo=P_interpol(k4[j],k_theo,P_theo4,Ntheo);
        pobs=P4[j];//printf("P4: %lf %lf %lf (%d,%d)\n",k4[j],ptheo,pobs,i,j);
        j++;
        if(j==NeffP4){j=0;modeP4=0;}
        }

        if(modeP2==1 && modeP0==0){
        ptheo=P_interpol(k2[j],k_theo,P_theo2,Ntheo);
        pobs=P2[j];//printf("P2: %lf %lf %lf (%d,%d)\n",k2[j],ptheo,pobs,i,j);
        j++;
        if(j==NeffP2){j=0;modeP2=0;}
        }

        if(modeP0==1){
        ptheo=P_interpol(k0[j],k_theo,P_theo0,Ntheo)-Pnoise;
        pobs=P0[j];//printf("P0: %lf %lf %lf (%d,%d)\n",k0[j],ptheo,pobs,i,j);
        j++;
        if(j==NeffP0){j=0;modeP0=0;}
        }

        difference[i]=ptheo-pobs;
      }
*/
ch2=0;
ch2SGC=0;
for(i=0;i<Ncov;i++)
{

   for(j=0;j<Ncov;j++)
      {
            ch2=ch2+difference[i]*1./cov[i+Ncov*j]*difference[j];

//if(i==j){printf("%d %lf %lf\n",i,difference[i],sqrt(cov[i+Ncov*j]));}

      }

}
//printf("chi2=%lf\n",ch2);
//exit(0);

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


if(Nchunks==2)//repetate for Nchunk=2
{
factor_sampling_mask=factor_sampling_mask_in;
modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

if(strcmp(path_to_mask2, "none") == 0)
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
if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){Neffmax=NeffP0SGC+1;}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){Neffmax=NeffP0SGC+1;}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){Neffmax=NeffP2SGC+1;}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){Neffmax=NeffP2SGC+1;}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){Neffmax=NeffP4SGC+1;}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){Neffmax=NeffP4SGC+1;}
*/

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


Ncov=NeffP0SGC*modeP0+NeffP2SGC*modeP2+NeffP4SGC*modeP4;

difference = (double*) calloc( Ncov, sizeof(double));
parameters1 =  (double*) calloc( dimension, sizeof(double));

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
parameters1[i]=parameters2[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<3){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}
}

if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==3 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=4 && i<=6){
parameters1[i]=parameters2[i1+offset];
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2[i1+offset];
i1++;
}

if(i==7 && strcmp(local_b2s2, "yes") == 0){
parameters1[i]=-4./7.*(parameters1[i-3]-1);
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2[i1+offset];
i1++;
}
if(i==8 && strcmp(local_b3nl, "yes") == 0){
parameters1[i]=32./315.*(parameters1[i-4]-1);
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2[i1+offset];
i1++;
}

if(i==9 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=0;
}

if(i==9 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){
parameters1[i]=parameters2[i1+offset];
i1++;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){
parameters1[i]=parameters1[i-1];
}

if(i==10 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){
parameters1[i]=0;
}


}

//for(i=0;i<dimension;i++){printf("SGC %d %lf\n",i,parameters1[i]);}

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps,type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,PnoiseSGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC,k0SGC,k2SGC,k4SGC, path_to_mask2,plan1, plan2, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1],spacing_dataSGC,spacing_theory);



free(parameters1);

if(strcmp(path_to_mask2, "none") == 0 )
{
i=-1;

        if(modeP0==1){
        for(j=0;j<NeffP0SGC;j++){
        ptheo=P_theo0[j]-PnoiseSGC;
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
        //ptheo=P_interpol(k4SGC[j],k_theo,P_theo4,Neffmax*factor_sampling_mask);
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
//        ptheo=P_interpol(k0SGC[j],k_theo,P_theo0,Neffmax*factor_sampling_mask)-PnoiseSGC;
Ninterpol=determine_N_singlearray(k_theo,k0SGC[j],Neffmax*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0SGC[j]<=0)
{
ptheo=0;
}
else{

if(interpolation_order==1)
{
w1=determine_w1_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso(P_theo0,Ninterpol,w1)-PnoiseSGC;
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(k_theo,k0SGC[j],Ninterpol,spacing_dataSGC);
//ptheo=P_interpol2_aniso_2ndorder(P_theo0,Ninterpol,w0,w1,w2)-PnoiseSGC;
}
ptheo=P_interpol_fast(k0SGC[j],P_theo0,Neffmax*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2)-PnoiseSGC;
}


        pobs=P0SGC[j];//printf("%lf %lf %lf\n",k0SGC[j],P0SGC[j],ptheo);
        j++;
        if(j==NeffP0SGC){j=0;modeP0=0;}
        }

        difference[i]=ptheo-pobs;
   }


}



/*
j=0;
    for(i=0;i<Ncov;i++)
    {

        if(modeP4==1 && modeP2==0 && modeP0==0){
        ptheo=P_interpol(k4SGC[j],k_theo,P_theo4,Ntheo);
        pobs=P4SGC[j];
        j++;
        if(j==NeffP4){j=0;modeP4=0;}
        }


        if(modeP2==1 && modeP0==0){
        ptheo=P_interpol(k2SGC[j],k_theo,P_theo2,Ntheo);
        pobs=P2SGC[j];
        j++;
        if(j==NeffP2){j=0;modeP2=0;}
        }

        if(modeP0==1){
        ptheo=P_interpol(k0SGC[j],k_theo,P_theo0,Ntheo)-PnoiseSGC;
        pobs=P0SGC[j];
        j++;
        if(j==NeffP0SGC){j=0;modeP0=0;}
        }


        difference[i]=ptheo-pobs;
   }
*/
ch2SGC=0;
for(i=0;i<Ncov;i++)
{

   for(j=0;j<Ncov;j++)
      {
            ch2SGC=ch2SGC+difference[i]*1./covSGC[i+Ncov*j]*difference[j];
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


ch2=ch2+ch2SGC;
}//Nchunks==2

//printf("%lf %lf\n",ch2,ch2SGC);
//exit(0);

return ch2;
}

void do_rsd_bao_mcmc(int nthreads,char *type_BAO_fit,char *type_of_analysis,char *fit_BAO,char *fit_RSD,double **Theory,int Ntheory,double *k_Plin,double *Plin,int N_Plin, double *k_Olin, double *Olin, int N_Olin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0bao, double *k0rsd, double *P0bao, double *P0rsd, double *errP0bao, double *errP0rsd,double Pnoise_rsd, int NeffP0bao, int NeffP0rsd, double *k2bao, double *k2rsd, double *P2bao, double *P2rsd, double *errP2bao, double *errP2rsd, int NeffP2bao, int NeffP2rsd, double *k4bao, double *k4rsd, double *P4bao, double *P4rsd, double *errP4bao, double *errP4rsd, int NeffP4bao, int NeffP4rsd, double *k11bao, double *k11rsd, double *k22bao, double *k22rsd, double *k33bao, double *k33rsd, double *B0bao, double *B0rsd, double *errB0bao, double *errB0rsd, double *Bnoise_bao, double *Bnoise_rsd, int NeffB0bao, int NeffB0rsd, double *k0baoSGC, double *k0rsdSGC, double *P0baoSGC, double *P0rsdSGC, double *errP0baoSGC, double *errP0rsdSGC,double Pnoise_rsdSGC, int NeffP0baoSGC,int NeffP0rsdSGC, double *k2baoSGC, double *k2rsdSGC, double *P2baoSGC, double *P2rsdSGC, double *errP2baoSGC, double *errP2rsdSGC, int NeffP2baoSGC, int NeffP2rsdSGC, double *k4baoSGC, double *k4rsdSGC, double *P4baoSGC, double *P4rsdSGC, double *errP4baoSGC, double *errP4rsdSGC, int NeffP4baoSGC, int NeffP4rsdSGC, double *k11baoSGC,double *k11rsdSGC, double *k22baoSGC,double *k22rsdSGC, double *k33baoSGC,double *k33rsdSGC,double *B0baoSGC,double *B0rsdSGC, double *errB0baoSGC, double *errB0rsdSGC, double *Bnoise_baoSGC, double *Bnoise_rsdSGC,int NeffB0baoSGC,int NeffB0rsdSGC, double *cov, double *covSGC, double alpha_min, double alpha_max, char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int  Npolynomial, char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit,char *sigma8_free,char *fog_free,char *fog_bs, int  Nchunks, char *path_output, char *identifier, char *do_plot, char *use_prop_cov, char *path_to_cov, int Nsteps, char *do_power_spectrum, char *do_bispectrum, double Sigma_smooth, char *spacing_dataNGC_bao, char *spacing_dataNGC_rsd, char *spacing_dataSGC_bao, char *spacing_dataSGC_rsd, char *spacing_theory,char *spacing_theory_rsd,char *type_of_analysis_BAO, char *type_of_analysis_FS)
{
double fraction;
int trial_mcmc;
long int N_max,N_print,N_burnout,j;
double **vector_buffer;
FILE *f;
double *Cov_prop;
double *vector_mean;
int N_Cov_prop,N_Cov_prop_bao,N_Cov_prop_rsd;
gsl_matrix *transform_inverse;
gsl_matrix *transform;
char name_file[2000];
long int *params_mcmc;
int modeP0bao,modeP2bao,modeP4bao;
int modeP0rsd,modeP2rsd,modeP4rsd;
int allsigmafixed,Nalphas,Nsigmas_tot,Nsigmas_free;
int factor_sampling_mask;
int i_thread;


factor_sampling_mask=10;//Sampling factor boost wrt to Neff max when the mask is applied. How much do I need to sample my model pre-mask-apply in order to have a satisfactory mask response? 10 times seems reasonable.

//ini BAO

if (nthreads<2){

  if(strcmp(type_of_analysis, "FSBAOISO") == 0){sprintf(name_file,"%s/mcmcFSBAOISO_output_%s.txt",path_output,identifier);}
  if(strcmp(type_of_analysis, "FSBAOANISO") == 0){sprintf(name_file,"%s/mcmcFSBAOANISO_output_%s.txt",path_output,identifier);}
  f=fopen(name_file,"w");
  if(f==NULL){printf("Error writing %s. Exiting now...\n",name_file);exit(0);}
  fclose(f);
} else {
  for (i_thread=0;i_thread<nthreads;i_thread++){

    if(strcmp(type_of_analysis, "FSBAOISO") == 0){sprintf(name_file,"%s/mcmcFSBAOISO_output_%s__%d.txt",path_output,identifier,i_thread+1);}
    if(strcmp(type_of_analysis, "FSBAOANISO") == 0){sprintf(name_file,"%s/mcmcFSBAOANISO_output_%s__%d.txt",path_output,identifier,i_thread+1);}
    f=fopen(name_file,"w");
    if(f==NULL){printf("Error writing %s. Exiting now...\n",name_file);exit(0);}
    fclose(f);
  }
  
  if(strcmp(type_of_analysis, "FSBAOISO") == 0){sprintf(name_file,"%s/mcmcFSBAOISO_output_%s.txt",path_output,identifier);}
  else if(strcmp(type_of_analysis, "FSBAOANISO") == 0){sprintf(name_file,"%s/mcmcFSBAOANISO_output_%s.txt",path_output,identifier);}
}

modeP0bao=0;
modeP2bao=0;
modeP4bao=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0bao=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2bao=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4bao=1;}

allsigmafixed=-1;
Nalphas=1;if(modeP0bao+modeP2bao+modeP4bao>1){Nalphas=2;}
Nsigmas_free=0;
if(strcmp(Sigma_independent, "yes") == 0 ){

    if(strcmp(Sigma_def_type, "para-perp") == 0)//This is the only possible case for BAOANISO
    {
       if(Sigma_type[0]>0){Nsigmas_free=Nsigmas_free+1;}
       if(Sigma_type[1]>0){Nsigmas_free=Nsigmas_free+1;}
    }

    if(strcmp(Sigma_def_type, "effective") == 0)
    {
        if(modeP0bao==1 && Sigma_type[0]>0){Nsigmas_free=Nsigmas_free+1;}
        if(modeP2bao==1 && Sigma_type[1]>0){Nsigmas_free=Nsigmas_free+1;}
        if(modeP4bao==1 && Sigma_type[2]>0){Nsigmas_free=Nsigmas_free+1;}
    }

}
if(strcmp(Sigma_independent, "no") == 0 ){
if(Sigma_type[0]>0){Nsigmas_free=1;}
if(Sigma_type[0]==0){Nsigmas_free=0;}
}
if(Nsigmas_free==0){allsigmafixed=0;}

//here the true value of Nalphas is always 2 when FSBAO
if(strcmp(type_of_analysis, "FSBAOISO") == 0){N_Cov_prop_bao=(Npolynomial+1)*Nchunks*(modeP0bao+modeP2bao+modeP4bao)+2+Nsigmas_free;}
if(strcmp(type_of_analysis, "FSBAOANISO") == 0){N_Cov_prop_bao=Nchunks*(1+Npolynomial*(modeP0bao+modeP2bao+modeP4bao))+2+Nsigmas_free+1;}
//if(strcmp(type_of_analysis, "FSBAOISO") == 0){N_Cov_prop_bao=(Npolynomial+1)*Nchunks*(modeP0bao+modeP2bao+modeP4bao)+Nalphas+Nsigmas_free;}
//if(strcmp(type_of_analysis, "FSBAOANISO") == 0){N_Cov_prop_bao=Nchunks*(1+Npolynomial*(modeP0bao+modeP2bao+modeP4bao))+Nalphas+Nsigmas_free+1;}

if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0){Nsigmas_tot=2;}
if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0){Nsigmas_tot=modeP0bao+modeP2bao+modeP4bao;}
if( strcmp(Sigma_independent, "no") == 0 ){Nsigmas_tot=1;}

//ini RSD
N_Cov_prop_rsd=3;//b1, b2, A

if(strcmp(local_b2s2, "no") == 0){N_Cov_prop_rsd++;}
if(strcmp(local_b3nl, "no") == 0){N_Cov_prop_rsd++;}
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){N_Cov_prop_rsd++;}
if(strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 ){N_Cov_prop_rsd++;}
N_Cov_prop_rsd=N_Cov_prop_rsd*Nchunks;//x2 
if(strcmp(sigma8_free, "yes") == 0){N_Cov_prop_rsd++;}
if(strcmp(RSD_fit, "yes") == 0){N_Cov_prop_rsd=N_Cov_prop_rsd+3;}//f, alpha-para, alpha-perp


N_Cov_prop=N_Cov_prop_rsd+N_Cov_prop_bao-2;
//printf("%d,%d,%d %s\n",N_Cov_prop,N_Cov_prop_bao,N_Cov_prop_rsd,use_prop_cov);
//exit(0);
//read prop. cov.
if(strcmp(use_prop_cov, "yes") == 0)//full mcmc run with proposal
{
trial_mcmc=0;//no trial
fraction=1;
Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

read_prop_cov(NULL,0,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);
transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);

generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse

free(Cov_prop);

mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,NULL,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,fit_RSD,k_Plin,Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao,k0rsd, P0bao, P0rsd, errP0bao,errP0rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, errP2bao,errP2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, errP4bao,errP4rsd, NeffP4bao,NeffP4rsd, k11bao,k11rsd, k22bao,k2rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, errP0baoSGC,errP0rsdSGC, NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, errP2baoSGC,errP2rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, errP4baoSGC,errP4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC, B0baoSGC,B0rsdSGC,Bnoise_baoSGC,Bnoise_rsdSGC, NeffB0baoSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent,ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot, Nsigmas_free,Theory,Ntheory, Pnoise_rsd,Pnoise_rsdSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , Sigma_smooth,factor_sampling_mask, spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS);
}
else//two iterations: first trial to get a proposal covariance, seccond full run with such proposal
{
trial_mcmc=1;//trial iteration
fraction=0.10;
params_mcmc=(long int*) calloc( 2, sizeof(long int ));
set_mcmc_parameters(params_mcmc);
N_burnout=params_mcmc[1];
free(params_mcmc);

N_max=(long int)(Nsteps*fraction);
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

//use proposed mean
set_proposal_mean(vector_mean,N_Cov_prop);

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);

generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse

vector_buffer = (double**) calloc(N_Cov_prop+2,sizeof(double*));

for(j=0;j<N_Cov_prop+2;j++)
{
   vector_buffer[j] = (double*) calloc(N_max,sizeof(double));
}

mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,fit_RSD,k_Plin,Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao,k0rsd, P0bao, P0rsd, errP0bao,errP0rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, errP2bao,errP2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, errP4bao,errP4rsd, NeffP4bao,NeffP4rsd, k11bao,k11rsd, k22bao,k2rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, errP0baoSGC,errP0rsdSGC, NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, errP2baoSGC,errP2rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, errP4baoSGC,errP4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC, B0baoSGC,B0rsdSGC,Bnoise_baoSGC,Bnoise_rsdSGC, NeffB0baoSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent,ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot, Nsigmas_free,Theory,Ntheory, Pnoise_rsd,Pnoise_rsdSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , Sigma_smooth,factor_sampling_mask, spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS);


Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized
//get covariance from buffer
read_prop_cov(vector_buffer,N_max,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);

freeTokens(vector_buffer,N_Cov_prop+2);

trial_mcmc=0;//no trial
fraction=1;

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse

free(Cov_prop);
mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,fit_RSD,k_Plin,Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao,k0rsd, P0bao, P0rsd, errP0bao,errP0rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, errP2bao,errP2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, errP4bao,errP4rsd, NeffP4bao,NeffP4rsd, k11bao,k11rsd, k22bao,k2rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, errP0baoSGC,errP0rsdSGC, NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, errP2baoSGC,errP2rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, errP4baoSGC,errP4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC, B0baoSGC,B0rsdSGC,Bnoise_baoSGC,Bnoise_rsdSGC, NeffB0baoSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent,ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot, Nsigmas_free,Theory,Ntheory, Pnoise_rsd,Pnoise_rsdSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , Sigma_smooth,factor_sampling_mask, spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS);

}


}


void do_rsd_mcmc(int nthreads, char *type_of_analysis, char *fit_BAO,double **Theory,int N_Plin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC,  char *path_to_mask2, char *spacing_maskSGC, double *k0, double *P0, double *errP0, double Pnoise, int NeffP0, double *k2, double *P2, double *errP2, int NeffP2, double *k4, double *P4, double *errP4, int NeffP4, double *k11, double *k22, double *k33, double *B0, double *errB0, double *Bnoise, int NeffB0, double *k0SGC, double *P0SGC, double *errP0SGC,double PnoiseSGC,int NeffP0SGC, double *k2SGC, double *P2SGC, double *errP2SGC,int NeffP2SGC, double *k4SGC, double *P4SGC, double *errP4SGC,int NeffP4SGC, double *k11SGC, double *k22SGC, double *k33SGC, double *B0SGC,double *errB0SGC, double *BnoiseSGC,int NeffB0SGC, double *cov, double *covSGC, double alpha_min, double alpha_max, char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free,char *fog_bs ,int Nchunks, char *path_output, char *identifier, char *do_plot, char *use_prop_cov, char *path_to_cov, long int Nsteps, char *do_power_spectrum, char *do_bispectrum, char *spacing_dataNGC,char *spacing_dataSGC, char *spacing_theory,char *type_of_analysis_BAO,char *type_of_analysis_FS)
{
double fraction;
int trial_mcmc;
long int N_max,N_print,N_burnout,j;
double **vector_buffer;
FILE *f;
double *Cov_prop;
double *vector_mean;
int N_Cov_prop;
gsl_matrix *transform_inverse;
gsl_matrix *transform;
char name_file[2000];
long int *params_mcmc;
int factor_sampling_mask;
int i_thread;

factor_sampling_mask=10;//Sampling factor boost wrt to Neff max when the mask is applied. How much do I need to sample my model pre-mask-apply in order to have a satisfactory mask response? 10 times seems reasonable.

if (nthreads<2){
  sprintf(name_file,"%s/mcmcFS_output_%s.txt",path_output,identifier);
  f=fopen(name_file,"w");
  if(f==NULL){printf("Error writing %s. Exiting now...\n",name_file);exit(0);}
  fclose(f);
} else {
  for (i_thread=0;i_thread<nthreads;i_thread++){
    sprintf(name_file,"%s/mcmcFS_output_%s__%d.txt",path_output,identifier,i_thread+1);
    f=fopen(name_file,"w");  
    if(f==NULL){printf("Error writing %s. Exiting now...\n",name_file);exit(0);}
    fclose(f);
  }
  sprintf(name_file,"%s/mcmcFS_output_%s.txt",path_output,identifier);
}
//Number of free parameters in each case. 
N_Cov_prop=3;//b1, b2, A

if(strcmp(local_b2s2, "no") == 0){N_Cov_prop++;}
if(strcmp(local_b3nl, "no") == 0){N_Cov_prop++;}
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){N_Cov_prop++;}
if(strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 ){N_Cov_prop++;}
N_Cov_prop=N_Cov_prop*Nchunks;//x2 
if(strcmp(sigma8_free, "yes") == 0){N_Cov_prop++;}
if(strcmp(RSD_fit, "yes") == 0){N_Cov_prop=N_Cov_prop+3;}//f, alpha-para, alpha-perp


if(strcmp(use_prop_cov, "yes") == 0)//full mcmc run with proposal
{
trial_mcmc=0;//no trial
fraction=1;
Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

//read prop. covariance
//read_prop_cov(vector_buffer,N_max,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);
read_prop_cov(NULL,0,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);
transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);

generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse

free(Cov_prop);
//mcmc_kernel();

//mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,NULL,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,NULL,NULL, N_Plin, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0, P0, errP0, NeffP0, k2, P2, errP2, NeffP2, k4, P4, errP4, NeffP4, k11, k22, k33, B0, errB0, Bnoise, NeffB0, k0SGC, P0SGC, errP0SGC, NeffP0SGC, k2SGC, P2SGC, errP2SGC, NeffP2SGC, k4SGC, P4SGC, errP4SGC, NeffP4SGC, k11SGC, k22SGC, k33SGC, B0SGC, BnoiseSGC, NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, 0,0, 0,Theory, Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , 0,factor_sampling_mask, spacing_dataNGC,spacing_dataSGC,spacing_theory);
mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,NULL,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, NULL,fit_BAO,NULL,NULL, 0, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  NULL,k0, NULL,P0, NULL,errP0, 0,NeffP0, NULL,k2, NULL,P2, NULL,errP2, 0,NeffP2, NULL,k4, NULL,P4, NULL,errP4, 0,NeffP4, NULL,k11, NULL,k22, NULL,k33, NULL,B0, NULL,errB0, NULL,Bnoise, 0,NeffB0, NULL,k0SGC, NULL,P0SGC, NULL,errP0SGC, 0,NeffP0SGC, NULL,k2SGC, NULL,P2SGC, NULL,errP2SGC, 0,NeffP2SGC, NULL,k4SGC, NULL,P4SGC, NULL,errP4SGC, 0,NeffP4SGC, NULL,k11SGC, NULL,k22SGC, NULL,k33SGC, NULL,B0SGC, NULL,BnoiseSGC, 0,NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, 0,0, 0,Theory,N_Plin,Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , 0,factor_sampling_mask, NULL,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS);



}
else//two iterations: first trial to get a proposal covariance, seccond full run with such proposal
{
trial_mcmc=1;//trial iteration
fraction=0.10;
params_mcmc=(long int*) calloc( 2, sizeof(long int ));
set_mcmc_parameters(params_mcmc);
N_burnout=params_mcmc[1];
free(params_mcmc);

N_max=(long int)(Nsteps*fraction);
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

//use proposed mean
set_proposal_mean(vector_mean,N_Cov_prop);

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);

generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse

vector_buffer = (double**) calloc(N_Cov_prop+2,sizeof(double*));

for(j=0;j<N_Cov_prop+2;j++)
{
   vector_buffer[j] = (double*) calloc(N_max,sizeof(double));
}
//mcmc_kernel();

//mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,NULL,NULL, N_Plin, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0, P0, errP0, NeffP0, k2, P2, errP2, NeffP2, k4, P4, errP4, NeffP4, k11, k22, k33, B0, errB0, Bnoise, NeffB0, k0SGC, P0SGC, errP0SGC, NeffP0SGC, k2SGC, P2SGC, errP2SGC, NeffP2SGC, k4SGC, P4SGC, errP4SGC, NeffP4SGC, k11SGC, k22SGC, k33SGC, B0SGC, BnoiseSGC, NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, 0,0, 0,Theory, Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs, 0,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory );
mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, NULL,fit_BAO,NULL,NULL, 0, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  NULL,k0, NULL,P0, NULL,errP0, 0,NeffP0, NULL,k2, NULL,P2, NULL,errP2, 0,NeffP2, NULL,k4, NULL,P4, NULL,errP4, 0,NeffP4, NULL,k11, NULL,k22, NULL,k33, NULL,B0, NULL,errB0, NULL,Bnoise, 0,NeffB0, NULL,k0SGC, NULL,P0SGC, NULL,errP0SGC, 0,NeffP0SGC, NULL,k2SGC, NULL,P2SGC, NULL,errP2SGC, 0,NeffP2SGC, NULL,k4SGC, NULL,P4SGC, NULL,errP4SGC, 0,NeffP4SGC, NULL,k11SGC, NULL,k22SGC, NULL,k33SGC, NULL,B0SGC, NULL,BnoiseSGC, 0,NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, 0,0, 0,Theory,N_Plin,Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , 0,factor_sampling_mask, NULL,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS);

Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

//get covariance from buffer
read_prop_cov(vector_buffer,N_max,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);

freeTokens(vector_buffer,N_Cov_prop+2);

trial_mcmc=0;//no trial
fraction=1;

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse

free(Cov_prop);

//mcmc_kernel();
//printf("%d %d %d\n",NeffP0,NeffP2,NeffP4);
//printf("%d %d %d\n",NeffP0SGC,NeffP2SGC,NeffP4SGC);
//mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,NULL,NULL, N_Plin, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0, P0, errP0, NeffP0, k2, P2, errP2, NeffP2, k4, P4, errP4, NeffP4, k11, k22, k33, B0, errB0, Bnoise, NeffB0, k0SGC, P0SGC, errP0SGC, NeffP0SGC, k2SGC, P2SGC, errP2SGC, NeffP2SGC, k4SGC, P4SGC, errP4SGC, NeffP4SGC, k11SGC, k22SGC, k33SGC, B0SGC, BnoiseSGC, NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, 0,0, 0,Theory, Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs, 0,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory );
mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, NULL,fit_BAO,NULL,NULL, 0, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  NULL,k0, NULL,P0, NULL,errP0, 0,NeffP0, NULL,k2, NULL,P2, NULL,errP2, 0,NeffP2, NULL,k4, NULL,P4, NULL,errP4, 0,NeffP4, NULL,k11, NULL,k22, NULL,k33, NULL,B0, NULL,errB0, NULL,Bnoise, 0,NeffB0, NULL,k0SGC, NULL,P0SGC, NULL,errP0SGC, 0,NeffP0SGC, NULL,k2SGC, NULL,P2SGC, NULL,errP2SGC, 0,NeffP2SGC, NULL,k4SGC, NULL,P4SGC, NULL,errP4SGC, 0,NeffP4SGC, NULL,k11SGC, NULL,k22SGC, NULL,k33SGC, NULL,B0SGC, NULL,BnoiseSGC, 0,NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, 0,0, 0,Theory,N_Plin,Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , 0,factor_sampling_mask, NULL,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS);

}





}



void make_a_rsd_plot(char *type_of_analysis,double *parameters2,double chi2_min, double **Theory,int N_Plin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC,double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0, double *P0, double *errP0, double Pnoise, double *k2, double *P2, double *errP2, double *k4, double *P4, double *errP4, double *k0SGC, double *P0SGC, double *errP0SGC, double PnoiseSGC, double *k2SGC, double *P2SGC, double *errP2SGC, double *k4SGC, double *P4SGC, double *errP4SGC, int NeffP0, int NeffP2, int NeffP4,int NeffP0SGC, int NeffP2SGC, int NeffP4SGC, double *cov, double *covSGC,  char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free, char *fog_bs, int Nchunks,int Npoints, int Ndof, char *path_output, char *identifier, fftw_plan plan1, fftw_plan plan2, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum, int factor_sampling_mask_in, char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory)
{
FILE *f1,*f2,*f12;
int open,NeffP;
char name1[2000],name2[2000],name12[2000];
char nameP0_1[2000],nameP0_2[2000],nameP0_12[2000];
char nameP2_1[2000],nameP2_2[2000],nameP2_12[2000];
char nameP4_1[2000],nameP4_2[2000],nameP4_12[2000];

int Neffmax,NeffmaxSGC;

double ptheo,ptheoSGC,ptheotot;
double P0tot,errP0tot;

double *parameters1,*ktheo,*Ptheo0,*Ptheo2,*Ptheo4;
double *ktheo0,*ktheo2,*ktheo4;
double *ktheoSGC,*Ptheo0SGC,*Ptheo2SGC,*Ptheo4SGC;
double *ktheo0SGC,*ktheo2SGC,*ktheo4SGC;
int i,j,l,i1;
int dimension;
int modeP0,modeP2,modeP4;
int Ntheo=N_Plin;
int offset;
double difference,difference_min;
int ik01_P0,ik01_P2,ik01_P4;
int combine;
int factor_sampling_mask;
int interpolation_order,shiftN;
double w1,w2,w0;
int Ninterpol;

interpolation_order=1;

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}


//for(i=0;i<11;i++){printf("%d %lf\n",i,parameters2[i]);}

factor_sampling_mask=factor_sampling_mask_in;

combine=0;
dimension=11;
offset=3;
if(strcmp(local_b2s2, "no") == 0){offset++;}//b2s2
if(strcmp(local_b3nl, "no") == 0){offset++;}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){offset++;}//fog_ps
if(strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 ){offset++;}//fog_bs
//if(strcmp(sigma8_free, "yes") == 0){offset++;}//s8 as free parameter
//if(strcmp(RSD_fit, "yes") == 0){offset=offset+3;}//f,apara,aperp as free parameter


modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

/*
ktheo = (double*) calloc( Ntheo, sizeof(double));

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
Ptheo0 = (double*) calloc( Ntheo, sizeof(double));
}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
Ptheo2 = (double*) calloc( Ntheo, sizeof(double));
}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
Ptheo4 = (double*) calloc( Ntheo, sizeof(double));
}
*/
if(strcmp(path_to_mask1, "none") == 0 )//no mask
{
if(NeffP0>0){ktheo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2>0){ktheo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4>0){ktheo4 = (double*) calloc( NeffP4, sizeof(double));}
if(NeffP0>0){Ptheo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2>0){Ptheo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4>0){Ptheo4 = (double*) calloc( NeffP4, sizeof(double));}
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
if(modeP0==1){Ptheo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2==1){Ptheo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4==1){Ptheo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}


//ktheo = (double*) calloc( Ntheo, sizeof(double));
//if(modeP0==1){Ptheo0 = (double*) calloc( Ntheo, sizeof(double));}
//if(modeP2==1){Ptheo2 = (double*) calloc( Ntheo, sizeof(double));}
//if(modeP4==1){Ptheo4 = (double*) calloc( Ntheo, sizeof(double));}
}


parameters1 = (double*) calloc( dimension, sizeof(double));
i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
parameters1[i]=parameters2[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<3){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}
}

if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==3 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=4 && i<=6){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==7 && strcmp(local_b2s2, "yes") == 0){
parameters1[i]=-4./7.*(parameters1[i-3]-1);
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==8 && strcmp(local_b3nl, "yes") == 0){
parameters1[i]=32./315.*(parameters1[i-4]-1);
}
if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==9 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=0;
}

if(i==9 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){
parameters1[i]=parameters1[i-1];
}

if(i==10 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){
parameters1[i]=0;
}


}

//for(i=0;i<dimension;i++){printf("NGC %d %lf\n",i,parameters1[i]);}

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, ktheo,ktheo0,ktheo2,ktheo4, Ptheo0, Ptheo2, Ptheo4,Pnoise,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC,k0,k2,k4, path_to_mask1,plan1, plan2, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1],spacing_dataNGC,spacing_theory);

free(parameters1);


if(Nchunks==2)
{

factor_sampling_mask=factor_sampling_mask_in;

	parameters1 = (double*) calloc( dimension, sizeof(double));

	if(strcmp(path_to_mask2, "none") == 0 )//no mask
	{
		if(NeffP0SGC>0){ktheo0SGC = (double*) calloc( NeffP0SGC, sizeof(double));}
		if(NeffP2SGC>0){ktheo2SGC = (double*) calloc( NeffP2SGC, sizeof(double));}
		if(NeffP4SGC>0){ktheo4SGC = (double*) calloc( NeffP4SGC, sizeof(double));}
		if(NeffP0SGC>0){Ptheo0SGC = (double*) calloc( NeffP0SGC, sizeof(double));}
		if(NeffP2SGC>0){Ptheo2SGC = (double*) calloc( NeffP2SGC, sizeof(double));}
		if(NeffP4SGC>0){Ptheo4SGC = (double*) calloc( NeffP4SGC, sizeof(double));}
	}
	else
	{
		//ktheoSGC = (double*) calloc( Ntheo, sizeof(double));
		//if(modeP0==1){Ptheo0SGC = (double*) calloc( Ntheo, sizeof(double));}
		//if(modeP2==1){Ptheo2SGC = (double*) calloc( Ntheo, sizeof(double));}
		//if(modeP4==1){Ptheo4SGC = (double*) calloc( Ntheo, sizeof(double));}
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
		if(modeP0==1){Ptheo0SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}
		if(modeP2==1){Ptheo2SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}
		if(modeP4==1){Ptheo4SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}





	}


	/*
	   ktheoSGC = (double*) calloc( Ntheo, sizeof(double));

	   if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
	   {
	   Ptheo0SGC = (double*) calloc( Ntheo, sizeof(double));
	   }
	   if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
	   {
	   Ptheo2SGC = (double*) calloc( Ntheo, sizeof(double));
	   }
	   if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
	   {
	   Ptheo4SGC = (double*) calloc( Ntheo, sizeof(double));
	   }
	   */
	i1=0;
	for(i=0;i<dimension;i++){

		if(strcmp(RSD_fit, "yes") == 0 && i<3){
			parameters1[i]=parameters2[i1];
			i1++;
		}

		if(strcmp(RSD_fit, "no") == 0 && i<3){
			if(i==0){parameters1[i]=1.;}
			if(i==1){parameters1[i]=1.;}
			if(i==2){parameters1[i]=0;}
		}

		if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
			parameters1[i]=parameters2[i1];
			i1++;
		}

		if(i==3 && strcmp(sigma8_free, "no") == 0 ){
			parameters1[i]=Theory[0][41];
		}

		if(i>=4 && i<=6){
			parameters1[i]=parameters2[i1+offset];
			i1++;
		}

		if(i==7 && strcmp(local_b2s2, "no") == 0){
			parameters1[i]=parameters2[i1+offset];
			i1++;
		}

		if(i==7 && strcmp(local_b2s2, "yes") == 0){
			parameters1[i]=-4./7.*(parameters1[i-3]-1);
		}

		if(i==8 && strcmp(local_b3nl, "no") == 0){
			parameters1[i]=parameters2[i1+offset];
			i1++;
		}

		if(i==8 && strcmp(local_b3nl, "yes") == 0){
			parameters1[i]=32./315.*(parameters1[i-4]-1);
		}
		if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
			parameters1[i]=parameters2[i1+offset];
			i1++;
		}

		if(i==9 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){
			parameters1[i]=0;
		}

		if(i==9 && strcmp(do_power_spectrum, "no") == 0){
			parameters1[i]=0;
		}

		if(i==10 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){
parameters1[i]=parameters2[i1+offset];
i1++;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){
parameters1[i]=parameters1[i-1];
}

if(i==10 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==10 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){
parameters1[i]=0;
}

}


//for(i=0;i<dimension;i++){printf("SGC %d %lf\n",i,parameters1[i]);}

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, ktheoSGC,ktheo0SGC,ktheo2SGC,ktheo4SGC, Ptheo0SGC, Ptheo2SGC, Ptheo4SGC,PnoiseSGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC,k0SGC,k2SGC,k4SGC, path_to_mask2,plan1, plan2, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4[NeffP4SGC-1],spacing_dataSGC,spacing_theory);

free(parameters1);



combine=0;
//Decide if combine or not here.if combine==0 then combine, otherwise no
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


sprintf(nameP0_1,"%s/Monopole1_FS_%s.txt",path_output,identifier);
sprintf(nameP0_2,"%s/Monopole2_FS_%s.txt",path_output,identifier);
sprintf(nameP0_12,"%s/Monopole12_FS_%s.txt",path_output,identifier);

sprintf(nameP2_1,"%s/Quadrupole1_FS_%s.txt",path_output,identifier);
sprintf(nameP2_2,"%s/Quadrupole2_FS_%s.txt",path_output,identifier);
sprintf(nameP2_12,"%s/Quadrupole12_FS_%s.txt",path_output,identifier);

sprintf(nameP4_1,"%s/Hexadecapole1_FS_%s.txt",path_output,identifier);
sprintf(nameP4_2,"%s/Hexadecapole2_FS_%s.txt",path_output,identifier);
sprintf(nameP4_12,"%s/Hexadecapole12_FS_%s.txt",path_output,identifier);


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

fprintf(f1,"#k \t Pdata \t err \t Pmodel\n");

if(Nchunks==2)
{
fprintf(f2,"#k \t Pdata \t err \t Pmodel\n");
if(combine==0){fprintf(f12,"#k \t Pdata \t err \t Pmodel\n");}
}

if(l==0){NeffP=NeffP0;}
if(l==2){NeffP=NeffP2;}
if(l==4){NeffP=NeffP4;}

for(i=0;i<NeffP;i++)
{

if(l==0)
{
if(strcmp(path_to_mask1, "none") == 0 ){ptheo=Ptheo0[i]-Pnoise;}
else{

Ninterpol=determine_N_singlearray(ktheo,k0[i],Neffmax*factor_sampling_mask,spacing_dataNGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k0[i]<=0)
{
ptheo=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(ktheo,k0[i],Ninterpol,spacing_dataNGC);
}
ptheo=P_interpol_fast(k0[i],Ptheo0,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2)-Pnoise;
}

//ptheo=P_interpol(k0[i],ktheo,Ptheo0,Neffmax*factor_sampling_mask)-Pnoise;

}
fprintf(f1,"%e %e %e %e\n",k0[i],P0[i],errP0[i],ptheo);
}

if(l==2)
{
if(strcmp(path_to_mask1, "none") == 0 ){ptheo=Ptheo2[i];}
else{

Ninterpol=determine_N_singlearray(ktheo,k2[i],Neffmax*factor_sampling_mask,spacing_dataNGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k2[i]<=0)
{
ptheo=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(ktheo,k2[i],Ninterpol,spacing_dataNGC);
}
ptheo=P_interpol_fast(k2[i],Ptheo2,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
}


//ptheo=P_interpol(k2[i],ktheo,Ptheo2,Neffmax*factor_sampling_mask);

}

fprintf(f1,"%e %e %e %e\n",k2[i],P2[i],errP2[i],ptheo);
}

if(l==4)
{
if(strcmp(path_to_mask1, "none") == 0 ){ptheo=Ptheo4[i];}
else{

Ninterpol=determine_N_singlearray(ktheo,k4[i],Neffmax*factor_sampling_mask,spacing_dataNGC);
if(Ninterpol>=Neffmax*factor_sampling_mask-shiftN || Ninterpol<0  || k4[i]<=0)
{
ptheo=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);
w1=determine_w1_2ndorder_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);
w2=determine_w2_2ndorder_singlearray(ktheo,k4[i],Ninterpol,spacing_dataNGC);
}
ptheo=P_interpol_fast(k4[i],Ptheo4,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2);
}


//ptheo=P_interpol(k4[i],ktheo,Ptheo4,Neffmax*factor_sampling_mask);

}
fprintf(f1,"%e %e %e %e\n",k4[i],P4[i],errP4[i],ptheo);
}



if(Nchunks==2)
{
//This part is written after, separated from this loop in case Neff is not NeffSGC
if(combine==0)
{


if(l==0)
{
if(strcmp(path_to_mask2, "none") == 0){ptheoSGC=Ptheo0SGC[i]-PnoiseSGC;}
else{

Ninterpol=determine_N_singlearray(ktheoSGC,k0SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k0SGC[i]<=0)
{
ptheoSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k0SGC[i],Ptheo0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2)-PnoiseSGC;
}


//ptheoSGC=P_interpol(k0SGC[i],ktheoSGC,Ptheo0SGC,NeffmaxSGC*factor_sampling_mask)-PnoiseSGC;

}
}

if(l==2)
{
if(strcmp(path_to_mask2, "none") == 0){ptheoSGC=Ptheo2SGC[i];}
else{

Ninterpol=determine_N_singlearray(ktheoSGC,k2SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k2SGC[i]<=0)
{
ptheoSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k2SGC[i],Ptheo2SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
}


//ptheoSGC=P_interpol(k2SGC[i],ktheoSGC,Ptheo2SGC,NeffmaxSGC*factor_sampling_mask);

}
}

if(l==4)
{
if(strcmp(path_to_mask2, "none") == 0){ptheoSGC=Ptheo4SGC[i];}
else{


Ninterpol=determine_N_singlearray(ktheo,k4SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k4SGC[i]<=0)
{
ptheoSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k4SGC[i],Ptheo4SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
}

//ptheoSGC=P_interpol(k4SGC[i],ktheoSGC,Ptheo4SGC,NeffmaxSGC*factor_sampling_mask);

}
}

if(l==0)
{

P0tot=(P0[i]*pow(errP0[ik01_P0],-2)+P0SGC[i]*pow(errP0SGC[ik01_P0],-2))/( pow(errP0[ik01_P0],-2)+pow(errP0SGC[ik01_P0],-2) );

ptheotot=(ptheo*pow(errP0[ik01_P0],-2)+ptheoSGC*pow(errP0SGC[ik01_P0],-2))/( pow(errP0[ik01_P0],-2)+pow(errP0SGC[ik01_P0],-2) );

errP0tot=errP0[i]*sqrt( pow(errP0[ik01_P0],-2)/(pow(errP0[ik01_P0],-2)+pow(errP0SGC[ik01_P0],-2)) );//errNS=errNGC*sqrt(ANGC/(ANGC+ASGC))

fprintf(f12,"%e %e %e %e\n",k0[i],P0tot,errP0tot,ptheotot);
}

if(l==2)
{
P0tot=(P2[i]*pow(errP2[ik01_P2],-2)+P2SGC[i]*pow(errP2SGC[ik01_P2],-2))/( pow(errP2[ik01_P2],-2)+pow(errP2SGC[ik01_P2],-2) );

ptheotot=(ptheo*pow(errP2[ik01_P2],-2)+ptheoSGC*pow(errP2SGC[ik01_P2],-2))/( pow(errP2[ik01_P2],-2)+pow(errP2SGC[ik01_P2],-2) );

errP0tot=errP2[i]*sqrt( pow(errP2[ik01_P2],-2)/(pow(errP2[ik01_P2],-2)+pow(errP2SGC[ik01_P2],-2)) );//errNS=errNGC*sqrt(ANGC/(ANGC+ASGC))

fprintf(f12,"%e %e %e %e\n",k2[i],P0tot,errP0tot,ptheotot);

}

if(l==4)
{
P0tot=(P4[i]*pow(errP4[ik01_P4],-2)+P4SGC[i]*pow(errP4SGC[ik01_P4],-2))/( pow(errP4[ik01_P4],-2)+pow(errP4SGC[ik01_P4],-2) );

ptheotot=(ptheo*pow(errP4[ik01_P4],-2)+ptheoSGC*pow(errP4SGC[ik01_P4],-2))/( pow(errP4[ik01_P4],-2)+pow(errP4SGC[ik01_P4],-2) );

errP0tot=errP4[i]*sqrt( pow(errP4[ik01_P4],-2)/(pow(errP4[ik01_P4],-2)+pow(errP4SGC[ik01_P4],-2)) );//errNS=errNGC*sqrt(ANGC/(ANGC+ASGC))

fprintf(f12,"%e %e %e %e\n",k4[i],P0tot,errP0tot,ptheotot);

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
if(strcmp(path_to_mask2, "none") != 0 ){

Ninterpol=determine_N_singlearray(ktheoSGC,k0SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k0SGC[i]<=0)
{
ptheoSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k0SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k0SGC[i],Ptheo0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2)-PnoiseSGC;
}

//ptheoSGC=P_interpol(k0SGC[i],ktheoSGC,Ptheo0SGC,NeffmaxSGC*factor_sampling_mask)-PnoiseSGC;

}
else{ptheoSGC=Ptheo0SGC[i]-PnoiseSGC;}
fprintf(f2,"%e %e %e %e\n",k0SGC[i],P0SGC[i],errP0SGC[i],ptheoSGC);
}

if(l==2)
{
if(strcmp(path_to_mask2, "none") != 0 ){

Ninterpol=determine_N_singlearray(ktheoSGC,k2SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k2SGC[i]<=0)
{
ptheoSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k2SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k2SGC[i],Ptheo2SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
}


//ptheoSGC=P_interpol(k2SGC[i],ktheoSGC,Ptheo2SGC,NeffmaxSGC*factor_sampling_mask);

}
else{ptheoSGC=Ptheo2SGC[i];}
fprintf(f2,"%e %e %e %e\n",k2SGC[i],P2SGC[i],errP2SGC[i],ptheoSGC);
}

if(l==4)
{
if(strcmp(path_to_mask2, "none") != 0 ){

Ninterpol=determine_N_singlearray(ktheoSGC,k4SGC[i],NeffmaxSGC*factor_sampling_mask,spacing_dataSGC);
if(Ninterpol>=NeffmaxSGC*factor_sampling_mask-shiftN || Ninterpol<0  || k4SGC[i]<=0)
{
ptheoSGC=0;
}
else{
if(interpolation_order==1)
{
w1=determine_w1_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
}
if(interpolation_order==2)
{
w0=determine_w0_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w1=determine_w1_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
w2=determine_w2_2ndorder_singlearray(ktheoSGC,k4SGC[i],Ninterpol,spacing_dataSGC);
}
ptheoSGC=P_interpol_fast(k4SGC[i],Ptheo4SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2);
}


//ptheoSGC=P_interpol(k4SGC[i],ktheoSGC,Ptheo4SGC,NeffmaxSGC*factor_sampling_mask);

}
else{ptheoSGC=Ptheo4SGC[i];}
fprintf(f2,"%e %e %e %e\n",k4SGC[i],P4SGC[i],errP4SGC[i],ptheoSGC);
}


}//for NeffP

}//if Nchunks=2

//write there the header and best-fitting parameters
fprintf(f1,"#");
if(strcmp(RSD_fit, "yes") == 0 && strcmp(sigma8_free, "yes") == 0){fprintf(f1,"apara \t aperp \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "yes") == 0 && strcmp(sigma8_free, "no") == 0){fprintf(f1,"apara \t aperp \t fs8 \t");}//f,apara,aperp as free parameter
if(strcmp(sigma8_free, "yes") == 0){fprintf(f1,"s8 \t");}
if(strcmp(sigma8_free, "yes") == 0){fprintf(f1,"b1_NGC \t b2_NGC \t Anoise_NGC\t");}
if(strcmp(sigma8_free, "no") == 0){fprintf(f1,"b1s8_NGC \t b2s8_NGC \t Anoise_NGC\t");}
if(strcmp(local_b2s2, "no") == 0){fprintf(f1,"b2s2_NGC \t");}//b2s2
if(strcmp(local_b3nl, "no") == 0){fprintf(f1,"b3nl_NGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f1,"sigmaP_NGC\t");}
fprintf(f1,"chi2 \t Npoints-Ndof\n");

//print here values on f1
i1=0;
fprintf(f1,"#");
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
if( strcmp(sigma8_free, "yes") == 0){fprintf(f1,"%lf\t",parameters2[i1]);}
if( strcmp(sigma8_free, "no") == 0 && i<2){fprintf(f1,"%lf\t",parameters2[i1]);}
if( strcmp(sigma8_free, "no") == 0 && i==2){fprintf(f1,"%lf\t",parameters2[i1]*Theory[0][41]);}
i1++;
}


if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(i>=4 && i<=6){
if( strcmp(sigma8_free, "no") == 0 && i<6){fprintf(f1,"%lf\t",parameters2[i1]*Theory[0][41]);}
if( strcmp(sigma8_free, "no") == 0 && i==6){fprintf(f1,"%lf\t",parameters2[i1]);}
if( strcmp(sigma8_free, "yes") == 0){fprintf(f1,"%lf\t",parameters2[i1]);}
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

}


fprintf(f1,"%lf \t %d-%d\n",chi2_min,Npoints,Ndof);

if(Nchunks==2)
{

fprintf(f2,"#");
if(strcmp(RSD_fit, "yes") == 0 && strcmp(sigma8_free, "yes") == 0){fprintf(f2,"apara \t aperp \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "yes") == 0 && strcmp(sigma8_free, "no") == 0){fprintf(f2,"apara \t aperp \t fs8 \t");}//f,apara,aperp as free parameter
if(strcmp(sigma8_free, "yes") == 0){fprintf(f2,"s8 \t");}
if(strcmp(sigma8_free, "yes") == 0){fprintf(f2,"b1_SGC \t b2_SGC \t Anoise_SGC\t");}
if(strcmp(sigma8_free, "no") == 0){fprintf(f2,"b1s8_SGC \t b2s8_SGC \t Anoise_SGC\t");}
if(strcmp(local_b2s2, "no") == 0){fprintf(f2,"b2s2_SGC \t");}//b2s2
if(strcmp(local_b3nl, "no") == 0){fprintf(f2,"b3nl_SGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f2,"sigmaP_SGC\t");}
fprintf(f2,"chi2 \t Npoints-Ndof\n");



i1=0;
fprintf(f2,"#");
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
if( strcmp(sigma8_free, "yes") == 0){fprintf(f2,"%lf\t",parameters2[i1]);}
if( strcmp(sigma8_free, "no") == 0 && i<2){fprintf(f2,"%lf\t",parameters2[i1]);}
if( strcmp(sigma8_free, "no") == 0 && i==2){fprintf(f2,"%lf\t",parameters2[i1]*Theory[0][41]);}
i1++;
}


if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
fprintf(f2,"%lf\t",parameters2[i1]);
i1++;
}

if(i>=4 && i<=6){
if( strcmp(sigma8_free, "no") == 0 && i<6){fprintf(f2,"%lf\t",parameters2[i1+offset]*Theory[0][41]);}
if( strcmp(sigma8_free, "no") == 0 && i==6){fprintf(f2,"%lf\t",parameters2[i1+offset]);}
if( strcmp(sigma8_free, "yes") == 0){fprintf(f2,"%lf\t",parameters2[i1+offset]);}
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
fprintf(f2,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
fprintf(f2,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
fprintf(f2,"%lf\t",parameters2[i1+offset]);
i1++;
}

}
fprintf(f2,"%lf \t %d-%d\n",chi2_min,Npoints,Ndof);




if(combine==0){

fprintf(f12,"#");
if(strcmp(RSD_fit, "yes") == 0 && strcmp(sigma8_free, "yes") == 0){fprintf(f12,"apara \t aperp \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "yes") == 0 && strcmp(sigma8_free, "no") == 0){fprintf(f12,"apara \t aperp \t fs8 \t");}//f,apara,aperp as free parameter
if(strcmp(sigma8_free, "yes") == 0){fprintf(f12,"s8 \t");}
if(strcmp(sigma8_free, "yes") == 0){fprintf(f12,"b1_NGC \t b2_NGC \t Anoise_NGC\t");}
if(strcmp(sigma8_free, "no") == 0){fprintf(f12,"b1s8_NGC \t b2s8_NGC \t Anoise_NGC\t");}
if(strcmp(local_b2s2, "no") == 0){fprintf(f12,"b2s2_NGC \t");}//b2s2
if(strcmp(local_b3nl, "no") == 0){fprintf(f12,"b3nl_NGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f12,"sigmaP_SGC\t");}

if(strcmp(sigma8_free, "yes") == 0){fprintf(f12,"b1_SGC \t b2_SGC \t Anoise_SGC\t");}
if(strcmp(sigma8_free, "no") == 0){fprintf(f12,"b1s8_SGC \t b2s8_SGC \t Anoise_SGC\t");}
if(strcmp(local_b2s2, "no") == 0){fprintf(f12,"b2s2_SGC \t");}//b2s2
if(strcmp(local_b3nl, "no") == 0){fprintf(f12,"b3nl_SGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f12,"sigmaP_SGC\t");}
fprintf(f12,"chi2 \t Npoints-Ndof\n");

i1=0;
fprintf(f12,"#");
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
if( strcmp(sigma8_free, "yes") == 0){fprintf(f12,"%lf\t",parameters2[i1]);}
if( strcmp(sigma8_free, "no") == 0 && i<2){fprintf(f12,"%lf\t",parameters2[i1]);}
if( strcmp(sigma8_free, "no") == 0 && i==2){fprintf(f12,"%lf\t",parameters2[i1]*Theory[0][41]);}
i1++;
}


if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if(i>=4 && i<=6){
if( strcmp(sigma8_free, "no") == 0 && i<6){fprintf(f12,"%lf\t",parameters2[i1+offset]*Theory[0][41]);}
if( strcmp(sigma8_free, "no") == 0 && i==6){fprintf(f12,"%lf\t",parameters2[i1+offset]);}
if( strcmp(sigma8_free, "yes") == 0){fprintf(f12,"%lf\t",parameters2[i1+offset]);}
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

}

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
//if( strcmp(sigma8_free, "yes") == 0){fprintf(f12,"%lf\t",parameters2[i1]);}
//if( strcmp(sigma8_free, "no") == 0 && i<2){fprintf(f12,"%lf\t",parameters2[i1]);}
//if( strcmp(sigma8_free, "no") == 0 && i==2){fprintf(f12,"%lf\t",parameters2[i1]*Theory[0][41]);}
i1++;
}


if(i==3 && strcmp(sigma8_free, "yes") == 0 ){
//fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if(i>=4 && i<=6){
if( strcmp(sigma8_free, "no") == 0 && i<6){fprintf(f12,"%lf\t",parameters2[i1+offset]*Theory[0][41]);}
if( strcmp(sigma8_free, "no") == 0 && i==6){fprintf(f12,"%lf\t",parameters2[i1+offset]);}
if( strcmp(sigma8_free, "yes") == 0){fprintf(f12,"%lf\t",parameters2[i1+offset]);}
i1++;
}

if(i==7 && strcmp(local_b2s2, "no") == 0){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==8 && strcmp(local_b3nl, "no") == 0){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==9 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

}


fprintf(f12,"%lf \t %d-%d\n",chi2_min,Npoints,Ndof);




}

}

fclose(f1);
if(Nchunks==2)
{
fclose(f2);
if(combine==0){fclose(f12);}
}


}//open

}//for l=0,2,4


if(modeP0==1){
if(strcmp(path_to_mask1, "none") == 0 ){free(ktheo0);}
free(Ptheo0);}
if(modeP2==1){
if(strcmp(path_to_mask1, "none") == 0 ){free(ktheo2);}
free(Ptheo2);}
if(modeP4==1){
if(strcmp(path_to_mask1, "none") == 0 ){free(ktheo4);}
free(Ptheo4);}

if(strcmp(path_to_mask1, "none") != 0 ){
free(ktheo);
}

if(Nchunks==2)
{
if(modeP0==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo0SGC);}
free(Ptheo0SGC);}
if(modeP2==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo2SGC);}
free(Ptheo2SGC);}
if(modeP4==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo4SGC);}
free(Ptheo4SGC);}
if(strcmp(path_to_mask2, "none") != 0 ){
free(ktheoSGC);
}
}


}
