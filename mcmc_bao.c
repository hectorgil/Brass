#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>
#include "functions.h"
#include "bao.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include "fftlog.h"
#include "priors.h"
#include "rsd.h"
#include "structures.h"
#include "integrals_rsd.h"
#include "cubature.h"
#define Pi (4.*atan(1.))

void set_parameters_convergence(double parameters_convergence[])
{
  double R_convergence=1.005;
  double steps=5e4;
  double chunks=4;

  parameters_convergence[0]=R_convergence;
  parameters_convergence[1]=steps;
  parameters_convergence[2]=chunks;
  if(parameters_convergence[2]<3){printf("Warning, chunks parameter must be larger than 2: chunks=%lf. Exiting now...\n",chunks);exit(0);}
}

void set_parameters_convergence_parallel(int nthreads, double parameters_convergence[])
{
  double R_convergence=1.005;
  double steps=5e4;
  double chunks=nthreads;

  parameters_convergence[0]=R_convergence;
  parameters_convergence[1]=steps;
  parameters_convergence[2]=chunks;
  if(parameters_convergence[2]<3){printf("Warning, chunks parameter must be larger than 2: chunks=%lf. You should increase the number of threads. Exiting now...\n",chunks);exit(0);}
}

void set_mcmc_parameters(long int params[])
{
long int N_print=10000;//buffer to print in file AND to check for convergence (which will require to read all the writting file)
long int Nburnin=10000;//number of firsts steps not considered
//long int N_print=1000;//buffer to print in file AND to check for convergence (which will require to read all the writting file)
//long int Nburnin=0;//1000;//number of firsts steps not considered

params[0]=N_print;
params[1]=Nburnin;
}


void do_log_file(int nthreads, char *path_output, char *name_file,int Nparams,long int lines,double time_run,long int abs_counter,long int j_run,long int Nburnout,char *identifier, double mean_params[],double min_params[])
{
//CHECK THIS FUNCTION FOR NEW TEMPLATE WITH SIGMAS
FILE *f;
long int i,j,l,l1,l2;
long int min_lines;
int chunks;
char name_out[2000];
double *R_convergence,*W_convergence,*B_convergence,*mean;
long int *Neff;
double ***param,**mean_chunk;
double **param_single;
int **weight;
int *weight_single;
long int eff_lines,Nefftot;
double chi2,chi2_min;
long int i_min,j_min;
double **error;
double Rtot;
double parameters_convergence[3];
error = (double **) calloc(Nparams, sizeof(double*));
for(i=0;i<Nparams;i++)
{
error[i]= (double *) calloc(Nparams, sizeof(double));
}

sprintf(name_out,"%s/logfilemcmc_%s.txt",path_output,identifier);
f=fopen(name_out,"w");
fprintf(f,"Number of steps: %ld\n",j_run-(Nburnout+1)*nthreads);
fprintf(f,"Number of accepted steps: %ld\n",abs_counter-(Nburnout+1)*nthreads);
fprintf(f,"Burn-in number of steps: %ld\n",Nburnout*nthreads);
fprintf(f,"Accepted fraction of steps: %lf\n",(abs_counter-(Nburnout+1)*nthreads)*1./(j_run-(Nburnout+1)*nthreads)*1.);
fclose(f);

set_parameters_convergence(parameters_convergence);
min_lines=(long int)(parameters_convergence[1]);
chunks=(int)(parameters_convergence[2]);

if(lines<min_lines*chunks){

//write there is not enough lines to perform an internal convergence test
f=fopen(name_out,"a");
fprintf(f,"Convergence Status: not enough steps for such test\n");
fclose(f);

weight_single = (int *) calloc(lines, sizeof(int));//absolute mean
param_single =  (double **) calloc(lines, sizeof(double*));
mean =  (double *) calloc(Nparams, sizeof(double));//absolute mean
for(l=0;l<lines;l++)
{
param_single[l]= (double *) calloc(Nparams, sizeof(double));
}

Nefftot=0;
f=fopen(name_file,"r");
if(f==NULL){printf("Error, could not read the file %s, which I just wrote! Exiting now....\n",name_file);exit(0);}
for(i=0;i<lines;i++)
{
      fscanf(f,"%d %lf ",&weight_single[i],&chi2);
      if(i==0){chi2_min=chi2;i_min=i;}
      if(chi2<chi2_min){chi2_min=chi2;i_min=i;}
      Nefftot=Nefftot+weight_single[i];
         for(l=0;l<Nparams;l++)
         {
            if(l==Nparams-1){fscanf(f,"%lf\n",&param_single[i][l]);mean[l]=mean[l]+param_single[i][l]*weight_single[i];}
            else{fscanf(f,"%lf ",&param_single[i][l]);mean[l]=mean[l]+param_single[i][l]*weight_single[i];}
         }

}
fclose(f);

for(l=0;l<Nparams;l++){min_params[l]=param_single[i_min][l];/*printf("log %d %lf\n",l,min_params[l]);*/}min_params[Nparams]=chi2_min;

for(l=0;l<Nparams;l++){mean_params[l]=mean[l]/Nefftot*1.;}

for(i=0;i<lines;i++)
{
 for(l1=0;l1<Nparams;l1++)
  {
        for(l2=l1;l2<Nparams;l2++)
        {
            error[l1][l2]=error[l1][l2]+(mean_params[l1]-param_single[i][l1])*(mean_params[l2]-param_single[i][l2])*weight_single[i]/(Nefftot*1.-1);
        }
  }
}
free(mean);
free(weight_single);
freeTokens(param_single,lines);

f=fopen(name_out,"a");
fprintf(f,"\nMean parameters\n");
for(l=0;l<Nparams;l++)
{
fprintf(f,"A[%ld]= %lf pm %lf\n",l,mean_params[l],sqrt(error[l][l]));
}

fprintf(f,"\nBest-fitting parameters\n");
for(l=0;l<Nparams;l++)
{
fprintf(f,"A[%ld]= %lf\n",l,min_params[l]);
}
fprintf(f,"chi2_min= %lf\n\n",chi2_min);

fprintf(f,"\nCross-Covariance parameters\n");
for(l1=0;l1<Nparams;l1++)
{
   for(l2=0;l2<Nparams;l2++)
   {
if(l2>=l1){fprintf(f,"%lf\t",error[l1][l2]/sqrt(error[l1][l1]*error[l2][l2]));}
else{fprintf(f,"%lf\t",error[l2][l1]/sqrt(error[l1][l1]*error[l2][l2]));}

if(l2==Nparams-1){fprintf(f,"\n");}

   }

}
fclose(f);

}
else
{
eff_lines=(long int)(lines*1./chunks*1.);

R_convergence =  (double *) calloc(Nparams, sizeof(double));
W_convergence =  (double *) calloc(Nparams, sizeof(double));
B_convergence =  (double *) calloc(Nparams, sizeof(double));
mean =  (double *) calloc(Nparams, sizeof(double));//absolute mean
Neff =  (long int *) calloc(chunks, sizeof(long int));

param=(double ***) calloc(chunks, sizeof(double**));
mean_chunk= (double **) calloc(chunks, sizeof(double*));
weight= (int **) calloc(chunks, sizeof(int *));

for(j=0;j<chunks;j++)
{
weight[j]= (int *) calloc(eff_lines, sizeof(int));
mean_chunk[j] = (double *) calloc(Nparams, sizeof(double));
param[j] = (double **) calloc(eff_lines, sizeof(double*));

  for(i=0;i<eff_lines;i++)
  {
      param[j][i] = (double *) calloc(Nparams, sizeof(double));
  }

}

f=fopen(name_file,"r");
if(f==NULL){printf("Error, could not read the file %s, which I just wrote! Exiting now....\n",name_file);exit(0);}
Nefftot=0;
for(j=0;j<chunks;j++)
{
      Neff[j]=0;
      for(i=0;i<eff_lines;i++)  
      {
      fscanf(f,"%d %lf ",&weight[j][i],&chi2);
      if(i==0 && j==0){chi2_min=chi2;i_min=i;j_min=j;}
      if(chi2<chi2_min){chi2_min=chi2;i_min=i;j_min=j;}
      Neff[j]=Neff[j]+weight[j][i];
      Nefftot=Nefftot+weight[j][i];
         for(l=0;l<Nparams;l++)
         {
           fscanf(f,"%lf ",&param[j][i][l]);mean[l]=mean[l]+param[j][i][l]*weight[j][i];mean_chunk[j][l]=mean_chunk[j][l]+param[j][i][l]*weight[j][i];
           if(l==Nparams-1){fscanf(f,"\n");}

         }
      }

}
fclose(f);


      for(l=0;l<Nparams;l++){min_params[l]=param[j_min][i_min][l];}min_params[Nparams]=chi2_min;

for(l=0;l<Nparams;l++)
{
   for(j=0;j<chunks;j++)
   {
     mean_chunk[j][l]=mean_chunk[j][l]/Neff[j]*1.;
   }
   mean[l]=mean[l]/Nefftot*1.;mean_params[l]=mean[l];
   B_convergence[l]=0;
   for(j=0;j<chunks;j++)
   {
     B_convergence[l]=B_convergence[l]+pow(mean_chunk[j][l]-mean[l],2)/(chunks-1.);
   }
}


for(l=0;l<Nparams;l++)
{
   W_convergence[l]=0;
   for(j=0;j<chunks;j++)
   {
      for(i=0;i<eff_lines;i++)
      {
         W_convergence[l]=W_convergence[l]+pow( param[j][i][l]-mean_chunk[j][l] ,2)*weight[j][i];
      }
  }
 W_convergence[l]=W_convergence[l]/((Nefftot-chunks)*1.);
}

f=fopen(name_out,"a");
fprintf(f,"Convergence Status\n");
Rtot=0;
for(l=0;l<Nparams;l++)
{
R_convergence[l]=( ( (Nefftot*1./chunks*1.)*1.-1.)/(Nefftot*1./chunks*1.)*1.*W_convergence[l]+B_convergence[l]*(1. + 1./chunks*1. )   )/(W_convergence[l]);
fprintf(f,"R[%ld]= %lf\t",l,R_convergence[l]);
Rtot=Rtot+R_convergence[l]/Nparams*1.;
}
fprintf(f,"\n<R>= %lf\n\n",Rtot);
fclose(f);

for(j=0;j<chunks;j++)
{

   for(i=0;i<eff_lines;i++)
   {

     for(l1=0;l1<Nparams;l1++)
     {
        for(l2=l1;l2<Nparams;l2++)
        {
            error[l1][l2]=error[l1][l2]+(mean[l1]-param[j][i][l1])*(mean[l2]-param[j][i][l2])*weight[j][i]/(Nefftot*1.-1);
        }   
     }
   }

}


f=fopen(name_out,"a");
fprintf(f,"Mean parameters\n");
for(l=0;l<Nparams;l++)
{
fprintf(f,"A[%ld]= %lf pm %lf\n",l,mean_params[l],sqrt(error[l][l]));
}

fprintf(f,"\nBest-fitting parameters\n");
for(l=0;l<Nparams;l++)
{
fprintf(f,"A[%ld]= %lf\n",l,min_params[l]);
}
fprintf(f,"chi2_min= %lf\n\n",chi2_min);

fprintf(f,"\nCross-Covariance parameters\n");
for(l1=0;l1<Nparams;l1++)
{
   for(l2=0;l2<Nparams;l2++)
   {
if(l2>=l1){fprintf(f,"%lf\t",error[l1][l2]/sqrt(error[l1][l1]*error[l2][l2]));}
else{fprintf(f,"%lf\t",error[l2][l1]/sqrt(error[l1][l1]*error[l2][l2]));}

if(l2==Nparams-1){fprintf(f,"\n");}

   }

}
fclose(f);

free(R_convergence); 
free(W_convergence); 
free(B_convergence);
free(mean);
free(Neff);

freeTokens2(param, chunks, eff_lines);
freeTokens(mean_chunk,chunks);
freeTokensInt(weight,chunks);
}



freeTokens(error,Nparams);
f=fopen(name_out,"a");
if(time_run>1.){

if(time_run>24.)
{
fprintf(f,"\nTime taken for mcmc chain: %lf days\n",time_run/24.);
}
else{
fprintf(f,"\nTime taken for mcmc chain: %lf hours\n",time_run);
}

}
else
{
fprintf(f,"\nTime taken for mcmc chain: %lf minutes\n",time_run*60.);
}
fclose(f);

//for(l=0;l<Nparams;l++){printf("logfin %d %lf\n",l,min_params[l]);}

}


/*
void do_Ptheo_multiple_aniso(char *type_BAO_fit, char *type_of_analysis, char *fit_BAO,int modeP0,int modeP2,int modeP4, double k_theo[],double k_theo0[], double k_theo2[],double k_theo4[], double P_theo0[], double P_theo2[], double P_theo4[],int NeffP0,int NeffP2,int NeffP4,int factor_for_sampling, double *parameters1,double *k_Plin,double *Plin,int Nlin, double *k_Olin, double *Olin, int NOlin, double Sigma_smooth,double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *spacing_mask, char *path_to_mask1, double k0[], double k2[], double k4[], int Npoly, fftw_plan plan1, fftw_plan plan2, double kmin_theo, double kmax_theo, double k0min , double k0max,double k2min, double k2max,double k4min, double k4max, int wiggle, char *spacing_data,char *spacing_theory)
{

int i,j;
double ptheo,ptheo2,ptheo4,k,junk;
double XMIN1[1]={-1};
double XMAX1[1]={+1};
double precision=1e-1;
int Nelements;
int NeffP;
int NeffP_max;
double kmax,kmin;
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

//if(kmax==k0max && kmin==k0min){NeffP=(NeffP0+25-2+(int)(k0min/(k0max-k0min)*(NeffP0-1.)))*factor_for_sampling;NeffP_max=NeffP0;}
//if(kmax==k2max && kmin==k2min){NeffP=(NeffP2+25-2+(int)(k2min/(k2max-k2min)*(NeffP2-1.)))*factor_for_sampling;NeffP_max=NeffP2;}
//if(kmax==k4max && kmin==k4min){NeffP=(NeffP4+25-2+(int)(k4min/(k4max-k4min)*(NeffP4-1.)))*factor_for_sampling;NeffP_max=NeffP4;}
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

NeffP=Nlin;

}

if(NeffP==0){printf("Warning, crossing intervals among P0 (%lf<k<%lf), P2 (%lf<k<%lf), P4 (%lf<k<%lf) make impossible to determine NeffP within mask application. Fix this\n",k0min,k0max,k2min,k2max,k4min,k4max);exit(0);}

//printf("NeffP=%d\n",NeffP);

NeffP0=NeffP0*modeP0;
NeffP2=NeffP2*modeP2;
NeffP4=NeffP4*modeP4;

        f_params *function_parameters;

        function_parameters = (f_params *) malloc(sizeof(f_params));

(*function_parameters).Nlin=Nlin;
(*function_parameters).NOlin=NOlin;
(*function_parameters).parameters1=parameters1;
(*function_parameters).k_Plin=k_Plin;
(*function_parameters).Plin=Plin;
(*function_parameters).k_Olin=k_Olin;
(*function_parameters).Olin=Olin;
(*function_parameters).Sigma_smooth=Sigma_smooth;
(*function_parameters).wiggle=wiggle;
(*function_parameters).Npoly=Npoly;

(*function_parameters).modeP0=modeP0;
(*function_parameters).modeP2=modeP2;
(*function_parameters).modeP4=modeP4;
(*function_parameters).spacing=spacing_theory;

//printf("%d %d %d\n",modeP0,modeP2,modeP4);
//printf("%lf %lf, %lf %lf, %lf %lf\n",kmin_data0,kmax_data0,kmin_data2,kmax_data2,kmin_data4,kmax_data4);

//printf("%d %s\n",Nlin,path_to_mask1);

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP0;}

if(modeP0==1)
{

    for(j=0;j<Nelements;j++)
    {
//       if(strcmp(path_to_mask1, "none") != 0){(*function_parameters).kinput=k_Plin[j];k_theo[j]=k_Plin[j];k=k_Plin[j];}
    if(strcmp(path_to_mask1, "none") != 0){

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"irregular") == 0){k_theo[j]=k_Plin[j];k=k_theo[j];(*function_parameters).kinput=k_theo[j];}

}
       else{k_theo0[j]=k0[j];(*function_parameters).kinput=k_theo0[j];k=k0[j];}

         (*function_parameters).mode=0;
         adapt_integrate(1,integralPaniso,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);
//ptheo=10000;

for(i=1;i<Npoly+1;i++){ptheo=ptheo+parameters1[5+i]*pow(k,2-i);}
                                          

         P_theo0[j]=ptheo;


    }

}

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP2;}

if(modeP2==1)
{

    for(j=0;j<Nelements;j++)
    {
//       if(strcmp(path_to_mask1, "none") != 0){(*function_parameters).kinput=k_Plin[j];k_theo[j]=k_Plin[j];k=k_Plin[j];}
    if(strcmp(path_to_mask1, "none") != 0){


//k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"irregular") == 0){k_theo[j]=k_Plin[j];k=k_theo[j];(*function_parameters).kinput=k_theo[j];}


}
       else{k_theo2[j]=k2[j];(*function_parameters).kinput=k_theo2[j];k=k2[j];}

         (*function_parameters).mode=2;
         adapt_integrate(1,integralPaniso,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);

//ptheo=10000;

for(i=1;i<Npoly+1;i++){ptheo=ptheo+parameters1[5+i+(Npoly)*modeP0]*pow(k,2-i);}

         P_theo2[j]=ptheo;

    }

}

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP4;}

if(modeP4==1)
{
    for(j=0;j<Nelements;j++)
    {
//       if(strcmp(path_to_mask1, "none") != 0){(*function_parameters).kinput=k_Plin[j];k_theo[j]=k_Plin[j];k=k_Plin[j];}
  if(strcmp(path_to_mask1, "none") != 0){

//k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"irregular") == 0){k_theo[j]=k_Plin[j];k=k_theo[j];(*function_parameters).kinput=k_theo[j];}

}
       else{k_theo4[j]=k4[j];(*function_parameters).kinput=k_theo4[j];k=k4[j];}

         (*function_parameters).mode=4;
         adapt_integrate(1,integralPaniso,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);

for(i=1;i<Npoly+1;i++){ptheo=ptheo+parameters1[5+i+(Npoly)*modeP0+(Npoly)*modeP2]*pow(k,2-i);}               

         P_theo4[j]=ptheo;

    }

}




free(function_parameters);
//exit(0);
//for(i=0;i<Nelements;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i],P_theo2[i]);}

if(strcmp(path_to_mask1, "none") != 0)
{
//apply_mask(type_of_analysis,modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,Nlin, pos, W0, W2, W4, W6, W8, Nmask, plan1, plan2, kmin, kmax, kmin_data0, kmax_data0, kmin_data2, kmax_data2, kmin_data4, kmax_data4);
apply_mask(type_of_analysis,modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,NeffP, pos, W0, W2, W4, W6, W8, Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,Sigma_smooth);

}
//for(i=0;i<Nelements;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i],P_theo2[i]);}

//printf("%lf %lf %lf\n",k_theo[300],P_theo0[300],P_theo2[300]);
//exit(0);
}
*/

void do_Ptheo_multiple_aniso(char *type_BAO_fit, char *type_of_analysis, char *fit_BAO,int modeP0,int modeP2,int modeP4, double k_theo[],double k_theo0[], double k_theo2[],double k_theo4[], double P_theo0[], double P_theo2[], double P_theo4[],int NeffP0,int NeffP2,int NeffP4,int factor_for_sampling, double *parameters1,double *k_Plin,double *Plin,int Nlin, double *k_Olin, double *Olin, int NOlin, double Sigma_smooth,double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *spacing_mask, char *path_to_mask1, double k0[], double k2[], double k4[], int Npoly, fftw_plan plan1, fftw_plan plan2, double kmin_theo, double kmax_theo, double k0min , double k0max,double k2min, double k2max,double k4min, double k4max, int wiggle, char *spacing_data,char *spacing_theory)
{

int i,j;
double ptheo,ptheo2,ptheo4,k,junk;
double XMIN1[1]={-1};
double XMAX1[1]={+1};
double precision=1e-1;
int Nelements;
int NeffP;
int NeffP_max;
double kmax,kmin;
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

NeffP=Nlin;

}
if(NeffP==0){printf("Warning, crossing intervals among P0 (%lf<k<%lf), P2 (%lf<k<%lf), P4 (%lf<k<%lf) make impossible to determine NeffP within mask application. Fix this\n",k0min,k0max,k2min,k2max,k4min,k4max);exit(0);}
NeffP0=NeffP0*modeP0;
NeffP2=NeffP2*modeP2;
NeffP4=NeffP4*modeP4;

        f_params *function_parameters;

        function_parameters = (f_params *) malloc(sizeof(f_params));

(*function_parameters).Nlin=Nlin;
(*function_parameters).NOlin=NOlin;
(*function_parameters).parameters1=parameters1;
(*function_parameters).k_Plin=k_Plin;
(*function_parameters).Plin=Plin;
(*function_parameters).k_Olin=k_Olin;
(*function_parameters).Olin=Olin;
(*function_parameters).Sigma_smooth=Sigma_smooth;
(*function_parameters).wiggle=wiggle;
(*function_parameters).Npoly=Npoly;
(*function_parameters).modeP0=modeP0;
(*function_parameters).modeP2=modeP2;
(*function_parameters).modeP4=modeP4;
(*function_parameters).spacing=spacing_theory;
if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP0;}

if(modeP0==1)
{

    for(j=0;j<Nelements;j++)
    {
    if(strcmp(path_to_mask1, "none") != 0){

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"irregular") == 0){k_theo[j]=k_Plin[j];k=k_theo[j];(*function_parameters).kinput=k_theo[j];}

}
       else{k_theo0[j]=k0[j];(*function_parameters).kinput=k_theo0[j];k=k0[j];}

         (*function_parameters).mode=0;
         adapt_integrate(1,integralPaniso,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);
for(i=1;i<Npoly+1;i++){ptheo=ptheo+parameters1[5+i]*pow(k,2-i);}


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

if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"irregular") == 0){k_theo[j]=k_Plin[j];k=k_theo[j];(*function_parameters).kinput=k_theo[j];}


}
       else{k_theo2[j]=k2[j];(*function_parameters).kinput=k_theo2[j];k=k2[j];}

         (*function_parameters).mode=2;
         adapt_integrate(1,integralPaniso,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);
for(i=1;i<Npoly+1;i++){ptheo=ptheo+parameters1[5+i+(Npoly)*modeP0]*pow(k,2-i);}

         P_theo2[j]=ptheo;

    }

}

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP4;}

if(modeP4==1)
{
    for(j=0;j<Nelements;j++)
    {
  if(strcmp(path_to_mask1, "none") != 0){
if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling) );k=k_theo[j];(*function_parameters).kinput=k_theo[j];}
if(strcmp(spacing_data,"irregular") == 0){k_theo[j]=k_Plin[j];k=k_theo[j];(*function_parameters).kinput=k_theo[j];}

}
       else{k_theo4[j]=k4[j];(*function_parameters).kinput=k_theo4[j];k=k4[j];}

         (*function_parameters).mode=4;
         adapt_integrate(1,integralPaniso,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);

for(i=1;i<Npoly+1;i++){ptheo=ptheo+parameters1[5+i+(Npoly)*modeP0+(Npoly)*modeP2]*pow(k,2-i);}

         P_theo4[j]=ptheo;

    }

}




free(function_parameters);
if(strcmp(path_to_mask1, "none") != 0)
{
apply_mask(type_of_analysis,modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,NeffP, pos, W0, W2, W4, W6, W8, Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,Sigma_smooth,1);
}
}



void do_Ptheo_multiple_iso(char *type_BAO_fit, char *type_of_analysis, char *fit_BAO,int modeP0,int modeP2,int modeP4, double k_theo[],double k_theo0[], double k_theo2[], double k_theo4[], double P_theo0[], double P_theo2[], double P_theo4[],int NeffP0, int NeffP2, int NeffP4,int factor_for_sampling, double *parameters1,double *k_Plin,double *Plin,int Nlin, double *k_Olin, double *Olin, int NOlin,double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *spacing_mask,char *path_to_mask1,double k0[], double k2[], double k4[], int Npoly, fftw_plan plan1, fftw_plan plan2, double kmin_theo , double kmax_theo, double k0min , double k0max,double k2min, double k2max,double k4min, double k4max, int wiggle, char *spacing_data, char *spacing_theory, double Sigma_smooth)
{
//if wiggle=0 no wiggle in the function. if wiggle=1 wiggle
int i,j;
double func0,func2,func4,olin_eff0,olin_eff2,olin_eff4;
double sigma_nl0,sigma_nl2,sigma_nl4;
double kprime0,kprime2,kprime4,arg0,arg2,arg4;
double alpha0,alpha2,alpha4;
double plin,kused;
int offset, Nelements, Nmultipoles;
double kmax,kmin;
int NeffP;
int NeffP_max;
int Nlin1;
double w0lin,w1lin,w2lin;
int interpolation_order,shiftN;

interpolation_order=1;

    if(interpolation_order==1){shiftN=1;}
    if(interpolation_order==2){shiftN=2;}


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
//if(kmax==k0max && kmin==k0min){NeffP=(NeffP0+25-2+(int)(k0min/(k0max-k0min)*(NeffP0-1.)))*factor_for_sampling;NeffP_max=NeffP0;}
//if(kmax==k2max && kmin==k2min){NeffP=(NeffP2+25-2+(int)(k2min/(k2max-k2min)*(NeffP2-1.)))*factor_for_sampling;NeffP_max=NeffP2;}
//if(kmax==k4max && kmin==k4min){NeffP=(NeffP4+25-2+(int)(k4min/(k4max-k4min)*(NeffP4-1.)))*factor_for_sampling;NeffP_max=NeffP4;}

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
NeffP=Nlin;
}



if(NeffP==0){printf("Warning, crossing intervals among P0 (%lf<k<%lf), P2 (%lf<k<%lf), P4 (%lf<k<%lf) make impossible to determine NeffP within mask application. Fix this\n",k0min,k0max,k2min,k2max,k4min,k4max);exit(0);}

//printf("NeffP=%d\n",NeffP);

NeffP0=NeffP0*modeP0;
NeffP2=NeffP2*modeP2;
NeffP4=NeffP4*modeP4;

Nmultipoles=0;
if(modeP0==1){Nmultipoles++;}
if(modeP2==1){Nmultipoles++;}
if(modeP4==1){Nmultipoles++;}

offset=0;
if(Nmultipoles==1){alpha0=parameters1[0];offset++;sigma_nl0=parameters1[1];offset++;}//note that alpha0 can be alpha0, alpha2 or alpha4 when Nmultipoles is eq. to 1 (idem for sigmanl0)
else
{
alpha0=pow(parameters1[0],1./3.)*pow(parameters1[1],2./3.);
alpha2=pow(parameters1[0],3./5.)*pow(parameters1[1],2./5.);
alpha4=pow(parameters1[0],5./7.)*pow(parameters1[1],2./7.);
offset=2;
if(modeP0==1){
              sigma_nl0=parameters1[2];offset++;
              if(modeP2==1){
                            sigma_nl2=parameters1[3];offset++;
                            if(modeP4==1){sigma_nl4=parameters1[4];offset++;}
              }else{sigma_nl4=parameters1[3];offset++;}

}
else{
          sigma_nl2=parameters1[2];offset++;
          sigma_nl4=parameters1[3];offset++;
}

}

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{

if(NeffP0>=NeffP2 && NeffP0>=NeffP4){Nelements=NeffP0;}
if(NeffP2>=NeffP0 && NeffP2>=NeffP4){Nelements=NeffP2;}
if(NeffP4>=NeffP0 && NeffP4>=NeffP2){Nelements=NeffP4;}

}

    for(j=0;j<Nelements;j++)
    {
//if(strcmp(path_to_mask1, "none") != 0){k_theo[j]=k_Plin[j];plin=Plin[j];kused=k_theo[j];}
if(strcmp(path_to_mask1, "none") != 0){


if( strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);kused=k_theo[j];}
if( strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling));kused=k_theo[j];}
if( strcmp(spacing_data,"log10") == 0){k_theo[j]=pow( 10, (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling));kused=k_theo[j];}
if( strcmp(spacing_data,"irregular") == 0){k_theo[j]=k_Plin[j];kused=k_theo[j];}



    Nlin1=determine_N_singlearray(k_Plin,kused,Nlin,spacing_theory);

if(interpolation_order==1){
w1lin=determine_w1_singlearray(k_Plin,kused,Nlin1,spacing_theory);
}
if(interpolation_order==2){
w0lin=determine_w0_2ndorder_singlearray(k_Plin,kused,Nlin1,spacing_theory);
w1lin=determine_w1_2ndorder_singlearray(k_Plin,kused,Nlin1,spacing_theory);
w2lin=determine_w2_2ndorder_singlearray(k_Plin,kused,Nlin1,spacing_theory);
}

plin=P_interpol_fast(kused,Plin,Nlin,spacing_theory,interpolation_order,Nlin1,w0lin,w1lin,w2lin);
//printf("%lf %d %lf %lf %s\n",kused,Nlin1,w1lin,plin,spacing_theory);//exit(0);
}
else{  

if(NeffP0>=NeffP2 && NeffP0>=NeffP4){kused=k0[j];/*plin=P_interpol(kused,k_Plin,Plin,Nlin);*/}
if(NeffP2>=NeffP0 && NeffP2>=NeffP4){kused=k2[j];/*plin=P_interpol(kused,k_Plin,Plin,Nlin);*/}
if(NeffP4>=NeffP0 && NeffP4>=NeffP2){kused=k4[j];/*plin=P_interpol(kused,k_Plin,Plin,Nlin);*/}


    Nlin1=determine_N_singlearray(k_Plin,kused,Nlin,spacing_theory);

if(interpolation_order==1){
w1lin=determine_w1_singlearray(k_Plin,kused,Nlin1,spacing_theory);
}
if(interpolation_order==2){
w0lin=determine_w0_2ndorder_singlearray(k_Plin,kused,Nlin1,spacing_theory);
w1lin=determine_w1_2ndorder_singlearray(k_Plin,kused,Nlin1,spacing_theory);
w2lin=determine_w2_2ndorder_singlearray(k_Plin,kused,Nlin1,spacing_theory);
}

plin=P_interpol_fast(kused,Plin,Nlin,spacing_theory,interpolation_order,Nlin1,w0lin,w1lin,w2lin);
//printf("%lf %d %lf %lf %s\n",kused,Nlin1,w1lin,plin,spacing_theory);//exit(0);


}


       func0=0;
       func2=0;
       func4=0;
       for(i=0;i<Npoly+1;i++)
       {
          if(i==0)
          {
            
              if(Nmultipoles==1){func0=plin*parameters1[2];}
              else{func0=plin*parameters1[offset];}
              
              if(Nmultipoles>1){
              func2=plin*parameters1[offset+Npoly+1];
              }

              if(Nmultipoles>2){
              func4=plin*parameters1[offset+(Npoly+1)*2];
              }

          }
          else
          {
             if(Nmultipoles==1){func0=func0+parameters1[2+i]*pow(kused,2-i);}
             else{func0=func0+parameters1[offset+i]*pow(kused,2-i);}

              if(Nmultipoles>1){
              func2=func2+parameters1[offset+i+Npoly+1]*pow(kused,2-i);
              }
              if(Nmultipoles>2){
              func4=func4+parameters1[offset+i+(Npoly+1)*2]*pow(kused,2-i);
              }

          }

       }


          kprime0=kused/alpha0;

          if(Nmultipoles>1){
if(strcmp(fit_BAO, "P02") == 0){kprime2=kused/alpha2;}
if(strcmp(fit_BAO, "P04") == 0){kprime2=kused/alpha4;}
if(strcmp(fit_BAO, "P24") == 0){kprime0=kused/alpha2;kprime2=kused/alpha4;}

          }
          if(Nmultipoles>2){
          kprime2=kused/alpha2;
          kprime4=kused/alpha4;
          }


  Nlin1=determine_N_singlearray(k_Olin,kprime0,NOlin,spacing_theory);
if(interpolation_order==1){
w1lin=determine_w1_singlearray(k_Olin,kprime0,Nlin1,spacing_theory);
}
if(interpolation_order==2){
w0lin=determine_w0_2ndorder_singlearray(k_Olin,kprime0,Nlin1,spacing_theory);
w1lin=determine_w1_2ndorder_singlearray(k_Olin,kprime0,Nlin1,spacing_theory);
w2lin=determine_w2_2ndorder_singlearray(k_Olin,kprime0,Nlin1,spacing_theory);
}
olin_eff0=P_interpol_fast(kprime0,Olin,NOlin,spacing_theory,interpolation_order,Nlin1,w0lin,w1lin,w2lin);


          if(Nmultipoles>1){          

  Nlin1=determine_N_singlearray(k_Olin,kprime2,NOlin,spacing_theory);
if(interpolation_order==1){
w1lin=determine_w1_singlearray(k_Olin,kprime2,Nlin1,spacing_theory);
}
if(interpolation_order==2){
w0lin=determine_w0_2ndorder_singlearray(k_Olin,kprime2,Nlin1,spacing_theory);
w1lin=determine_w1_2ndorder_singlearray(k_Olin,kprime2,Nlin1,spacing_theory);
w2lin=determine_w2_2ndorder_singlearray(k_Olin,kprime2,Nlin1,spacing_theory);
}
olin_eff2=P_interpol_fast(kprime2,Olin,NOlin,spacing_theory,interpolation_order,Nlin1,w0lin,w1lin,w2lin);



          }
          if(Nmultipoles>2){
  
  Nlin1=determine_N_singlearray(k_Olin,kprime4,NOlin,spacing_theory);
if(interpolation_order==1){
w1lin=determine_w1_singlearray(k_Olin,kprime4,Nlin1,spacing_theory);
}
if(interpolation_order==2){
w0lin=determine_w0_2ndorder_singlearray(k_Olin,kprime4,Nlin1,spacing_theory);
w1lin=determine_w1_2ndorder_singlearray(k_Olin,kprime4,Nlin1,spacing_theory);
w2lin=determine_w2_2ndorder_singlearray(k_Olin,kprime4,Nlin1,spacing_theory);
}
olin_eff4=P_interpol_fast(kprime4,Olin,NOlin,spacing_theory,interpolation_order,Nlin1,w0lin,w1lin,w2lin);



        }
if(wiggle==0){olin_eff0=1;olin_eff2=1;olin_eff4=1;}

          arg0=-0.5*sigma_nl0*sigma_nl0*kused*kused;
          arg2=-0.5*sigma_nl2*sigma_nl2*kused*kused;
          arg4=-0.5*sigma_nl4*sigma_nl4*kused*kused;

       
if(strcmp(fit_BAO, "P0") == 0){
          P_theo0[j]=func0*(1.+(olin_eff0-1)*exp(arg0));if(strcmp(path_to_mask1, "none") == 0){k_theo0[j]=k0[j];}
}
if(strcmp(fit_BAO, "P2") == 0){
          P_theo2[j]=func0*(1.+(olin_eff0-1)*exp(arg0));if(strcmp(path_to_mask1, "none") == 0){k_theo2[j]=k2[j];}
}
if(strcmp(fit_BAO, "P4") == 0){
          P_theo4[j]=func0*(1.+(olin_eff0-1)*exp(arg0));if(strcmp(path_to_mask1, "none") == 0){k_theo4[j]=k4[j];}
}
if(strcmp(fit_BAO, "P02") == 0){

if(strcmp(path_to_mask1, "none") != 0){
          P_theo0[j]=func0*(1.+(olin_eff0-1)*exp(arg0));
          P_theo2[j]=func2*(1.+(olin_eff2-1)*exp(arg2));
}
else{
if(j<NeffP0){P_theo0[j]=func0*(1.+(olin_eff0-1)*exp(arg0));k_theo0[j]=k0[j];}
if(j<NeffP2){P_theo2[j]=func2*(1.+(olin_eff2-1)*exp(arg2));k_theo2[j]=k2[j];}
}



}
if(strcmp(fit_BAO, "P04") == 0){

if(strcmp(path_to_mask1, "none") != 0){
          P_theo0[j]=func0*(1.+(olin_eff0-1)*exp(arg0));
          P_theo4[j]=func2*(1.+(olin_eff2-1)*exp(arg4));
}
else{
if(j<NeffP0){P_theo0[j]=func0*(1.+(olin_eff0-1)*exp(arg0));k_theo0[j]=k0[j];}
if(j<NeffP4){P_theo4[j]=func2*(1.+(olin_eff2-1)*exp(arg4));k_theo4[j]=k4[j];}
}
}
if(strcmp(fit_BAO, "P24") == 0){

if(strcmp(path_to_mask1, "none") != 0){
          P_theo2[j]=func0*(1.+(olin_eff0-1)*exp(arg2));
          P_theo4[j]=func2*(1.+(olin_eff2-1)*exp(arg4));
}
else{
if(j<NeffP2){P_theo2[j]=func0*(1.+(olin_eff0-1)*exp(arg2));k_theo2[j]=k2[j];}
if(j<NeffP4){P_theo4[j]=func2*(1.+(olin_eff2-1)*exp(arg4));k_theo4[j]=k4[j];}
}

}
if(strcmp(fit_BAO, "P024") == 0){
if(strcmp(path_to_mask1, "none") != 0){
          P_theo0[j]=func0*(1.+(olin_eff0-1)*exp(arg0));
          P_theo2[j]=func2*(1.+(olin_eff2-1)*exp(arg2));
          P_theo4[j]=func4*(1.+(olin_eff4-1)*exp(arg4));
}
else{
if(j<NeffP0){P_theo0[j]=func0*(1.+(olin_eff0-1)*exp(arg0));k_theo0[j]=k0[j];}
if(j<NeffP2){P_theo2[j]=func2*(1.+(olin_eff2-1)*exp(arg2));k_theo2[j]=k2[j];}
if(j<NeffP4){P_theo4[j]=func4*(1.+(olin_eff4-1)*exp(arg4));k_theo4[j]=k4[j];}

}
}

    
     }

//printf("%d\n",Nelements);
//for(i=0;i<Nelements;i++){printf("%lf %lf %lf %lf\n",k_theo0[i],k_theo2[i],P_theo0[i],P_theo2[i]);}
//for(i=0;i<Nelements;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i],P_theo2[i]);}
//exit(0);

if(strcmp(path_to_mask1, "none") != 0)
{
apply_mask(type_of_analysis,modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,NeffP, pos, W0, W2, W4, W6, W8, Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,Sigma_smooth,1);
}

//for(i=0;i<Nelements;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i],P_theo2[i]);}
//exit(0);

}

int get_convergence(char *filename, int Nparams,long int lines)
{
//dividir la cadena en 6 chunks de al menos 1e5 y comprobar la convergancia entre estos chunks independientes
int convergence,chunks,j,l;
long int i,min_lines,eff_lines;
FILE *f;
int **weight;
double ***param;
double **mean_chunk;
double *R_convergence,*W_convergence,*B_convergence,*mean;
long int *Neff;
long int Nefftot;
double Rlimit,Rtot;
double parameters_convergence[3];

set_parameters_convergence(parameters_convergence);
Rlimit=parameters_convergence[0];
min_lines=(long int)(parameters_convergence[1]);
chunks=(int)(parameters_convergence[2]);

convergence=0;
if(lines<min_lines*chunks){printf("Not enough steps for convergence test yet: lines=%ld, min_lines=%ld, chunks=%d, thread_id=%d\n",lines,min_lines,chunks,omp_get_thread_num());return convergence;}

eff_lines=(long int)(lines*1./chunks*1.);

R_convergence =  (double *) calloc(Nparams, sizeof(double));
W_convergence =  (double *) calloc(Nparams, sizeof(double));
B_convergence =  (double *) calloc(Nparams, sizeof(double));
mean =  (double *) calloc(Nparams, sizeof(double));//absolute mean
Neff =  (long int *) calloc(chunks, sizeof(long int));

param=(double ***) calloc(chunks, sizeof(double**));
mean_chunk= (double **) calloc(chunks, sizeof(double*));
weight= (int **) calloc(chunks, sizeof(int *));

for(j=0;j<chunks;j++)
{
weight[j]= (int *) calloc(eff_lines, sizeof(int));
mean_chunk[j] = (double *) calloc(Nparams, sizeof(double));
param[j] = (double **) calloc(eff_lines, sizeof(double*));
 
  for(i=0;i<eff_lines;i++)
  {
      param[j][i] = (double *) calloc(Nparams, sizeof(double));
  }

}

f=fopen(filename,"r");
Nefftot=0;
for(j=0;j<chunks;j++)
{  
  for(i=0;i<eff_lines;i++)
  {
      fscanf(f,"%d %*f ",&weight[j][i]);
      Neff[j]=Neff[j]+weight[j][i];
      Nefftot=Nefftot+weight[j][i];
      for(l=0;l<Nparams;l++)
      {
         fscanf(f,"%lf ",&param[j][i][l]);mean[l]=mean[l]+param[j][i][l]*weight[j][i];mean_chunk[j][l]=mean_chunk[j][l]+param[j][i][l]*weight[j][i];
         if(l==Nparams-1){fscanf(f,"\n");}

      }    
  }

}
fclose(f);

for(l=0;l<Nparams;l++)
{  
   for(j=0;j<chunks;j++)
   {
     mean_chunk[j][l]=mean_chunk[j][l]/Neff[j]*1.;    
   }
   mean[l]=mean[l]/Nefftot*1.;
   B_convergence[l]=0;
   for(j=0;j<chunks;j++)
   {
     B_convergence[l]=B_convergence[l]+pow(mean_chunk[j][l]-mean[l],2)/(chunks-1.);
   }
   
}

for(l=0;l<Nparams;l++)
{
   W_convergence[l]=0;
   for(j=0;j<chunks;j++)
   {

      for(i=0;i<eff_lines;i++)
      {
         W_convergence[l]=W_convergence[l]+pow(param[j][i][l]-mean_chunk[j][l],2)*weight[j][i];
      }
  }
 W_convergence[l]=W_convergence[l]/((Nefftot-chunks)*1.);
}

printf("\n===Convergence Status===\n\n");
Rtot=0;
for(l=0;l<Nparams;l++)
{
R_convergence[l]=( ( (Nefftot*1./chunks*1.)-1.)/(Nefftot*1./chunks*1.)*1.*W_convergence[l]+B_convergence[l]*(1. + 1./chunks*1. )   )/(W_convergence[l]);
printf("R[%d]= %lf\n",l,R_convergence[l]);
Rtot=Rtot+R_convergence[l]/Nparams*1.;
}

printf("\n<R>= %lf\n",Rtot);


//Condition on the individual convergence of parameters
for(l=0;l<Nparams;l++)
{
if(R_convergence[l]<Rlimit){convergence++;}
}

if(convergence==Nparams){printf("Convergence reached!\n");convergence=1;}
else{
if(Nparams-convergence>1){printf("Convergence not reached, %d parameters missing to converge.\n",Nparams-convergence);}
else{printf("Convergence not reached, %d parameter missing to converge.\n",Nparams-convergence);}
convergence=0;}


free(R_convergence);
free(W_convergence);
free(B_convergence);
free(mean);
free(Neff);

freeTokens2(param, chunks, eff_lines);
freeTokens(mean_chunk,chunks);
freeTokensInt(weight,chunks);

return convergence;
}

int get_convergence_parallel(int nthreads, char *identifier, int Nparams)
{
//obtener nthreads chunks de las nthreads cadenas de al menos 1e5 y comprobar la convergancia entre estos chunks independientes
int convergence,chunks,j,l,ithread;
long int i,min_lines,eff_lines;
long int lines[nthreads];
FILE *f;
char filename[2000];
int **weight;
double ***param;
double **mean_chunk;
double *R_convergence,*W_convergence,*B_convergence,*mean;
long int *Neff;
long int Nefftot;
double Rlimit,Rtot;
double parameters_convergence[3];

set_parameters_convergence_parallel(nthreads,parameters_convergence);
Rlimit=parameters_convergence[0];
min_lines=(long int)(parameters_convergence[1]);
chunks=(int)(parameters_convergence[2]);

//get the number of lines of each chain file
for (ithread=0;ithread<nthreads;ithread++)
{
  sprintf(filename,"%s__%d.txt",identifier,ithread+1);
  lines[ithread]=countlinesLI(filename);
} 

convergence=0;
if(getMin(lines, nthreads)<min_lines){printf("Not enough steps for convergence test yet: lines=%ld, min_lines=%ld, chunks=%d, thread_id=%d\n",getMin(lines, nthreads),min_lines,chunks,omp_get_thread_num());return convergence;}

eff_lines=getMin(lines, nthreads);

R_convergence =  (double *) calloc(Nparams, sizeof(double));
W_convergence =  (double *) calloc(Nparams, sizeof(double));
B_convergence =  (double *) calloc(Nparams, sizeof(double));
mean =  (double *) calloc(Nparams, sizeof(double));//absolute mean
Neff =  (long int *) calloc(chunks, sizeof(long int));

param=(double ***) calloc(chunks, sizeof(double**));
mean_chunk= (double **) calloc(chunks, sizeof(double*));
weight= (int **) calloc(chunks, sizeof(int *));

for(j=0;j<chunks;j++)
{
weight[j]= (int *) calloc(eff_lines, sizeof(int));
mean_chunk[j] = (double *) calloc(Nparams, sizeof(double));
param[j] = (double **) calloc(eff_lines, sizeof(double*));
 
  for(i=0;i<eff_lines;i++)
  {
      param[j][i] = (double *) calloc(Nparams, sizeof(double));
  }

}


Nefftot=0;
for(j=0;j<chunks;j++)
{  
  sprintf(filename,"%s__%d.txt",identifier,j+1);
  f=fopen(filename,"r");
  skiplines(f,lines[j]-eff_lines);
  for(i=0;i<eff_lines;i++)
  {
      fscanf(f,"%d %*f ",&weight[j][i]);
      Neff[j]=Neff[j]+weight[j][i];
      Nefftot=Nefftot+weight[j][i];
      for(l=0;l<Nparams;l++)
      {
         fscanf(f,"%lf ",&param[j][i][l]);mean[l]=mean[l]+param[j][i][l]*weight[j][i];mean_chunk[j][l]=mean_chunk[j][l]+param[j][i][l]*weight[j][i];
         if(l==Nparams-1){fscanf(f,"\n");}

      }    
  }
  fclose(f);
}


for(l=0;l<Nparams;l++)
{  
   for(j=0;j<chunks;j++)
   {
     mean_chunk[j][l]=mean_chunk[j][l]/Neff[j]*1.;    
   }
   mean[l]=mean[l]/Nefftot*1.;
   B_convergence[l]=0;
   for(j=0;j<chunks;j++)
   {
     B_convergence[l]=B_convergence[l]+pow(mean_chunk[j][l]-mean[l],2)/(chunks-1.);
   }
   
}

for(l=0;l<Nparams;l++)
{
   W_convergence[l]=0;
   for(j=0;j<chunks;j++)
   {

      for(i=0;i<eff_lines;i++)
      {
         W_convergence[l]=W_convergence[l]+pow(param[j][i][l]-mean_chunk[j][l],2)*weight[j][i];
      }
  }
 W_convergence[l]=W_convergence[l]/((Nefftot-chunks)*1.);
}

printf("\n===Convergence Status===\n\n");
Rtot=0;
for(l=0;l<Nparams;l++)
{
R_convergence[l]=( ( (Nefftot*1./chunks*1.)-1.)/(Nefftot*1./chunks*1.)*1.*W_convergence[l]+B_convergence[l]*(1. + 1./chunks*1. )   )/(W_convergence[l]);
printf("R[%d]= %lf\n",l,R_convergence[l]);
Rtot=Rtot+R_convergence[l]/Nparams*1.;
}

printf("\n<R>= %lf\n",Rtot);


//Condition on the individual convergence of parameters
for(l=0;l<Nparams;l++)
{
if(R_convergence[l]<Rlimit){convergence++;}
}

if(convergence==Nparams){printf("Convergence reached!\n");convergence=1;}
else{
if(Nparams-convergence>1){printf("Convergence not reached, %d parameters missing to converge.\n",Nparams-convergence);}
else{printf("Convergence not reached, %d parameter missing to converge.\n",Nparams-convergence);}
convergence=0;}


free(R_convergence);
free(W_convergence);
free(B_convergence);
free(mean);
free(Neff);

freeTokens2(param, chunks, eff_lines);
freeTokens(mean_chunk,chunks);
freeTokensInt(weight,chunks);

return convergence;
}



void read_prop_cov(double **vector_buffer,long int N_max, int trial_mcmc, char *path_to_cov,double Cov_prop[],double vector_cov_prop[],int N_Cov_prop)
{
// fill in Cov_prop AND vector_cov_prop
long int Nsteps;
long int weighted_steps;
int i,j;
long int l;
double **elements;
int *weight;
char *ending;
FILE *f;



ending = strrchr(path_to_cov, '.');
if ((trial_mcmc==0) && (strcmp(ending, ".cov") == 0)){
  f=fopen(path_to_cov,"r");
  fscanf(f,"%*s %*s \n");
  for(j=0;j<N_Cov_prop;j++)
  {
     fscanf(f,"%lf\n", &vector_cov_prop[j]);
  }
  fscanf(f,"%*s %*s \n");
  for(j=0;j<N_Cov_prop;j++)
  {
     for(i=0;i<N_Cov_prop;i++)
     {
       fscanf(f,"%lf\t", &Cov_prop[j+i*N_Cov_prop]);
     }
     fscanf(f,"\n");
  }

  fclose(f);
  for(j=0;j<N_Cov_prop;j++)
  {
   printf("Element %d x=%lf pm %lf\n",j,vector_cov_prop[j],sqrt(Cov_prop[j+j*N_Cov_prop]));
  }
  
}
else {
  if(trial_mcmc==0)
  {
    Nsteps=countlinesLI(path_to_cov);
    elements = (double**) calloc(N_Cov_prop,sizeof(double*));
    weight = (int*) calloc(Nsteps,sizeof(int));

    for(j=0;j<N_Cov_prop;j++)
    {
       elements[j] = (double*) calloc(Nsteps,sizeof(double));
    }

    f=fopen(path_to_cov,"r");

    weighted_steps=0;
    for(l=0;l<Nsteps;l++)
    {
        for(j=-1;j<N_Cov_prop;j++)
        {
           if(j==-1){fscanf(f,"%d %*f\t", &weight[l]);}
           if(j>=0 && j<N_Cov_prop-1){fscanf(f,"%lf\t", &elements[j][l]);}
           if(j==N_Cov_prop-1){fscanf(f,"%lf\n", &elements[j][l]);} 
           if(j>=0){vector_cov_prop[j]=vector_cov_prop[j]+elements[j][l]*weight[l];}
        }
        weighted_steps=weighted_steps+weight[l];

    }
    fclose(f);
  }
  else
  {
  Nsteps=N_max;
  weighted_steps=0;
  for(l=0;l<Nsteps;l++)
  {
      for(j=0;j<N_Cov_prop;j++)
      {
         vector_cov_prop[j]=vector_cov_prop[j]+vector_buffer[j+2][l]*vector_buffer[0][l];
      }
      weighted_steps=weighted_steps+(long int)(vector_buffer[0][l]);


  }

  }

  for(j=0;j<N_Cov_prop;j++)
  {
     for(i=0;i<N_Cov_prop;i++)
     {
       Cov_prop[j+i*N_Cov_prop]=0;
       for(l=0;l<Nsteps;l++)
       {
  if(trial_mcmc==0){Cov_prop[j+i*N_Cov_prop]=Cov_prop[j+i*N_Cov_prop]+weight[l]*(vector_cov_prop[i]/weighted_steps*1.-elements[i][l])*(vector_cov_prop[j]/weighted_steps*1.-elements[j][l])/(weighted_steps*1.-1.);}
  else{Cov_prop[j+i*N_Cov_prop]=Cov_prop[j+i*N_Cov_prop]+vector_buffer[0][l]*(vector_cov_prop[i]/weighted_steps*1.-vector_buffer[i+2][l])*(vector_cov_prop[j]/weighted_steps*1.-vector_buffer[j+2][l])/(weighted_steps*1.-1.);}

       }
     }
  }


     for(j=0;j<N_Cov_prop;j++)
     {
       vector_cov_prop[j]=vector_cov_prop[j]/weighted_steps*1.;
       printf("Element %d x=%lf pm %lf\n",j,vector_cov_prop[j],sqrt(Cov_prop[j+j*N_Cov_prop]));
     }

  if(trial_mcmc==0){
  freeTokens(elements,N_Cov_prop);
  free(weight);
  }
}
}

void write_prop_cov(char *path_to_chain, int N_cov)
{
// infer cov AND vector_cov from chain file and write them into new mcmc...1234.cov file
long int Nsteps;
long int weighted_steps;
int i,j;
long int l;
double *cov, *vector_cov;
double **elements;
int *weight;
char name_wo_extension[2000], path_to_cov[2000];
FILE *f;

// read the chain file and store in memory
Nsteps=countlinesLI(path_to_chain);
elements = (double**) calloc(N_cov,sizeof(double*));
weight = (int*) calloc(Nsteps,sizeof(int));
cov = (double*) calloc(N_cov*N_cov, sizeof(double));
vector_cov = (double*) calloc(N_cov, sizeof(double));//zero-inizialized

for(j=0;j<N_cov;j++)
{
   elements[j] = (double*) calloc(Nsteps,sizeof(double));
}

f=fopen(path_to_chain,"r");

weighted_steps=0;
for(l=0;l<Nsteps;l++)
{
    for(j=-1;j<N_cov;j++)
    {
       if(j==-1){fscanf(f,"%d %*f\t", &weight[l]);}
       if(j>=0 && j<N_cov-1){fscanf(f,"%lf\t", &elements[j][l]);}
       if(j==N_cov-1){fscanf(f,"%lf\n", &elements[j][l]);} 
       if(j>=0){vector_cov[j]=vector_cov[j]+elements[j][l]*weight[l];}
    }
    weighted_steps=weighted_steps+weight[l];

}
fclose(f);

// from the stored steps, compute covariance

for(j=0;j<N_cov;j++)
{
   for(i=0;i<N_cov;i++)
   {
     cov[j+i*N_cov]=0;
     for(l=0;l<Nsteps;l++)
     {
        cov[j+i*N_cov]=cov[j+i*N_cov]+weight[l]*(vector_cov[i]/weighted_steps*1.-elements[i][l])*(vector_cov[j]/weighted_steps*1.-elements[j][l])/(weighted_steps*1.-1.);
     }
   }
}

freeTokens(elements,N_cov);
free(weight);


// write the mean values and the covariance into a file
strcpy(name_wo_extension, path_to_chain);
strtok(name_wo_extension, ".");
sprintf(path_to_cov,"%s.cov",name_wo_extension);
f=fopen(path_to_cov,"w");


fprintf(f,"#Mean parameters: \n");
for(j=0;j<N_cov;j++)
{
 vector_cov[j]=vector_cov[j]/weighted_steps*1.;
 fprintf(f,"%lf\n",vector_cov[j]);
}
fprintf(f,"#Covariance Matrix: \n");
for(j=0;j<N_cov;j++)
{
   for(i=0;i<N_cov;i++)
   {
     fprintf(f,"%lf\t",cov[j+i*N_cov]);
   }
   fprintf(f,"\n");
}

fclose(f);
}


void generate_rotation_matrix(int n,double *Cov_prop, double *vector_mean, gsl_matrix transform_inverse[], gsl_matrix transform[],int trial_mcmc)
{
int i,j;
double eval_i;
double tij,tinverseij;
int s;
double *error;


if(trial_mcmc==1)
{
error= (double*) calloc(n,sizeof(double));
//call struct for inizializing errors
set_proposal_error(error,n);

          for(i=0;i<n;i++)
          {
          for(j=0;j<n;j++)
          {
                 gsl_matrix_set (transform, i, j, 0);
                 if(i==j){gsl_matrix_set (transform, i, j, error[i] );}
                 gsl_matrix_set (transform_inverse, i, j, 0);
                 if(i==j){gsl_matrix_set (transform_inverse, i, j, 1./(error[i]));}

          }
          }
//Identity matrix with sqrt(variance) in diagonal. 
free(error);
}
else
{


  gsl_matrix_view m = gsl_matrix_view_array (Cov_prop, n, n);

  gsl_vector *eval = gsl_vector_alloc (n);//autovalores
  gsl_matrix *evec = gsl_matrix_alloc (n, n);//autovectores

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);

  gsl_matrix * L;

  gsl_eigen_symmv (&m.matrix, eval, evec, w);

  gsl_eigen_symmv_free (w);

  gsl_eigen_symmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_ASC);


           L = gsl_matrix_alloc (n, n);

          for(i=0;i<n;i++)
          {
          for(j=0;j<n;j++)
          {
                 gsl_matrix_set (L, i, j, 0);

          }
          }
  for (i = 0; i < n; i++)
      {
         eval_i = gsl_vector_get (eval, i);
           gsl_vector_view evec_i = gsl_matrix_column (evec, i);

        gsl_matrix_set (L, i, i, eval_i);//Unitaria amb autovalors a la diagonal

      }


          for(i=0;i<n;i++)
          {
          for(j=0;j<n;j++)
          {

                tij=sqrt(gsl_vector_get (eval, j))*gsl_matrix_get (evec, i, j);//Rotation matrix^-1= S*R^-1, where S is sqrt(L), where L is unity with eigenvalues in diagonal. 
                tinverseij=1./sqrt(gsl_vector_get (eval, i))*gsl_matrix_get (evec, j, i);//Rotation matrix^-1= S*R^-1, where S is sqrt(L), where L is unity with eigenvalues in diagonal. 

                gsl_matrix_set (transform, i, j, tij);
                gsl_matrix_set (transform_inverse, i, j, tinverseij);

          }
          }

  gsl_matrix_free (evec);
  gsl_vector_free (eval);

printf("Proposal Covariance Ready...\n");

          for(i=0;i<n;i++)
          {
          for(j=0;j<n;j++)
          {

if(j==n-1){printf("%e\n", gsl_matrix_get (transform_inverse, i, j) );}
else{printf("%e\t", gsl_matrix_get (transform_inverse, i, j) );}

         }

         }

  gsl_matrix_free(L);
}

}

//void mcmc_kernel(int nthreads, char *type_BAO_fit,char *type_of_analysis,int trial_mcmc,double **vector_buffer,double fraction, int N_Cov_prop, gsl_matrix *transform, gsl_matrix *transform_inverse, double *vector_mean, char *name_file_output_mcmc, char *fit_BAO,double *k_Plin,double *Plin,int N_Plin, double *k_Olin, double *Olin, int N_Olin, double *pos, double *W0, double *W2, double *W4,double *W6, double *W8,int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC,  double *k0, double *P0, double *errP0, int NeffP0, double *k2, double *P2, double *errP2, int NeffP2, double *k4, double *P4, double *errP4, int NeffP4, double *k11, double *k22, double *k33, double *B0, double *errB0, double *Bnoise, int NeffB0, double *k0SGC, double *P0SGC, double *errP0SGC,int NeffP0SGC, double *k2SGC, double *P2SGC, double *errP2SGC,int NeffP2SGC, double *k4SGC, double *P4SGC, double *errP4SGC,int NeffP4SGC, double *k11SGC, double *k22SGC, double *k33SGC,double *B0SGC, double *BnoiseSGC,int NeffB0SGC, double *cov, double *covSGC, double alpha_min, double alpha_max, char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks, char *path_output, char *identifier, char *do_plot, long int Nsteps, char *do_power_spectrum, char *do_bispectrum,int Nalphas,int Nsigmas_tot, int Nsigmas_free, double **Theory, double Pnoise, double PnoiseSGC, char *ptmodel_ps,  char *rsdmodel_ps,  char *fogmodel_ps,  char *ptmodel_bs,   char *local_b2s2,   char *local_b3nl,  char *RSD_fit,  char *sigma8_free,  char *fog_free,  char *fog_bs, double Sigma_smooth,int factor_sampling_mask,char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory )

void mcmc_kernel(int nthreads, char *type_BAO_fit,char *type_of_analysis,int trial_mcmc,double **vector_buffer,double fraction, int N_Cov_prop, gsl_matrix *transform, gsl_matrix *transform_inverse, double *vector_mean, char *name_file_output_mcmc, char *fit_BAO,char *fit_RSD,double *k_Plin,double *Plin,int N_Plin, double *k_Olin, double *Olin, int N_Olin, double *pos, double *W0, double *W2, double *W4,double *W6, double *W8,int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC,  double *k0bao, double *k0rsd, double *P0bao, double *P0rsd, double *errP0bao, double *errP0rsd, int NeffP0bao, int NeffP0rsd, double *k2bao, double *k2rsd, double *P2bao, double *P2rsd, double *errP2bao, double *errP2rsd, int NeffP2bao, int NeffP2rsd, double *k4bao, double *k4rsd, double *P4bao, double *P4rsd, double *errP4bao, double *errP4rsd, int NeffP4bao, int NeffP4rsd, double *k11bao, double *k11rsd, double *k22bao, double *k22rsd, double *k33bao, double *k33rsd, double *B0bao, double *B0rsd, double *errB0bao, double *errB0rsd, double *Bnoise_bao, double *Bnoise_rsd, int NeffB0bao, int NeffB0rsd, double *k0baoSGC, double *k0rsdSGC, double *P0baoSGC, double *P0rsdSGC, double *errP0baoSGC, double *errP0rsdSGC,int NeffP0baoSGC,int NeffP0rsdSGC, double *k2baoSGC,double *k2rsdSGC, double *P2baoSGC,double *P2rsdSGC, double *errP2baoSGC,double *errP2rsdSGC,int NeffP2baoSGC,int NeffP2rsdSGC, double *k4baoSGC, double *k4rsdSGC, double *P4baoSGC,double *P4rsdSGC, double *errP4baoSGC,double *errP4rsdSGC,int NeffP4baoSGC,int NeffP4rsdSGC, double *k11baoSGC,double *k11rsdSGC, double *k22baoSGC,double *k22rsdSGC, double *k33baoSGC,double *k33rsdSGC,double *B0baoSGC,double *B0rsdSGC, double *Bnoise_baoSGC, double *Bnoise_rsdSGC,int NeffB0baoSGC,int NeffB0rsdSGC, double *cov, double *covSGC, double alpha_min, double alpha_max, char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks, char *path_output, char *identifier, char *do_plot, long int Nsteps, char *do_power_spectrum, char *do_bispectrum,int Nalphas,int Nsigmas_tot, int Nsigmas_free, double **Theory, int Ntheory, double Pnoise, double PnoiseSGC, char *ptmodel_ps,  char *rsdmodel_ps,  char *fogmodel_ps,  char *ptmodel_bs,   char *local_b2s2,   char *local_b3nl,  char *RSD_fit,  char *sigma8_free,  char *fog_free,  char *fog_bs, double Sigma_smooth,int factor_sampling_mask,char *spacing_dataNGC_bao,char *spacing_dataNGC_rsd,char *spacing_dataSGC_bao,char *spacing_dataSGC_rsd,char *spacing_theory_bao,char *spacing_theory_rsd,char *type_of_analysis_BAO,char *type_of_analysis_FS )
{
int i_thread;
double time_run;
long int time_ini,time_final;
FILE *f;
double *Dwhite,*Doriginal,*step,*input_vector_trial,*input_vector_accepted;
int weight,convergence_private,convergence_shared,weight_print;
double normalize_step;
long int l,i,j,j_run_shared,j_run_private,counter,abs_counter_shared,abs_counter_private,i1,i2;
long int Nburnout,N_print,N_max,lines;
double L_current,L_proposed,acceptance,random_acceptance;
double **vector_print;
//int Nplan;
//double params[5];
int Nplanbao,Nplanrsd;
double paramsbao[5],paramsrsd[5];
double *priors_low;
double *priors_high;
double kmax_Nlin,kmin_Nlin;
double kmax_Ntheory,kmin_Ntheory;
double kmax_databao,kmin_databao;
double kmax_datarsd,kmin_datarsd;
long int *params_mcmc;
double *mean_vector, *min_vector;
double prior;
int dimension,dimensionbao,dimensionrsd;
int Ndof,Npoints,NpointsSGC;
int Npointsbao,Npointsrsd,NpointsbaoSGC,NpointsrsdSGC;
//double *parameters_plot;
double *parameters2;
double *parameters2_bao,*parameters2_rsd;
int modeP0bao,modeP2bao,modeP4bao;
int modeP0rsd,modeP2rsd,modeP4rsd;
char name_wo_extension[2000];
char name_file_thread[2000];
FILE *urandom;
unsigned short randstate[3];
double *seed;
int baoiso_shift;
//fftw_complex *a_pointer;
//fftw_complex *b_pointer;

fftw_complex *a_pointerbao;
fftw_complex *b_pointerbao;

fftw_complex *a_pointerrsd;
fftw_complex *b_pointerrsd;

fftw_plan plan1bao,plan1rsd,plan2bao,plan2rsd;

modeP0bao=0;modeP2bao=0;modeP4bao=0;Npointsbao=0;
if(strcmp(type_of_analysis_BAO,"yes")==0){
if(strcmp(fit_BAO, "P0") == 0){Npointsbao=NeffP0bao;modeP0bao=1;}
if(strcmp(fit_BAO, "P2") == 0){Npointsbao=NeffP2bao;modeP2bao=1;}
if(strcmp(fit_BAO, "P4") == 0){Npointsbao=NeffP4bao;modeP4bao=1;}
if(strcmp(fit_BAO, "P02") == 0){Npointsbao=NeffP0bao+NeffP2bao;modeP0bao=1;modeP2bao=1;}
if(strcmp(fit_BAO, "P04") == 0){Npointsbao=NeffP0bao+NeffP4bao;modeP0bao=1;modeP4bao=1;}
if(strcmp(fit_BAO, "P24") == 0){Npointsbao=NeffP2bao+NeffP4bao;modeP2bao=1;modeP4bao=1;}
if(strcmp(fit_BAO, "P024") == 0){Npointsbao=NeffP0bao+NeffP2bao+NeffP4bao;modeP0bao=1;modeP2bao=1;modeP4bao=1;}}

modeP0rsd=0;modeP2rsd=0;modeP4rsd=0;Npointsrsd=0;
if(strcmp(type_of_analysis_FS,"yes")==0){
if(strcmp(fit_RSD, "P0") == 0){Npointsrsd=NeffP0rsd;modeP0rsd=1;}
if(strcmp(fit_RSD, "P2") == 0){Npointsrsd=NeffP2rsd;modeP2rsd=1;}
if(strcmp(fit_RSD, "P4") == 0){Npointsrsd=NeffP4rsd;modeP4rsd=1;}
if(strcmp(fit_RSD, "P02") == 0){Npointsrsd=NeffP0rsd+NeffP2rsd;modeP0rsd=1;modeP2rsd=1;}
if(strcmp(fit_RSD, "P04") == 0){Npointsrsd=NeffP0rsd+NeffP4rsd;modeP0rsd=1;modeP4rsd=1;}
if(strcmp(fit_RSD, "P24") == 0){Npointsrsd=NeffP2rsd+NeffP4rsd;modeP2rsd=1;modeP4rsd=1;}
if(strcmp(fit_RSD, "P024") == 0){Npointsrsd=NeffP0rsd+NeffP2rsd+NeffP4rsd;modeP0rsd=1;modeP2rsd=1;modeP4rsd=1;}}
Npoints=Npointsrsd+Npointsbao;

NpointsSGC=0;
if(Nchunks==2)
{
NpointsbaoSGC=0;
if(strcmp(type_of_analysis_BAO,"yes")==0){
if(strcmp(fit_BAO, "P0") == 0){NpointsbaoSGC=NeffP0baoSGC;}
if(strcmp(fit_BAO, "P2") == 0){NpointsbaoSGC=NeffP2baoSGC;}
if(strcmp(fit_BAO, "P4") == 0){NpointsbaoSGC=NeffP4baoSGC;}
if(strcmp(fit_BAO, "P02") == 0){NpointsbaoSGC=NeffP0baoSGC+NeffP2baoSGC;}
if(strcmp(fit_BAO, "P04") == 0){NpointsbaoSGC=NeffP0baoSGC+NeffP4baoSGC;}
if(strcmp(fit_BAO, "P24") == 0){NpointsbaoSGC=NeffP2baoSGC+NeffP4baoSGC;}
if(strcmp(fit_BAO, "P024") == 0){NpointsbaoSGC=NeffP0baoSGC+NeffP2baoSGC+NeffP4baoSGC;}}

NpointsrsdSGC=0;
if(strcmp(type_of_analysis_FS,"yes")==0){
if(strcmp(fit_RSD, "P0") == 0){NpointsrsdSGC=NeffP0rsdSGC;}
if(strcmp(fit_RSD, "P2") == 0){NpointsrsdSGC=NeffP2rsdSGC;}
if(strcmp(fit_RSD, "P4") == 0){NpointsrsdSGC=NeffP4rsdSGC;}
if(strcmp(fit_RSD, "P02") == 0){NpointsrsdSGC=NeffP0rsdSGC+NeffP2rsdSGC;}
if(strcmp(fit_RSD, "P04") == 0){NpointsrsdSGC=NeffP0rsdSGC+NeffP4rsdSGC;}
if(strcmp(fit_RSD, "P24") == 0){NpointsrsdSGC=NeffP2rsdSGC+NeffP4rsdSGC;}
if(strcmp(fit_RSD, "P024") == 0){NpointsrsdSGC=NeffP0rsdSGC+NeffP2rsdSGC+NeffP4rsdSGC;}}

NpointsSGC=NpointsbaoSGC+NpointsrsdSGC;


}

//For BAO-only
if(strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0){
dimension=Nchunks*(Npolynomial+1)*(modeP0bao+modeP2bao+modeP4bao)+Nalphas+Nsigmas_tot;//for parameters plot
dimensionbao=Nchunks*(Npolynomial+1)*(modeP0bao+modeP2bao+modeP4bao)+Nalphas+Nsigmas_tot;//for parameters plot
//parameters2 =  (double*) calloc( dimension, sizeof(double));

if( strcmp(type_of_analysis, "FSBAOISO") == 0 && Nalphas==1){dimension=dimension+1;}//we add space for an extra alpha when FS-BAO and P0,P2,P4

//For Bispectrum to be included here

}

if(strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
dimension=Nchunks*(1+Npolynomial*(modeP0bao+modeP2bao+modeP4bao))+Nalphas+Nsigmas_tot+1;//for parameters plot
dimensionbao=Nchunks*(1+Npolynomial*(modeP0bao+modeP2bao+modeP4bao))+Nalphas+Nsigmas_tot+1;//for parameters plot
//parameters2 =  (double*) calloc( dimension, sizeof(double));
}

if(strcmp(type_of_analysis, "FS") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
//For RSD-only
dimension=3;//b1,b2,A
if(strcmp(local_b2s2, "no") == 0){dimension++;}//b2s2
if(strcmp(local_b3nl, "no") == 0){dimension++;}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){dimension++;}//fog_ps
if(strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 ){dimension++;}//fog_bs
dimension=dimension*Nchunks;//x2
if(strcmp(sigma8_free, "yes") == 0){dimension++;}//s8 as free parameter
if(strcmp(RSD_fit, "yes") == 0){dimension=dimension+3;}//f,apara,aperp as free parameter
dimensionrsd=dimension;
}

if(strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){dimension=dimensionrsd+dimensionbao-2;}

//printf("dimension %d\n",dimension);
//exit(0);
/*
 *  FFTW plans ARE NOT thread-safe and MUST BE created BEFORE the multi-thread starts
 */
//BAO
//
kmax_databao=-1;
kmin_databao=999999;
if(strcmp(type_of_analysis_BAO,"yes")==0){
if(modeP4bao==1)
{
if(kmin_databao>k4bao[0]){kmin_databao=k4bao[0];}
if(kmax_databao<k4bao[NeffP4bao-1]){kmax_databao=k4bao[NeffP4bao-1];}
if(Nchunks==2){
if(kmin_databao>k4baoSGC[0]){kmin_databao=k4baoSGC[0];}
if(kmax_databao<k4baoSGC[NeffP4baoSGC-1]){kmax_databao=k4baoSGC[NeffP4baoSGC-1];}}
}
if(modeP2bao==1)
{
if(kmin_databao>k2bao[0]){kmin_databao=k2bao[0];}
if(kmax_databao<k2bao[NeffP2bao-1]){kmax_databao=k2bao[NeffP2bao-1];}
if(Nchunks==2){
if(kmin_databao>k2baoSGC[0]){kmin_databao=k2baoSGC[0];}
if(kmax_databao<k2baoSGC[NeffP2baoSGC-1]){kmax_databao=k2baoSGC[NeffP2baoSGC-1];}}
}
if(modeP0bao==1)
{
if(kmin_databao>k0bao[0]){kmin_databao=k0bao[0];}
if(kmax_databao<k0bao[NeffP0bao-1]){kmax_databao=k0bao[NeffP0bao-1];}
if(Nchunks==2){
if(kmin_databao>k0baoSGC[0]){kmin_databao=k0baoSGC[0];}
if(kmax_databao<k0baoSGC[NeffP0baoSGC-1]){kmax_databao=k0baoSGC[NeffP0baoSGC-1];}}
}
if(kmin_databao==999999 || kmax_databao==-1){printf("Error with values of kmin,kmax for data-BAO. Exiting now...\n");exit(0);}
}


//RSD
kmax_datarsd=-1;
kmin_datarsd=999999;
if(strcmp(type_of_analysis_FS,"yes")==0){
if(modeP4rsd==1)
{
if(kmin_datarsd>k4rsd[0]){kmin_datarsd=k4rsd[0];}
if(kmax_datarsd<k4rsd[NeffP4rsd-1]){kmax_datarsd=k4rsd[NeffP4rsd-1];}
if(Nchunks==2){
if(kmin_datarsd>k4rsdSGC[0]){kmin_datarsd=k4rsdSGC[0];}
if(kmax_datarsd<k4rsdSGC[NeffP4rsdSGC-1]){kmax_datarsd=k4rsdSGC[NeffP4rsdSGC-1];}}
}
if(modeP2rsd==1)
{
if(kmin_datarsd>k2rsd[0]){kmin_datarsd=k2rsd[0];}
if(kmax_datarsd<k2rsd[NeffP2rsd-1]){kmax_datarsd=k2rsd[NeffP2rsd-1];}
if(Nchunks==2){
if(kmin_datarsd>k2rsdSGC[0]){kmin_datarsd=k2rsdSGC[0];}
if(kmax_datarsd<k2rsdSGC[NeffP2rsdSGC-1]){kmax_datarsd=k2rsdSGC[NeffP2rsdSGC-1];}}
}
if(modeP0rsd==1)
{
if(kmin_datarsd>k0rsd[0]){kmin_datarsd=k0rsd[0];}
if(kmax_datarsd<k0rsd[NeffP0rsd-1]){kmax_datarsd=k0rsd[NeffP0rsd-1];}
if(Nchunks==2){
if(kmin_datarsd>k0rsdSGC[0]){kmin_datarsd=k0rsdSGC[0];}
if(kmax_datarsd<k0rsdSGC[NeffP0rsdSGC-1]){kmax_datarsd=k0rsdSGC[NeffP0rsdSGC-1];}}
}
if(kmin_datarsd==999999 || kmax_datarsd==-1){printf("Error with values of kmin,kmax for data-RSD. Exiting now...\n");exit(0);}
}

if( strcmp(type_of_analysis, "BAOISO") == 0 ||  strcmp(type_of_analysis, "BAOANISO") == 0){
kmin_Nlin=k_Plin[0];
kmax_Nlin=k_Plin[N_Plin-1];
set_mask_params(paramsbao,kmin_Nlin,kmax_Nlin,N_Plin,kmin_databao,kmax_databao);

Nplanbao=(int)(paramsbao[2]);

    /*fftw_plan*/ plan1bao = fftw_plan_dft_1d(Nplanbao,  a_pointerbao,  b_pointerbao,  -1, FFTW_ESTIMATE);//forward plan
    /*fftw_plan*/ plan2bao = fftw_plan_dft_1d(Nplanbao,  b_pointerbao,  b_pointerbao, +1, FFTW_ESTIMATE);//reverse plan

}

if(  strcmp(type_of_analysis, "FS") == 0 ){
kmin_Ntheory=Theory[0][0];
kmax_Ntheory=Theory[Ntheory-1][0];
set_mask_params(paramsrsd,kmin_Ntheory,kmax_Ntheory,Ntheory,kmin_datarsd,kmax_datarsd);

Nplanrsd=(int)(paramsrsd[2]);

    /*fftw_plan*/ plan1rsd = fftw_plan_dft_1d(Nplanrsd,  a_pointerrsd,  b_pointerrsd,  -1, FFTW_ESTIMATE);//forward plan
    /*fftw_plan*/ plan2rsd = fftw_plan_dft_1d(Nplanrsd,  b_pointerrsd,  b_pointerrsd, +1, FFTW_ESTIMATE);//reverse plan


}

if( strcmp(type_of_analysis, "FSBAOISO") == 0 ||  strcmp(type_of_analysis, "FSBAOANISO") == 0){
kmin_Nlin=k_Plin[0];
kmax_Nlin=k_Plin[N_Plin-1];
kmin_Ntheory=Theory[0][0];
kmax_Ntheory=Theory[Ntheory-1][0];

//if(kmax_Nlin<kmax_Ntheory){set_mask_params(params,kmin_Ntheory,kmax_Ntheory,Ntheory,kmin_datarsd,kmax_datarsd);}
//else{set_mask_params(params,kmin_Nlin,kmax_Nlin,N_Plin,kmin_databao,kmax_databao);}
set_mask_params(paramsrsd,kmin_Ntheory,kmax_Ntheory,Ntheory,kmin_datarsd,kmax_datarsd);
set_mask_params(paramsbao,kmin_Nlin,kmax_Nlin,N_Plin,kmin_databao,kmax_databao);

Nplanbao=(int)(paramsbao[2]);

    /*fftw_plan*/ plan1bao = fftw_plan_dft_1d(Nplanbao,  a_pointerbao,  b_pointerbao,  -1, FFTW_ESTIMATE);//forward plan
    /*fftw_plan*/ plan2bao = fftw_plan_dft_1d(Nplanbao,  b_pointerbao,  b_pointerbao, +1, FFTW_ESTIMATE);//reverse plan

Nplanrsd=(int)(paramsrsd[2]);

    /*fftw_plan*/ plan1rsd = fftw_plan_dft_1d(Nplanrsd,  a_pointerrsd,  b_pointerrsd,  -1, FFTW_ESTIMATE);//forward plan
    /*fftw_plan*/ plan2rsd = fftw_plan_dft_1d(Nplanrsd,  b_pointerrsd,  b_pointerrsd, +1, FFTW_ESTIMATE);//reverse plan


}


priors_low=(double*) calloc(N_Cov_prop, sizeof(double));
priors_high=(double*) calloc(N_Cov_prop, sizeof(double));

set_mcmc_priors(alpha_min,alpha_max,priors_low,priors_high,N_Cov_prop);

for(i=0;i<N_Cov_prop;i++)
{
if(priors_low[i]>=priors_high[i]){printf("Error. Non-phyiscal prior range: prior_low[%ld] (%lf) >=prior_high[%ld] (%lf). Exiting now...\n",i,priors_low[i],i,priors_high[i]);printf("Ncov=%d\n",N_Cov_prop);exit(0);}
}


params_mcmc=(long int*) calloc( 2, sizeof(long int ));
set_mcmc_parameters(params_mcmc);
N_print=params_mcmc[0];
Nburnout=params_mcmc[1];
free(params_mcmc);

if(Nsteps<N_print && trial_mcmc==0){printf("Warning: Unusual low value of maximum number of steps (%ld) with relation to printing buffer (%ld). Exiting now...\n",Nsteps,N_print);exit(0);}

N_max=(long int)(Nsteps*fraction);
Nburnout=(long int)(Nburnout*fraction);
abs_counter_shared = 0;
j_run_shared = 0;
convergence_shared = 0;
if(trial_mcmc==1){normalize_step=1.9/sqrt(N_Cov_prop*1.);}//fraction of the step wrt the typical rms
else{normalize_step=1.9/sqrt(1.*N_Cov_prop);}

time_ini=time(NULL);

srand48(time(NULL)*123456789);
//srand48(123456789);
seed = (double *) calloc(nthreads*3,sizeof(double));

for(i=0;i<nthreads*3;i++)seed[i]=drand48();

#pragma omp parallel for \
  private (urandom,randstate,i_thread,  f, Dwhite, Doriginal, step, parameters2,parameters2_bao,parameters2_rsd, input_vector_trial, \
           input_vector_accepted, weight, convergence_private, weight_print, \
           l, i, j, j_run_private, counter, abs_counter_private, i1, i2,baoiso_shift, vector_print, \
           L_current, L_proposed, acceptance, random_acceptance, prior, name_file_thread) \
  shared (plan1bao,plan2bao,plan1rsd,plan2rsd,seed,nthreads,type_BAO_fit, type_of_analysis, trial_mcmc, N_Cov_prop, transform, transform_inverse, \
          name_file_output_mcmc, fit_BAO,fit_RSD, k_Plin, Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4, W6, \
          W8, Nmask, path_to_mask1, spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, \
          path_to_mask2, spacing_maskSGC, k0bao,k0rsd, P0bao,P0rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, NeffP4bao,NeffP4rsd, \
          k11bao,k11rsd, k22bao,k22rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, \
          NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC, B0baoSGC,B0rsdSGC, Bnoise_baoSGC,Bnoise_rsdSGC, NeffB0baoSGC,NeffB0rsdSGC, \
          cov, covSGC, Sigma_def_type, Sigma_independent, ffactor, Sigma_type, \
          Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, Nsteps, \
          do_power_spectrum, do_bispectrum, Nalphas, Nsigmas_tot, Nsigmas_free, Theory, Pnoise, PnoiseSGC, \
          ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl, RSD_fit, sigma8_free, fog_free, \
          fog_bs, Sigma_smooth, factor_sampling_mask, spacing_dataNGC_rsd,spacing_dataNGC_bao, spacing_dataSGC_bao,spacing_dataSGC_rsd, spacing_theory_bao,spacing_theory_rsd, \
          vector_buffer, vector_mean, abs_counter_shared, convergence_shared, j_run_shared, \
          Nburnout, N_print, N_max, priors_low, priors_high, dimension,dimensionbao,dimensionrsd, \
          modeP0bao,modeP0rsd, modeP2bao,modeP2rsd, modeP4bao,modeP4rsd, name_wo_extension, normalize_step, lines,type_of_analysis_BAO,type_of_analysis_FS) \
  schedule (static,1) \
  num_threads (nthreads)

for (i_thread=0;i_thread<nthreads;i_thread++)
{

printf("Starting mcmc chain with index %d on thread %d\n", i_thread, omp_get_thread_num());


Dwhite = (double*) calloc( N_Cov_prop, sizeof(double));
Doriginal = (double*) calloc( N_Cov_prop, sizeof(double));
step = (double*) calloc( N_Cov_prop, sizeof(double));
input_vector_trial = (double*) calloc( N_Cov_prop, sizeof(double));
input_vector_accepted = (double*) calloc( N_Cov_prop, sizeof(double));
parameters2 =  (double*) calloc( dimension, sizeof(double));

if(strcmp(type_of_analysis_BAO,"yes")==0 && strcmp(type_of_analysis_FS,"yes")==0)
{
parameters2_bao =  (double*) calloc( dimensionbao, sizeof(double));
parameters2_rsd =  (double*) calloc( dimensionrsd, sizeof(double));
}

//if(trial_mcmc==1){normalize_step=1.0;}//fraction of the step wrt the typical rms
//else{normalize_step=1.0;}


//srand48(time(NULL)*123456789*(i_thread+1));//HGM: I have included here i_thread as a source for the random seed
//srand48_r(time(NULL)*123456789*(i_thread+1));//HGM: I have included here i_thread as a source for the random seed

   
      urandom = fopen ("/dev/urandom", "r");
      if(urandom==NULL){
      printf("Warning, /dev/urandom fail when reading..., getting the seed from drand48 with time as inizialization...\n");
randstate[0]=(unsigned short)(seed[i_thread*3]*65535);
randstate[1]=(unsigned short)(seed[i_thread*3+1]*65535);
randstate[2]=(unsigned short)(seed[i_thread*3+2]*65535);
}
      else{
      setvbuf (urandom, NULL, _IONBF, 0);  // turn off buffering
    
       randstate[0] = (fgetc (urandom) << 8) | fgetc (urandom);
       randstate[1] = (fgetc (urandom) << 8) | fgetc (urandom);
       randstate[2] = (fgetc (urandom) << 8) | fgetc (urandom);
       fclose (urandom);
     }

/*
printf("i_thread=%d %hu %hu %hu\n",i_thread,randstate[0],randstate[1],randstate[2]);
random_acceptance=erand48(randstate);
printf("%d, A %lf\n",i_thread,random_acceptance);
random_acceptance=erand48(randstate);
printf("%d, B %lf\n",i_thread,random_acceptance);
*/

//initial condiction
for(i=0;i<N_Cov_prop;i++)
{
input_vector_accepted[i]=vector_mean[i];
}

//for BAO
if(strcmp(type_of_analysis, "BAOISO") == 0){

parameters2[0]=input_vector_accepted[0];i1=1;i2=1;
if(modeP0bao+modeP2bao+modeP4bao>1){parameters2[1]=input_vector_accepted[1];i1=2;i2=2;}
if(Nsigmas_free==0)
{
if(modeP0bao+modeP2bao+modeP4bao==1)
{
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
  parameters2[1]=Sigma_nl_mean[0];
  parameters2[2]=Sigma_nl_mean[1];
  i1=3;i2=2;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0 ){parameters2[1]=Sigma_nl_mean[0];i1=2;}
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
    if(modeP0bao==1){parameters2[1]=Sigma_nl_mean[0];}
    if(modeP2bao==1){parameters2[1]=Sigma_nl_mean[1];}
    if(modeP4bao==1){parameters2[1]=Sigma_nl_mean[2];}
    i1=2;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
   parameters2[1]=Sigma_nl_mean[0];
   i1=2;
  }
}
if(modeP0bao+modeP2bao+modeP4bao>1){
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=Sigma_nl_mean[1];
  i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];parameters2[4]=Sigma_nl_mean[2];i1=5;}
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==0){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];i1=4;}
    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[2];i1=4;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[1];parameters2[3]=Sigma_nl_mean[2];i1=4;}
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0){parameters2[2]=Sigma_nl_mean[0];i1=3;}
}
}
//for(i=Nalphas+Nsigmas_tot-Nsigmas_free;i<N_Cov_prop;i++){parameters2[i+Nsigmas_tot-Nsigmas_free]=input_vector_accepted[i];}//set the rest of parameters
for(i=i1;i<dimension;i++){parameters2[i]=input_vector_accepted[i2];i2++;}//set the rest of parameters


L_current=chi2_bao(type_BAO_fit,type_of_analysis, fit_BAO ,parameters2, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao, P0bao,NeffP0bao, k2bao, P2bao,NeffP2bao, k4bao, P4bao,NeffP4bao, k0baoSGC, P0baoSGC,NeffP0baoSGC, k2baoSGC, P2baoSGC,NeffP2baoSGC, k4baoSGC, P4baoSGC, NeffP4baoSGC, cov, covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, Nchunks, plan1bao, plan2bao, do_power_spectrum, do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);
//printf("chi2-baoiso kernel=%lf\n",L_current);

//exit(0);
}

if(strcmp(type_of_analysis, "BAOANISO") == 0){

parameters2[0]=input_vector_accepted[0];
parameters2[1]=input_vector_accepted[1];
i1=2;i2=2;
if(Nsigmas_free==0)
{

  if(strcmp(Sigma_independent, "yes") == 0)//Nsigma_tot=2
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=Sigma_nl_mean[1];
  parameters2[4]=input_vector_accepted[2];
  i1=5;i2=3;

  }
  if(strcmp(Sigma_independent, "no") == 0)//Nsigma_tot=1
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=input_vector_accepted[2];
  i1=4;i2=3;
  }
 
}

for(i=i1;i<dimension;i++){parameters2[i]=input_vector_accepted[i2];i2++;}//set the rest of parameters
//for(i=0;i<dimension;i++){printf("%ld, %lf\n",i,parameters2[i]);}

//exit(0);
L_current=chi2_bao(type_BAO_fit,type_of_analysis, fit_BAO ,parameters2, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao, P0bao,NeffP0bao, k2bao, P2bao,NeffP2bao, k4bao, P4bao,NeffP4bao, k0baoSGC, P0baoSGC,NeffP0baoSGC, k2baoSGC, P2baoSGC,NeffP2baoSGC, k4baoSGC, P4baoSGC, NeffP4baoSGC, cov, covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, Nchunks, plan1bao, plan2bao, do_power_spectrum, do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);
//printf("chi2-bao aniso kernel=%lf\n",L_current);
//exit(0);
}
//exit(0);

if(strcmp(type_of_analysis, "FS") == 0){

//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",input_vector_accepted[0],input_vector_accepted[1],input_vector_accepted[2],input_vector_accepted[3],input_vector_accepted[4],input_vector_accepted[5],input_vector_accepted[6],input_vector_accepted[7],input_vector_accepted[8],input_vector_accepted[9],input_vector_accepted[10]);

//exit(0);
//L_current=chi2_rsd();

L_current=chi2_rsd_mcmc(type_of_analysis,input_vector_accepted, Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0rsd, P0rsd, Pnoise, k2rsd, P2rsd, k4rsd, P4rsd, k0rsdSGC, P0rsdSGC, PnoiseSGC, k2rsdSGC, P2rsdSGC, k4rsdSGC, P4rsdSGC, NeffP0rsd, NeffP2rsd, NeffP4rsd,NeffP0rsdSGC,NeffP2rsdSGC,NeffP4rsdSGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1rsd, plan2rsd, fit_RSD, do_power_spectrum, do_bispectrum,factor_sampling_mask,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd);

//if(i_thread==0){printf("chi2 FS kernel=%lf (%d)\n",L_current,i_thread);}
//printf("chi2 FS kernel=%lf\n",L_current);
//exit(0);

}

if(strcmp(type_of_analysis, "FSBAOISO") == 0){

parameters2[0]=input_vector_accepted[0];//apara
parameters2[1]=input_vector_accepted[1];//aperp
parameters2[2]=input_vector_accepted[2];//f

//parameters2_bao[0]=input_vector_accepted[0];//apara
//parameters2_bao[1]=input_vector_accepted[1];//aperp

if(modeP0bao+modeP2bao+modeP4bao>1){
parameters2_bao[0]=input_vector_accepted[0];//apara
parameters2_bao[1]=input_vector_accepted[1];//aperp
baoiso_shift=1;
}
else
{
baoiso_shift=2;
if(modeP0bao==1){parameters2_bao[0]=pow(input_vector_accepted[0],1./3.)*pow(input_vector_accepted[1],2./3.);  }
if(modeP2bao==1){parameters2_bao[0]=pow(input_vector_accepted[0],3./5.)*pow(input_vector_accepted[1],2./5.);  }
if(modeP4bao==1){parameters2_bao[0]=pow(input_vector_accepted[0],5./7.)*pow(input_vector_accepted[1],2./7.);  }
}

i1=2;i2=3;
if(Nsigmas_free==0)
{

if(modeP0bao+modeP2bao+modeP4bao==1)
{

  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
  parameters2_bao[1]=Sigma_nl_mean[0];
  parameters2_bao[2]=Sigma_nl_mean[1];

  parameters2[3]=Sigma_nl_mean[0];
  parameters2[4]=Sigma_nl_mean[1];

  i1=4;//??
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
     parameters2[3]=Sigma_nl_mean[0];
     parameters2_bao[1]=Sigma_nl_mean[0];
     i1=3;//??

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  { 
    if(modeP0bao==1){parameters2_bao[1]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[0];}
    if(modeP2bao==1){parameters2_bao[1]=Sigma_nl_mean[1];parameters2[3]=Sigma_nl_mean[1];}
    if(modeP4bao==1){parameters2_bao[1]=Sigma_nl_mean[2];parameters2[3]=Sigma_nl_mean[2];}
    i1=3;//??
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
   parameters2[3]=Sigma_nl_mean[0];
   parameters2_bao[1]=Sigma_nl_mean[0];
   i1=3;//??
  }


}
if(modeP0bao+modeP2bao+modeP4bao>1){

  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2_bao[2]=Sigma_nl_mean[0];
  parameters2_bao[3]=Sigma_nl_mean[1];
  parameters2[3]=Sigma_nl_mean[0];
  parameters2[4]=Sigma_nl_mean[1];

  i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[3]=Sigma_nl_mean[0];
  parameters2_bao[2]=Sigma_nl_mean[0];

  i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[1];
         parameters2_bao[4]=Sigma_nl_mean[2];   
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[1];
         parameters2[5]=Sigma_nl_mean[2];
         i1=5;}

    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[1];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[1];
         i1=4;}

    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[2];
         i1=4;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[1];
         parameters2_bao[3]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[1];
         parameters2[4]=Sigma_nl_mean[2];
         i1=4;}
    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2[3]=Sigma_nl_mean[0];
         i1=3;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[1];
         parameters2[3]=Sigma_nl_mean[1];
         i1=3;}
    if(modeP0bao==0 && modeP2bao==0 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[2];
         i1=3;}

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0){
          parameters2[3]=Sigma_nl_mean[0];parameters2_bao[2]=Sigma_nl_mean[0];
          i1=3;}
}
}

for(i=i1+1;i<dimensionbao+baoiso_shift;i++){parameters2[i]=input_vector_accepted[i2];parameters2_bao[i-baoiso_shift]=input_vector_accepted[i2];i2++;}//set the rest of parameters (check we need dimension+1)

parameters2_rsd[0]=input_vector_accepted[0];
parameters2_rsd[1]=input_vector_accepted[1];
parameters2_rsd[2]=input_vector_accepted[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=input_vector_accepted[i2];parameters2[i-3+dimensionbao+1]=input_vector_accepted[i2];i2++;}

L_current=chi2_bao_rsd(type_BAO_fit,type_of_analysis,fit_BAO,fit_RSD,parameters2_bao,parameters2_rsd, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin,Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao,k0rsd, P0bao,P0rsd,Pnoise,NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd,NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd,NeffP4bao,NeffP4rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC,PnoiseSGC,NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC,NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC,cov,covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1bao, plan2bao, plan1rsd, plan2rsd, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory_bao,spacing_theory_rsd);

//for(i=0;i<dimensionbao;i++){printf("BAO-%d=%lf\n",i,parameters2_bao[i]);}
//for(i=0;i<dimensionrsd;i++){printf("RSD-%d=%lf\n",i,parameters2_rsd[i]);}
//for(i=0;i<dimension;i++){printf("BAO+RSD-%d=%lf\n",i,parameters2[i]);}
//exit(0);
//printf("FS-BAOISO %lf\n",L_current);

}

if(strcmp(type_of_analysis, "FSBAOANISO") == 0){

parameters2[0]=input_vector_accepted[0];//apara
parameters2[1]=input_vector_accepted[1];//aperp
parameters2[2]=input_vector_accepted[2];//f

parameters2_bao[0]=input_vector_accepted[0];//apara
parameters2_bao[1]=input_vector_accepted[1];//aperp

i1=2;i2=3;
if(Nsigmas_free==0)
{

  if(strcmp(Sigma_independent, "yes") == 0)//Nsigma_tot=2
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2_bao[3]=Sigma_nl_mean[1];//sigmaperp
  parameters2[3]=Sigma_nl_mean[0];//sigmapara
  parameters2[4]=Sigma_nl_mean[1];//sigmaperp

  i1=4;i2=3;

  }
  if(strcmp(Sigma_independent, "no") == 0)//Nsigma_tot=1
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2[3]=Sigma_nl_mean[0];//sigmapara

  i1=3;i2=3;
  }

}

for(i=i1+1;i<dimensionbao+1;i++){parameters2[i]=input_vector_accepted[i2];parameters2_bao[i-1]=input_vector_accepted[i2];i2++;}//set the rest of parameters (check we need dimension+1)

parameters2_rsd[0]=input_vector_accepted[0];
parameters2_rsd[1]=input_vector_accepted[1];
parameters2_rsd[2]=input_vector_accepted[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=input_vector_accepted[i2];parameters2[i-3+dimensionbao+1]=input_vector_accepted[i2];i2++;}

//for(i=0;i<dimensionbao;i++){printf("BAO-%d=%lf\n",i,parameters2_bao[i]);}
//for(i=0;i<dimensionrsd;i++){printf("RSD-%d=%lf\n",i,parameters2_rsd[i]);}
//for(i=0;i<dimension;i++){printf("BAO+RSD-%d=%lf\n",i,parameters2[i]);}


L_current=chi2_bao_rsd(type_BAO_fit,type_of_analysis,fit_BAO,fit_RSD,parameters2_bao,parameters2_rsd, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin,Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao,k0rsd, P0bao,P0rsd,Pnoise,NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd,NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd,NeffP4bao,NeffP4rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC,PnoiseSGC,NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC,NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC,cov,covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1bao, plan2bao, plan1rsd, plan2rsd, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory_bao,spacing_theory_rsd);

//printf("FS-BAOANISO %lf\n",//L_current);

}

if(trial_mcmc==0)
{
vector_print = (double**) calloc(N_Cov_prop+2,sizeof(double*));

for(j=0;j<N_Cov_prop+2;j++)
{
   vector_print[j] = (double*) calloc(N_print,sizeof(double));

}

}

weight=1;
j_run_private=-1;
counter=0;//printing counter
abs_counter_private=0;//absolute counter
convergence_private=0;

do
{
j_run_private++;

if(j_run_private==0)
{
//for(i=0;i<N_Cov_prop;i++){step[i]=0.5*drand48();}
for(i=0;i<N_Cov_prop;i++){step[i]=0.5*erand48(randstate);}

}
else
{
//for(i=0;i<N_Cov_prop;i++){step[i]=drand48();}
for(i=0;i<N_Cov_prop;i++){step[i]=erand48(randstate);}
}
//initialize_vector(Doriginal);
for(i=0;i<N_Cov_prop;i++)
{
Doriginal[i]=input_vector_accepted[i];
//printf("%ld %lf\n",i,Doriginal[i]);
}
      //Rotation
       for(i=0;i<N_Cov_prop;i++)
       {
           Dwhite[i]=0;
           for(j=0;j<N_Cov_prop;j++)
           {
                Dwhite[i]=Dwhite[i]+(Doriginal[j]-vector_mean[j])*(gsl_matrix_get(transform_inverse, i, j));
  //              printf("%ld %ld %lf\n",i,j,gsl_matrix_get(transform_inverse, i, j));
           }

       }
 //make step
   for(i=0;i<N_Cov_prop;i++){Dwhite[i]=Dwhite[i]+(step[i]*2-1.)*normalize_step;}

 //Rotate back
 for(i=0;i<N_Cov_prop;i++)
       {
           input_vector_trial[i]=0;
           for(j=0;j<N_Cov_prop;j++)
           {
                input_vector_trial[i]=input_vector_trial[i]+Dwhite[j]*(gsl_matrix_get(transform, i, j));
           }
                input_vector_trial[i]=input_vector_trial[i]+vector_mean[i];
       }

if(strcmp(type_of_analysis, "BAOISO") == 0){

parameters2[0]=input_vector_trial[0];i1=1;i2=1;
if(modeP0bao+modeP2bao+modeP4bao>1){parameters2[1]=input_vector_trial[1];i1=2;i2=2;}
if(Nsigmas_free==0)
{
if(modeP0bao+modeP2bao+modeP4bao==1)
{
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
  parameters2[1]=Sigma_nl_mean[0];
  parameters2[2]=Sigma_nl_mean[1];i1=3;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0 ){parameters2[1]=Sigma_nl_mean[0];i1=2;}
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
    if(modeP0bao==1){parameters2[1]=Sigma_nl_mean[0];}
    if(modeP2bao==1){parameters2[1]=Sigma_nl_mean[1];}
    if(modeP4bao==1){parameters2[1]=Sigma_nl_mean[2];}
    i1=2;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
   parameters2[1]=Sigma_nl_mean[0];i1=2;
  }
}
if(modeP0bao+modeP2bao+modeP4bao>1){
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=Sigma_nl_mean[1];
  i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];parameters2[4]=Sigma_nl_mean[2];i1=5;}
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==0){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];i1=4;}
    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[2];i1=4;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[1];parameters2[3]=Sigma_nl_mean[2];i1=4;}
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0){parameters2[2]=Sigma_nl_mean[0];i1=3;}
}
}

//for(i=Nalphas+Nsigmas_tot-Nsigmas_free;i<N_Cov_prop;i++){parameters2[i+Nsigmas_tot-Nsigmas_free]=input_vector_trial[i];}//set the rest of parameters
for(i=i1;i<dimension;i++){parameters2[i]=input_vector_trial[i2];i2++;}//set the rest of parameters

//for(i=0;i<dimension;i++){printf("%d, %lf\n",i,parameters2[i]);}

L_proposed=chi2_bao(type_BAO_fit,type_of_analysis, fit_BAO ,parameters2, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao, P0bao,NeffP0bao, k2bao, P2bao,NeffP2bao, k4bao, P4bao,NeffP4bao, k0baoSGC, P0baoSGC,NeffP0baoSGC, k2baoSGC, P2baoSGC,NeffP2baoSGC, k4baoSGC, P4baoSGC, NeffP4baoSGC, cov, covSGC,  Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, Nchunks, plan1bao, plan2bao, do_power_spectrum, do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);
//printf("chi2-baoiso =%lf, %d\n",L_proposed,abs_counter_private);

//exit(0);
/*
for(i=0;i<Nsigmas_tot;i++)
{
if(Sigma_type[i]==1){printf("prior2 %lf (%lf,%lf,%lf)\n",gauss(parameters2[Nalphas+i],Sigma_nl_mean[i],Sigma_nl_stddev[i]),parameters2[Nalphas+i],Sigma_nl_mean[i],Sigma_nl_stddev[i]);}
}
*/

}

if(strcmp(type_of_analysis, "BAOANISO") == 0){

parameters2[0]=input_vector_trial[0];
parameters2[1]=input_vector_trial[1];
i1=2;i2=2;
if(Nsigmas_free==0)
{

  if(strcmp(Sigma_independent, "yes") == 0)//Nsigma_tot=2
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=Sigma_nl_mean[1];
  parameters2[4]=input_vector_trial[2];
  i1=5;i2=3;

  }
  if(strcmp(Sigma_independent, "no") == 0)//Nsigma_tot=1
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=input_vector_trial[2];
  i1=4;i2=3;
  }

}

for(i=i1;i<dimension;i++){parameters2[i]=input_vector_trial[i2];i2++;}//set the rest of parameters

//;for(i=0;i<dimension;i++){printf("%ld, %lf\n",i,parameters2[i]);}

L_proposed=chi2_bao(type_BAO_fit,type_of_analysis, fit_BAO ,parameters2, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao, P0bao,NeffP0bao, k2bao, P2bao,NeffP2bao, k4bao, P4bao,NeffP4bao, k0baoSGC, P0baoSGC,NeffP0baoSGC, k2baoSGC, P2baoSGC,NeffP2baoSGC, k4baoSGC, P4baoSGC, NeffP4baoSGC, cov, covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, Nchunks, plan1bao, plan2bao, do_power_spectrum, do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);
//printf("chi2-baoaniso=%lf, %d\n",L_proposed,abs_counter_private);
//exit(0);
}



if(strcmp(type_of_analysis, "FS") == 0){
//L_proposed=chi2_rsd();

//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",input_vector_trial[0],input_vector_trial[1],input_vector_trial[2],input_vector_trial[3],input_vector_trial[4],input_vector_trial[5],input_vector_trial[6],input_vector_trial[7],input_vector_trial[8],input_vector_trial[9],input_vector_trial[10]);

L_proposed=chi2_rsd_mcmc(type_of_analysis,input_vector_trial, Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0rsd, P0rsd, Pnoise, k2rsd, P2rsd, k4rsd, P4rsd, k0rsdSGC, P0rsdSGC, PnoiseSGC, k2rsdSGC, P2rsdSGC, k4rsdSGC, P4rsdSGC, NeffP0rsd, NeffP2rsd, NeffP4rsd,NeffP0rsdSGC,NeffP2rsdSGC,NeffP4rsdSGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1rsd, plan2rsd, fit_RSD, do_power_spectrum, do_bispectrum,factor_sampling_mask,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd);

//if(i_thread==0){printf("chi2 FS kernel=%lf (%d)\n",L_proposed,i_thread);}
//printf("chi2-FS=%lf\n",L_proposed);

}

if(strcmp(type_of_analysis, "FSBAOISO") == 0){

parameters2[0]=input_vector_trial[0];//apara
parameters2[1]=input_vector_trial[1];//aperp
parameters2[2]=input_vector_trial[2];//f

if(modeP0bao+modeP2bao+modeP4bao>1){
parameters2_bao[0]=input_vector_trial[0];//apara
parameters2_bao[1]=input_vector_trial[1];//aperp
baoiso_shift=1;
}
else
{
baoiso_shift=2;
if(modeP0bao==1){parameters2_bao[0]=pow(input_vector_trial[0],1./3.)*pow(input_vector_trial[1],2./3.);  }
if(modeP2bao==1){parameters2_bao[0]=pow(input_vector_trial[0],3./5.)*pow(input_vector_trial[1],2./5.);  }
if(modeP4bao==1){parameters2_bao[0]=pow(input_vector_trial[0],5./7.)*pow(input_vector_trial[1],2./7.);  }
}

i1=2;i2=3;

if(Nsigmas_free==0)
{

if(modeP0bao+modeP2bao+modeP4bao==1)
{

  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
  parameters2_bao[1]=Sigma_nl_mean[0];
  parameters2_bao[2]=Sigma_nl_mean[1];

  parameters2[3]=Sigma_nl_mean[0];
  parameters2[4]=Sigma_nl_mean[1];

  i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
     parameters2[3]=Sigma_nl_mean[0];
     parameters2_bao[1]=Sigma_nl_mean[0];
     i1=3;

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
    if(modeP0bao==1){parameters2_bao[1]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[0];}
    if(modeP2bao==1){parameters2_bao[1]=Sigma_nl_mean[1];parameters2[3]=Sigma_nl_mean[1];}
    if(modeP4bao==1){parameters2_bao[1]=Sigma_nl_mean[2];parameters2[3]=Sigma_nl_mean[2];}
    i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
   parameters2[3]=Sigma_nl_mean[0];
   parameters2_bao[1]=Sigma_nl_mean[0];
   i1=3;
  }


}
if(modeP0bao+modeP2bao+modeP4bao>1){




  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2_bao[2]=Sigma_nl_mean[0];
  parameters2_bao[3]=Sigma_nl_mean[1];
  parameters2[3]=Sigma_nl_mean[0];
  parameters2[4]=Sigma_nl_mean[1];

  i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[3]=Sigma_nl_mean[0];
  parameters2_bao[2]=Sigma_nl_mean[0];

  i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[1];
         parameters2_bao[4]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[1];
         parameters2[5]=Sigma_nl_mean[2];
         i1=5;}
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[1];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[1];
         i1=4;}

    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[2];
         i1=4;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[1];
         parameters2_bao[3]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[1];
         parameters2[4]=Sigma_nl_mean[2];
         i1=4;}
    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2[3]=Sigma_nl_mean[0];
         i1=3;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[1];
         parameters2[3]=Sigma_nl_mean[1];
         i1=3;}
    if(modeP0bao==0 && modeP2bao==0 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[2];
         i1=3;}

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0){
          parameters2[3]=Sigma_nl_mean[0];parameters2_bao[2]=Sigma_nl_mean[0];
          i1=3;}
}
}
for(i=i1+1;i<dimensionbao+baoiso_shift;i++){parameters2[i]=input_vector_trial[i2];parameters2_bao[i-baoiso_shift]=input_vector_trial[i2];i2++;}//set the rest of parameters (check we need dimension+1)

parameters2_rsd[0]=input_vector_trial[0];
parameters2_rsd[1]=input_vector_trial[1];
parameters2_rsd[2]=input_vector_trial[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=input_vector_trial[i2];parameters2[i-3+dimensionbao+1]=input_vector_trial[i2];i2++;}


L_proposed=chi2_bao_rsd(type_BAO_fit,type_of_analysis,fit_BAO,fit_RSD,parameters2_bao,parameters2_rsd, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin,Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao,k0rsd, P0bao,P0rsd,Pnoise,NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd,NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd,NeffP4bao,NeffP4rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC,PnoiseSGC,NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC,NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC,cov,covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1bao, plan2bao,plan1rsd, plan2rsd, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory_bao,spacing_theory_rsd);
//printf("chi2-FSbaoiso =%lf, %d\n",L_proposed,abs_counter_private);
}

if(strcmp(type_of_analysis, "FSBAOANISO") == 0){

parameters2[0]=input_vector_trial[0];//apara
parameters2[1]=input_vector_trial[1];//aperp
parameters2[2]=input_vector_trial[2];//f

parameters2_bao[0]=input_vector_trial[0];//apara
parameters2_bao[1]=input_vector_trial[1];//aperp

i1=2;i2=3;
if(Nsigmas_free==0)
{

  if(strcmp(Sigma_independent, "yes") == 0)//Nsigma_tot=2
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2_bao[3]=Sigma_nl_mean[1];//sigmaperp
  parameters2[3]=Sigma_nl_mean[0];//sigmapara
  parameters2[4]=Sigma_nl_mean[1];//sigmaperp

  i1=4;i2=3;

  }
  if(strcmp(Sigma_independent, "no") == 0)//Nsigma_tot=1
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2[3]=Sigma_nl_mean[0];//sigmapara

  i1=3;i2=3;
  }

}

for(i=i1+1;i<dimensionbao+1;i++){parameters2[i]=input_vector_trial[i2];parameters2_bao[i-1]=input_vector_trial[i2];i2++;}//set the rest of parameters (check we need dimension+1)

parameters2_rsd[0]=input_vector_trial[0];
parameters2_rsd[1]=input_vector_trial[1];
parameters2_rsd[2]=input_vector_trial[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=input_vector_trial[i2];parameters2[i-3+dimensionbao+1]=input_vector_trial[i2];i2++;}

//for(i=0;i<dimensionrsd;i++){printf("%d %lf\n",i,parameters2_rsd[i]);}

L_proposed=chi2_bao_rsd(type_BAO_fit,type_of_analysis,fit_BAO,fit_RSD,parameters2_bao,parameters2_rsd, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin,Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao,k0rsd, P0bao,P0rsd,Pnoise,NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd,NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd,NeffP4bao,NeffP4rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC,PnoiseSGC,NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC,NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC,cov,covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1bao, plan2bao, plan1rsd, plan2rsd, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory_bao,spacing_theory_rsd);
//printf("chi2-FSbaoaniso =%lf, %d\n",L_proposed,abs_counter_private);

}



acceptance=exp(-0.5*(L_proposed-L_current));
//random_acceptance=drand48();
random_acceptance=erand48(randstate);

//printf("%lf > %lf => accepted\n",acceptance,random_acceptance);

for(i=0;i<N_Cov_prop;i++)
{
if(input_vector_trial[i]>priors_high[i]){acceptance=0;}
if(input_vector_trial[i]<priors_low[i]){acceptance=0;}
//if(acceptance==0){printf("out by prior\n");}
}


//if(trial_mcmc==1){
//printf("Lcurrent=%lf, Ltrial=%lf eff counter: %ld (%d) accpetance=%lf, random=%lf\n",L_current,L_proposed,abs_counter_private,weight,acceptance,random_acceptance);
//for(i=0;i<N_Cov_prop;i++){printf("%d, %lf\n",i,input_vector_trial[i]);}
//}
//exit(0);


if(acceptance>random_acceptance)//accepted
{
//print old point with its weight
for(i=0;i<N_Cov_prop;i++){input_vector_accepted[i]=input_vector_trial[i];}

//exit(0);
//BAO only
if(strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){
//previous step is written with its weight

prior=0;
for(i=0;i<Nsigmas_tot;i++)
{
if(Sigma_type[i]==1){prior=prior+gauss(input_vector_trial[Nalphas+i],Sigma_nl_mean[i],Sigma_nl_stddev[i]);/*printf("prior2 %lf (%lf,%lf,%lf)\n",gauss(input_vector_trial[Nalphas+i],Sigma_nl_mean[i],Sigma_nl_stddev[i]),input_vector_trial[Nalphas+i],Sigma_nl_mean[i],Sigma_nl_stddev[i]);*/}
}

}

if(strcmp(type_of_analysis, "FS") == 0){
//nothing for RSD, 
prior=0;
}

if(trial_mcmc==0)
{
vector_print[0][counter]=weight;
vector_print[1][counter]=L_current-prior;
}
else
{

if(abs_counter_private>Nburnout){
vector_buffer[0][i_thread*N_max/nthreads + abs_counter_private - Nburnout-1]=weight;
vector_buffer[1][i_thread*N_max/nthreads + abs_counter_private - Nburnout-1]=L_current-prior;
}

}
for(i=2;i<N_Cov_prop+2;i++)
{
if(trial_mcmc==0){
vector_print[i][counter]=Doriginal[i-2];//print the old one with the corresponding weight

}
else{

if(abs_counter_private>Nburnout){
vector_buffer[i][i_thread*N_max/nthreads + abs_counter_private - Nburnout-1]=Doriginal[i-2];
}

}

}

if(abs_counter_private>Nburnout){//count only once the minimum number of steps has been reached
counter++;
}

abs_counter_private++;

if(counter==N_print && trial_mcmc==0)
{

//volcado en archivo
if (nthreads<2){
  f=fopen(name_file_output_mcmc,"a");
}
else{
  strcpy(name_wo_extension, name_file_output_mcmc);
  strtok(name_wo_extension, ".");
  sprintf(name_file_thread,"%s__%d.txt",name_wo_extension,i_thread+1);
  f=fopen(name_file_thread,"a");
}
for(l=0;l<counter;l++)
{
   weight_print=(int)(vector_print[0][l]);
   fprintf(f,"%d \t %e \t",weight_print,vector_print[1][l]);

   for(i=0;i<N_Cov_prop;i++)
   {
      if(i==N_Cov_prop-1){fprintf(f,"%e\n",vector_print[i+2][l]);}
      else{fprintf(f,"%e\t",vector_print[i+2][l]);}
   }

}

fclose(f);

freeTokens(vector_print,N_Cov_prop+2);

vector_print = (double**) calloc(N_Cov_prop+2,sizeof(double*));

for(j=0;j<N_Cov_prop+2;j++)
{
   vector_print[j] = (double*) calloc(N_print,sizeof(double));
}

convergence_private=0;

if(nthreads<2)
{
  //cada x iteraciones comprobar convergence y decidir si se sale. Numero maximo salir, numero minimo no salir nunca
  if(trial_mcmc==0 && abs_counter_private>Nburnout)//check convergence once the file is written
  {
  lines=abs_counter_private-(Nburnout+1)*nthreads;
  convergence_private=get_convergence(name_file_output_mcmc,N_Cov_prop,lines);//if convergence is 1 out
  }
}
else 
{
  if(trial_mcmc==0 && abs_counter_private>Nburnout)//check convergence once all files are written
  {
  convergence_private=get_convergence_parallel(nthreads,name_wo_extension,N_Cov_prop);//if convergence is 1 out
  }
}  
counter=0;
}

//update L_current and weight before leaving
L_current=L_proposed;
weight=1;

}
else
{//up-weight old
weight++;
if(weight>=10000){printf("Warning: Weight values > 50000. abs_count=%ld, chi2=%lf. Exiting now...\n",abs_counter_private,L_current);exit(0);} 
//if(weight>=N_print){printf("Warning: Weight values > %d (N_print). abs_count=%ld, chi2=%lf. Exiting now...\n",weight,abs_counter_private,L_current);exit(0);}
}


if(convergence_private==1){convergence_shared=1;break;}//out of do-while loop once convergence has been reached

}while((abs_counter_private-Nburnout-1<N_max/nthreads) && (convergence_shared==0));

abs_counter_shared += abs_counter_private;
j_run_shared += j_run_private;

free(Dwhite);
free(Doriginal);
free(step);
free(input_vector_accepted);
free(input_vector_trial);
free(parameters2);
if(strcmp(type_of_analysis_BAO,"yes")==0 && strcmp(type_of_analysis_FS,"yes")==0)
{
free(parameters2_bao);
free(parameters2_rsd);
}

if(counter!=0 && trial_mcmc==0)//free if hasn't been freed before (in the case of exit by abs_counter_shared>=Nsteps only)
{
freeTokens(vector_print,N_Cov_prop+2);
}

}

if(convergence_shared==0){

if(trial_mcmc==0)printf("Chain exited by number of steps: %ld > %ld + %d * %ld\n",abs_counter_shared,N_max,nthreads,Nburnout);

}
else{printf("Chain exited by convergence.\n");}

if(trial_mcmc==0){printf("%lf acceptance rate\n",(abs_counter_shared-(Nburnout+1)*nthreads)*1./(j_run_shared-(Nburnout+1)*nthreads)*1.);}

time_final=time(NULL);
//write log file and get minimum

if(trial_mcmc==0){

time_run=(time_final-time_ini)/(60.*60.);//in hours
lines=abs_counter_shared-(Nburnout+1)*nthreads;//check

if (nthreads>1){
concatenate_files(nthreads,name_wo_extension);
}

parameters2 =  (double*) calloc( dimension, sizeof(double));
if(strcmp(type_of_analysis_BAO,"yes")==0 && strcmp(type_of_analysis_FS,"yes")==0)
{
parameters2_bao =  (double*) calloc( dimensionbao, sizeof(double));
parameters2_rsd =  (double*) calloc( dimensionrsd, sizeof(double));

}
mean_vector=(double*) calloc(N_Cov_prop, sizeof(double));
min_vector=(double*) calloc(N_Cov_prop+1, sizeof(double));
do_log_file(nthreads,path_output,name_file_output_mcmc,N_Cov_prop,lines,time_run,abs_counter_shared,j_run_shared,Nburnout,identifier,mean_vector,min_vector);
write_prop_cov(name_file_output_mcmc,N_Cov_prop);

//for(i=0;i<N_Cov_prop+1;i++){printf("%d %lf\n",i,min_vector[i]);}
//exit(0);
    if( strcmp(do_plot, "yes") == 0 )
    {

Ndof=N_Cov_prop;
//parameters_plot =  (double*) calloc( dimension, sizeof(double));

/*
 for(j=0;j<dimension;j++)
 {

//For BAO-only
if(strcmp(type_of_analysis, "BAOISO") == 0){

     if( strcmp(Sigma_nl_type, "fixed") == 0 )
     {
          if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P4") == 0 )
          {
             if(j<1){parameters_plot[j]=min_vector[j];}
             if(j==1){parameters_plot[j]=Sigma_nl_mean;}
             if(j>1){parameters_plot[j]=min_vector[j-1];}
          }
          else
          {
             if(j<2){parameters_plot[j]=min_vector[j];}
             if(j==2){parameters_plot[j]=Sigma_nl_mean;}
             if(j>2){parameters_plot[j]=min_vector[j-1];}
          }
    
     }
     else
     {               
       parameters_plot[j]=min_vector[j];
     }
}

if(strcmp(type_of_analysis, "BAOANISO") == 0){


}
if(strcmp(type_of_analysis, "FS") == 0){  
//For RSD-here 
}

}
*/


//for BAO
if(strcmp(type_of_analysis, "BAOISO") == 0){
parameters2[0]=min_vector[0];i1=1;i2=1;
if(modeP0bao+modeP2bao+modeP4bao>1){parameters2[1]=min_vector[1];i1=2;i2=2;}
if(Nsigmas_free==0)
{
if(modeP0bao+modeP2bao+modeP4bao==1)
{
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
  parameters2[1]=Sigma_nl_mean[0];
  parameters2[2]=Sigma_nl_mean[1];i1=3;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0 ){parameters2[1]=Sigma_nl_mean[0];i1=2;}
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
    if(modeP0bao==1){parameters2[1]=Sigma_nl_mean[0];}
    if(modeP2bao==1){parameters2[1]=Sigma_nl_mean[1];}
    if(modeP4bao==1){parameters2[1]=Sigma_nl_mean[2];}
    i1=2;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
   parameters2[1]=Sigma_nl_mean[0];i1=2;
  }
}
if(modeP0bao+modeP2bao+modeP4bao>1){
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=Sigma_nl_mean[1];i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];parameters2[4]=Sigma_nl_mean[2];i1=5;}
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==0){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[1];i1=4;}
    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[2];i1=4;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==1){parameters2[2]=Sigma_nl_mean[1];parameters2[3]=Sigma_nl_mean[2];i1=4;}
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0){parameters2[2]=Sigma_nl_mean[0];i1=3;}
}
}
//for(i=Nalphas+Nsigmas_tot-Nsigmas_free;i<N_Cov_prop;i++){parameters2[i+Nsigmas_tot-Nsigmas_free]=min_vector[i];}//set the rest of parameters

//for(i=i1;i<N_Cov_prop+1;i++){parameters2[i]=min_vector[i2];i2++;}//set the rest of parameters BUG HERE
for(i=i1;i<dimension;i++){parameters2[i]=min_vector[i2];i2++;}//set the rest of parameters


//printf("Ncov+1=%d, dimension=%d\n",N_Cov_prop+1,dimension);

//for(i=0;i<20;i++){printf("\n%d %lf\n",i,parameters2[i]);}

             make_a_bao_plot(type_BAO_fit,type_of_analysis,fit_BAO,parameters2,min_vector[N_Cov_prop], k0bao,P0bao,errP0bao,NeffP0bao,k0baoSGC,P0baoSGC,errP0baoSGC,NeffP0baoSGC, k2bao,P2bao,errP2bao,NeffP2bao,k2baoSGC,P2baoSGC,errP2baoSGC,NeffP2baoSGC,k4bao,P4bao,errP4bao,NeffP4bao,k4baoSGC,P4baoSGC,errP4baoSGC,NeffP4baoSGC, k_Plin, Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev,  Npolynomial, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier,plan1bao,plan2bao,Nalphas,Nsigmas_tot, Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);
}

if(strcmp(type_of_analysis, "BAOANISO") == 0){
parameters2[0]=min_vector[0];
parameters2[1]=min_vector[1];
i1=2;i2=2;
if(Nsigmas_free==0)
{
  if(strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=Sigma_nl_mean[1];
  parameters2[4]=min_vector[2];
  i1=5;i2=3;

  }
  if(strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[2]=Sigma_nl_mean[0];
  parameters2[3]=min_vector[2];
  i1=4;i2=3;

  }

}

//for(i=i1;i<N_Cov_prop+1;i++){parameters2[i]=min_vector[i2];i2++;}//set the rest of parameters BUG HERE???
for(i=i1;i<dimension;i++){parameters2[i]=min_vector[i2];i2++;}//set the rest of parameters 

             make_a_bao_plot(type_BAO_fit,type_of_analysis,fit_BAO,parameters2,min_vector[N_Cov_prop], k0bao,P0bao,errP0bao,NeffP0bao,k0baoSGC,P0baoSGC,errP0baoSGC,NeffP0baoSGC, k2bao,P2bao,errP2bao,NeffP2bao,k2baoSGC,P2baoSGC,errP2baoSGC,NeffP2baoSGC,k4bao,P4bao,errP4bao,NeffP4bao,k4baoSGC,P4baoSGC,errP4baoSGC,NeffP4baoSGC, k_Plin, Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev,  Npolynomial, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier,plan1bao,plan2bao,Nalphas,Nsigmas_tot, Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);
}





//For RSD
if(strcmp(type_of_analysis, "FS") == 0){


//make_a_rsd_plot(type_of_analysis,min_vector,min_vector[N_Cov_prop], Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0rsd, P0rsd, errP0rsd, Pnoise, k2rsd, P2rsd,errP2rsd, k4rsd, P4rsd,errP4rsd, k0rsdSGC, P0rsdSGC,errP0rsdSGC, PnoiseSGC, k2rsdSGC, P2rsdSGC,errP2rsdSGC, k4rsdSGC, P4rsdSGC,errP4rsdSGC, NeffP0rsd, NeffP2rsd, NeffP4rsd,NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier, plan1, plan2, fit_RSD, do_power_spectrum, do_bispectrum,factor_sampling_mask,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd);

i2=0;
for(i=0;i<dimension;i++){parameters2[i]=min_vector[i2];i2++;}//set the rest of parameters 
make_a_rsd_plot(type_of_analysis,parameters2,min_vector[N_Cov_prop], Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0rsd, P0rsd, errP0rsd, Pnoise, k2rsd, P2rsd,errP2rsd, k4rsd, P4rsd,errP4rsd, k0rsdSGC, P0rsdSGC,errP0rsdSGC, PnoiseSGC, k2rsdSGC, P2rsdSGC,errP2rsdSGC, k4rsdSGC, P4rsdSGC,errP4rsdSGC, NeffP0rsd, NeffP2rsd, NeffP4rsd,NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier, plan1rsd, plan2rsd, fit_RSD, do_power_spectrum, do_bispectrum,factor_sampling_mask,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd);


}

//for BAO+RSD
if(strcmp(type_of_analysis, "FSBAOISO") == 0){

parameters2[0]=min_vector[0];//apara
parameters2[1]=min_vector[1];//aperp
parameters2[2]=min_vector[2];//f

//parameters2_bao[0]=min_vector[0];//apara
//parameters2_bao[1]=min_vector[1];//aperp

if(modeP0bao+modeP2bao+modeP4bao>1){
parameters2_bao[0]=min_vector[0];//apara
parameters2_bao[1]=min_vector[1];//aperp
baoiso_shift=1;
}
else
{
baoiso_shift=2;
if(modeP0bao==1){parameters2_bao[0]=pow(min_vector[0],1./3.)*pow(min_vector[1],2./3.);  }
if(modeP2bao==1){parameters2_bao[0]=pow(min_vector[0],3./5.)*pow(min_vector[1],2./5.);  }
if(modeP4bao==1){parameters2_bao[0]=pow(min_vector[0],5./7.)*pow(min_vector[1],2./7.);  }
}


i1=2;i2=3;
if(Nsigmas_free==0)
{

if(modeP0bao+modeP2bao+modeP4bao==1)
{

  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
  parameters2_bao[1]=Sigma_nl_mean[0];
  parameters2_bao[2]=Sigma_nl_mean[1];

  parameters2[3]=Sigma_nl_mean[0];
  parameters2[4]=Sigma_nl_mean[1];

  i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
     parameters2[3]=Sigma_nl_mean[0];
     parameters2_bao[1]=Sigma_nl_mean[0];
     i1=3;

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0 )
  {
    if(modeP0bao==1){parameters2_bao[1]=Sigma_nl_mean[0];parameters2[3]=Sigma_nl_mean[0];}
    if(modeP2bao==1){parameters2_bao[1]=Sigma_nl_mean[1];parameters2[3]=Sigma_nl_mean[1];}
    if(modeP4bao==1){parameters2_bao[1]=Sigma_nl_mean[2];parameters2[3]=Sigma_nl_mean[2];}
    i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0 )
  {
   parameters2[3]=Sigma_nl_mean[0];
   parameters2_bao[1]=Sigma_nl_mean[0];
   i1=3;
  }


}
if(modeP0bao+modeP2bao+modeP4bao>1){


  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
  parameters2_bao[2]=Sigma_nl_mean[0];
  parameters2_bao[3]=Sigma_nl_mean[1];
  parameters2[3]=Sigma_nl_mean[0];
  parameters2[4]=Sigma_nl_mean[1];

  i1=4;
  }
  if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "no") == 0)
  {
  parameters2[3]=Sigma_nl_mean[0];
  parameters2_bao[2]=Sigma_nl_mean[0];

  i1=3;
  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0)
  {
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[1];
         parameters2_bao[4]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[1];
         parameters2[5]=Sigma_nl_mean[2];
         i1=5;}
    if(modeP0bao==1 && modeP2bao==1 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[1];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[1];
         i1=4;}

    if(modeP0bao==1 && modeP2bao==0 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2_bao[3]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[0];
         parameters2[4]=Sigma_nl_mean[2];
         i1=4;}
   if(modeP0bao==1 && modeP2bao==0 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[0];
         parameters2[3]=Sigma_nl_mean[0];
         i1=3;}
    if(modeP0bao==0 && modeP2bao==1 && modeP4bao==0){
         parameters2_bao[2]=Sigma_nl_mean[1];
         parameters2[3]=Sigma_nl_mean[1];
         i1=3;}
    if(modeP0bao==0 && modeP2bao==0 && modeP4bao==1){
         parameters2_bao[2]=Sigma_nl_mean[2];
         parameters2[3]=Sigma_nl_mean[2];
         i1=3;}

  }
  if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "no") == 0){
          parameters2[3]=Sigma_nl_mean[0];parameters2_bao[2]=Sigma_nl_mean[0];
          i1=3;}
}
}
for(i=i1+1;i<dimensionbao+baoiso_shift;i++){parameters2[i]=min_vector[i2];parameters2_bao[i-baoiso_shift]=min_vector[i2];i2++;}//set the rest of parameters (check we need dimension+1)

parameters2_rsd[0]=min_vector[0];
parameters2_rsd[1]=min_vector[1];
parameters2_rsd[2]=min_vector[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=min_vector[i2];parameters2[i-3+dimensionbao+1]=min_vector[i2];i2++;}

//for(i=0;i<dimensionbao;i++){printf("BAO: %d %lf\n",i,parameters2_bao[i]);}
//for(i=0;i<dimensionrsd;i++){printf("RSD: %d %lf\n",i,parameters2_rsd[i]);}

    make_a_bao_plot(type_BAO_fit,type_of_analysis,fit_BAO,parameters2_bao,min_vector[N_Cov_prop], k0bao,P0bao,errP0bao,NeffP0bao,k0baoSGC,P0baoSGC,errP0baoSGC,NeffP0baoSGC, k2bao,P2bao,errP2bao,NeffP2bao,k2baoSGC,P2baoSGC,errP2baoSGC,NeffP2baoSGC,k4bao,P4bao,errP4bao,NeffP4bao,k4baoSGC,P4baoSGC,errP4baoSGC,NeffP4baoSGC, k_Plin, Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev,  Npolynomial, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier,plan1bao,plan2bao,Nalphas,Nsigmas_tot, Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);

make_a_rsd_plot(type_of_analysis,parameters2_rsd,min_vector[N_Cov_prop], Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0rsd, P0rsd, errP0rsd, Pnoise, k2rsd, P2rsd,errP2rsd, k4rsd, P4rsd,errP4rsd, k0rsdSGC, P0rsdSGC,errP0rsdSGC, PnoiseSGC, k2rsdSGC, P2rsdSGC,errP2rsdSGC, k4rsdSGC, P4rsdSGC,errP4rsdSGC, NeffP0rsd, NeffP2rsd, NeffP4rsd,NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier, plan1rsd, plan2rsd, fit_RSD, do_power_spectrum, do_bispectrum,factor_sampling_mask,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd);
}

if(strcmp(type_of_analysis, "FSBAOANISO") == 0){

parameters2[0]=min_vector[0];//apara
parameters2[1]=min_vector[1];//aperp
parameters2[2]=min_vector[2];//f

parameters2_bao[0]=min_vector[0];//apara
parameters2_bao[1]=min_vector[1];//aperp

i1=2;i2=3;
if(Nsigmas_free==0)
{

  if(strcmp(Sigma_independent, "yes") == 0)//Nsigma_tot=2
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2_bao[3]=Sigma_nl_mean[1];//sigmaperp
  parameters2[3]=Sigma_nl_mean[0];//sigmapara
  parameters2[4]=Sigma_nl_mean[1];//sigmaperp

  i1=4;i2=3;

  }
  if(strcmp(Sigma_independent, "no") == 0)//Nsigma_tot=1
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2[3]=Sigma_nl_mean[0];//sigmapara

  i1=3;i2=3;
  }

}

for(i=i1+1;i<dimensionbao+1;i++){parameters2[i]=min_vector[i2];parameters2_bao[i-1]=min_vector[i2];i2++;}//set the rest of parameters (check we need dimension+1)

parameters2_rsd[0]=min_vector[0];
parameters2_rsd[1]=min_vector[1];
parameters2_rsd[2]=min_vector[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=min_vector[i2];parameters2[i-3+dimensionbao+1]=min_vector[i2];i2++;}

//for(i=0;i<dimensionbao;i++){printf("BAO: %d %lf\n",i,parameters2_bao[i]);}
//for(i=0;i<dimensionrsd;i++){printf("RSD: %d %lf\n",i,parameters2_rsd[i]);}


 make_a_bao_plot(type_BAO_fit,type_of_analysis,fit_BAO,parameters2_bao,min_vector[N_Cov_prop], k0bao,P0bao,errP0bao,NeffP0bao,k0baoSGC,P0baoSGC,errP0baoSGC,NeffP0baoSGC, k2bao,P2bao,errP2bao,NeffP2bao,k2baoSGC,P2baoSGC,errP2baoSGC,NeffP2baoSGC,k4bao,P4bao,errP4bao,NeffP4bao,k4baoSGC,P4baoSGC,errP4baoSGC,NeffP4baoSGC, k_Plin, Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0,W2,W4,W6,W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev,  Npolynomial, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier,plan1bao,plan2bao,Nalphas,Nsigmas_tot, Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory_bao);

make_a_rsd_plot(type_of_analysis,parameters2_rsd,min_vector[N_Cov_prop], Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0rsd, P0rsd, errP0rsd, Pnoise, k2rsd, P2rsd,errP2rsd, k4rsd, P4rsd,errP4rsd, k0rsdSGC, P0rsdSGC,errP0rsdSGC, PnoiseSGC, k2rsdSGC, P2rsdSGC,errP2rsdSGC, k4rsdSGC, P4rsdSGC,errP4rsdSGC, NeffP0rsd, NeffP2rsd, NeffP4rsd,NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier, plan1rsd, plan2rsd, fit_RSD, do_power_spectrum, do_bispectrum,factor_sampling_mask,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd);


}

//free(parameters_plot);

    }

free(mean_vector);
free(min_vector);
free(parameters2);
if(strcmp(type_of_analysis_BAO,"yes")==0 && strcmp(type_of_analysis_FS,"yes")==0)
{
free(parameters2_bao);
free(parameters2_rsd);
}

}


/*
free(Dwhite);
free(Doriginal);
free(step);
free(input_vector_accepted);
free(input_vector_trial);
*/
gsl_matrix_free(transform_inverse);
gsl_matrix_free(transform);
if(strcmp(type_of_analysis_BAO,"yes")==0)
{
fftw_destroy_plan(plan1bao);
fftw_destroy_plan(plan2bao);
}
if(strcmp(type_of_analysis_FS,"yes")==0)
{
fftw_destroy_plan(plan1rsd);
fftw_destroy_plan(plan2rsd);
}
free(priors_low);
free(priors_high);
free(seed);
}

void do_bao_mcmc(int nthreads, char *type_BAO_fit,char *type_of_analysis,char *fit_BAO,double *k_Plin,double *Plin,int N_Plin, double *k_Olin, double *Olin, int N_Olin, double *pos, double *W0, double *W2, double *W4,double *W6, double *W8,int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC,  double *k0, double *P0, double *errP0, int NeffP0, double *k2, double *P2, double *errP2, int NeffP2, double *k4, double *P4, double *errP4, int NeffP4, double *k11, double *k22, double *k33, double *B0, double *errB0, double *Bnoise, int NeffB0, double *k0SGC, double *P0SGC, double *errP0SGC,int NeffP0SGC, double *k2SGC, double *P2SGC, double *errP2SGC,int NeffP2SGC, double *k4SGC, double *P4SGC, double *errP4SGC,int NeffP4SGC, double *k11SGC, double *k22SGC, double *k33SGC,double *B0SGC, double *errB0SGC, double *BnoiseSGC,int NeffB0SGC, double *cov, double *covSGC, double alpha_min, double alpha_max,  char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int Npolynomial, int Nchunks, char *path_output, char *identifier, char *do_plot, char *use_prop_cov, char *path_to_cov, long int Nsteps, char *do_power_spectrum, char *do_bispectrum, double Sigma_smooth,char *spacing_dataNGC,char *spacing_dataSGC, char *spacing_theory,char *type_of_analysis_BAO,char *type_of_analysis_FS)
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
int modeP0,modeP2,modeP4;
int allsigmafixed,Nalphas,Nsigmas_tot,Nsigmas_free;
int factor_sampling_mask;
int i_thread;

factor_sampling_mask=10;//Sampling factor boost wrt to Neff max when the mask is applied. How much do I need to sample my model pre-mask-apply in order to have a satisfactory mask response? 10 times seems reasonable.

if (nthreads<2){
 
  if(strcmp(type_of_analysis, "BAOISO") == 0){sprintf(name_file,"%s/mcmcBAOISO_output_%s.txt",path_output,identifier);}
  else if(strcmp(type_of_analysis, "BAOANISO") == 0){sprintf(name_file,"%s/mcmcBAOANISO_output_%s.txt",path_output,identifier);}
  else{sprintf(name_file,"%s/mcmcFS_output_%s.txt",path_output,identifier);}
  f=fopen(name_file,"w");
  if(f==NULL){printf("Error writing %s. Exiting now...\n",name_file);exit(0);}
  fclose(f);
} else {
  for (i_thread=0;i_thread<nthreads;i_thread++){
    
    if(strcmp(type_of_analysis, "BAOISO") == 0){sprintf(name_file,"%s/mcmcBAOISO_output_%s__%d.txt",path_output,identifier,i_thread+1);}
    else if(strcmp(type_of_analysis, "BAOANISO") == 0){sprintf(name_file,"%s/mcmcBAOANISO_output_%s__%d.txt",path_output,identifier,i_thread+1);}
    else{sprintf(name_file,"%s/mcmcFS_output_%s__%d.txt",path_output,identifier,i_thread+1);}
    f=fopen(name_file,"w");  
    if(f==NULL){printf("Error writing %s. Exiting now...\n",name_file);exit(0);}
    fclose(f);
  }
  
  if(strcmp(type_of_analysis, "BAOISO") == 0){sprintf(name_file,"%s/mcmcBAOISO_output_%s.txt",path_output,identifier);}
  else if(strcmp(type_of_analysis, "BAOANISO") == 0){sprintf(name_file,"%s/mcmcBAOANISO_output_%s.txt",path_output,identifier);}
  else{sprintf(name_file,"%s/mcmcFS_output_%s.txt",path_output,identifier);}
}

/*
//Number of free(varied) parameters in each case. 
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P4") == 0)
{
if(strcmp(Sigma_nl_type, "fixed")==0){N_Cov_prop=(Npolynomial+1)*Nchunks+1;}// sigma nl is not varied
else{N_Cov_prop=(Npolynomial+1)*Nchunks+2;}// sigma nl, alpha_para, alpha_perp are varied
}
if(strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0)
{
if(strcmp(Sigma_nl_type, "fixed")==0){N_Cov_prop=(Npolynomial+1)*Nchunks*2+2;}// sigma nl is not varied
else{N_Cov_prop=(Npolynomial+1)*Nchunks*2+3;}// sigma nl, alpha_para, alpha_perp are varied
}

if(strcmp(fit_BAO, "P024") == 0)
{
if(strcmp(Sigma_nl_type, "fixed")==0){N_Cov_prop=(Npolynomial+1)*Nchunks*3+2;}// sigma nl is not varied
else{N_Cov_prop=(Npolynomial+1)*Nchunks*3+3;}// sigma nl, alpha_para, alpha_perp are varied
}
*/
modeP0=0;
modeP2=0;
modeP4=0;
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

allsigmafixed=-1;
Nalphas=1;if(modeP0+modeP2+modeP4>1){Nalphas=2;}
Nsigmas_free=0;
if(strcmp(Sigma_independent, "yes") == 0 ){

    if(strcmp(Sigma_def_type, "para-perp") == 0)//This is the only possible case for BAOANISO
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
if(Sigma_type[0]>0){Nsigmas_free=1;}
if(Sigma_type[0]==0){Nsigmas_free=0;}
}
if(Nsigmas_free==0){allsigmafixed=0;}

//old baoaniso
//N_Cov_prop=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+Nalphas+Nsigmas_free;//number free  parameters
//if(strcmp(type_of_analysis, "BAOANISO") == 0){N_Cov_prop++;}//add beta as free parameter

if(strcmp(type_of_analysis, "BAOISO") == 0){N_Cov_prop=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+Nalphas+Nsigmas_free;}
if(strcmp(type_of_analysis, "BAOANISO") == 0){N_Cov_prop=Nchunks*(1+Npolynomial*(modeP0+modeP2+modeP4))+Nalphas+Nsigmas_free+1;}


//printf("%d\n",N_Cov_prop++);
//exit(0);
if(strcmp(Sigma_def_type, "para-perp") == 0 && strcmp(Sigma_independent, "yes") == 0){Nsigmas_tot=2;}
if(strcmp(Sigma_def_type, "effective") == 0 && strcmp(Sigma_independent, "yes") == 0){Nsigmas_tot=modeP0+modeP2+modeP4;}
if( strcmp(Sigma_independent, "no") == 0 ){Nsigmas_tot=1;}
//dimension=(Npolynomial+1)*Nchunks*(modeP0+modeP2+modeP4)+Nalphas+Nsigmas_tot;//number total parameters


//if(NeffB0>0){N_Cov_prop=N_Cov_prop+3;}//change this in the future


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

//mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,NULL, fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO, k_Plin, Plin, N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0, P0,errP0,NeffP0, k2, P2,errP2,NeffP2, k4, P4,errP4,NeffP4,k11,k22,k33,B0,errB0,Bnoise,NeffB0, k0SGC, P0SGC,errP0SGC,NeffP0SGC, k2SGC, P2SGC,errP2SGC,NeffP2SGC, k4SGC, P4SGC,errP4SGC,NeffP4SGC,k11SGC,k22SGC,k33SGC,B0SGC,BnoiseSGC,NeffB0SGC, cov, covSGC, alpha_min,alpha_max, Sigma_def_type, Sigma_independent,  ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks,  path_output, identifier, do_plot,Nsteps,do_power_spectrum,do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free, NULL, 0, 0, NULL, NULL,NULL,NULL,NULL, NULL, NULL,NULL,NULL,NULL,Sigma_smooth,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory);
mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,NULL, fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,NULL, k_Plin, Plin, N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0,NULL, P0,NULL,errP0,NULL,NeffP0,0, k2,NULL, P2,NULL,errP2,NULL,NeffP2,0, k4,NULL, P4,NULL,errP4,NULL,NeffP4,0,k11,NULL,k22,NULL,k33,NULL,B0,NULL,errB0,NULL,Bnoise,NULL,NeffB0,0, k0SGC,NULL, P0SGC,NULL,errP0SGC,NULL,NeffP0SGC,0, k2SGC,NULL, P2SGC,NULL,errP2SGC,NULL,NeffP2SGC,0, k4SGC,NULL, P4SGC,NULL,errP4SGC,NULL,NeffP4SGC,0,k11SGC,NULL,k22SGC,NULL,k33SGC,NULL,B0SGC,NULL,BnoiseSGC,NULL,NeffB0SGC,0, cov, covSGC, alpha_min,alpha_max, Sigma_def_type, Sigma_independent,  ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks,  path_output, identifier, do_plot,Nsteps,do_power_spectrum,do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free, NULL,0, 0, 0, NULL, NULL,NULL,NULL,NULL, NULL, NULL,NULL,NULL,NULL,Sigma_smooth,factor_sampling_mask,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,NULL,type_of_analysis_BAO,type_of_analysis_FS);


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

//mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer, fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO, k_Plin, Plin, N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0, P0,errP0,NeffP0, k2, P2,errP2,NeffP2, k4, P4,errP4,NeffP4,k11,k22,k33,B0,errB0,Bnoise,NeffB0, k0SGC, P0SGC,errP0SGC,NeffP0SGC, k2SGC, P2SGC,errP2SGC,NeffP2SGC, k4SGC, P4SGC,errP4SGC,NeffP4SGC,k11SGC,k22SGC,k33SGC,B0SGC,BnoiseSGC,NeffB0SGC, cov, covSGC, alpha_min,alpha_max, Sigma_def_type, Sigma_independent,  ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot,Nsteps,do_power_spectrum,do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free,NULL, 0, 0, NULL, NULL,NULL,NULL,NULL, NULL, NULL,NULL,NULL,NULL,Sigma_smooth,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory);
mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer, fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,NULL, k_Plin, Plin, N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0,NULL, P0,NULL,errP0,NULL,NeffP0,0, k2,NULL, P2,NULL,errP2,NULL,NeffP2,0, k4,NULL, P4,NULL,errP4,NULL,NeffP4,0,k11,NULL,k22,NULL,k33,NULL,B0,NULL,errB0,NULL,Bnoise,NULL,NeffB0,0, k0SGC,NULL, P0SGC,NULL,errP0SGC,NULL,NeffP0SGC,0, k2SGC,NULL, P2SGC,NULL,errP2SGC,NULL,NeffP2SGC,0, k4SGC,NULL, P4SGC,NULL,errP4SGC,NULL,NeffP4SGC,0,k11SGC,NULL,k22SGC,NULL,k33SGC,NULL,B0SGC,NULL,BnoiseSGC,NULL,NeffB0SGC,0, cov, covSGC, alpha_min,alpha_max, Sigma_def_type, Sigma_independent,  ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks,  path_output, identifier, do_plot,Nsteps,do_power_spectrum,do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free, NULL,0, 0, 0, NULL, NULL,NULL,NULL,NULL, NULL, NULL,NULL,NULL,NULL,Sigma_smooth,factor_sampling_mask,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,NULL,type_of_analysis_BAO,type_of_analysis_FS);

free(vector_mean);

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

//mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop,transform, transform_inverse,vector_mean, name_file, fit_BAO,k_Plin,Plin,N_Plin,k_Olin,Olin,N_Olin,pos,W0,W2,W4,W6,W8,Nmask,path_to_mask1,spacing_maskNGC,posSGC,W0SGC,W2SGC,W4SGC, W6SGC, W8SGC,NmaskSGC, path_to_mask2,spacing_maskSGC, k0, P0, errP0,NeffP0,k2,P2, errP2, NeffP2, k4,P4, errP4,NeffP4, k11,k22, k33, B0, errB0, Bnoise, NeffB0, k0SGC, P0SGC, errP0SGC,NeffP0SGC, k2SGC, P2SGC, errP2SGC,NeffP2SGC, k4SGC, P4SGC, errP4SGC,NeffP4SGC, k11SGC, k22SGC, k33SGC,B0SGC, BnoiseSGC,NeffB0SGC, cov,covSGC, alpha_min, alpha_max,  Sigma_def_type, Sigma_independent,  ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free, NULL, 0, 0, NULL, NULL,NULL,NULL,NULL, NULL, NULL,NULL,NULL,NULL, Sigma_smooth,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory);
mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer, fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,NULL, k_Plin, Plin, N_Plin,k_Olin,Olin,N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0,NULL, P0,NULL,errP0,NULL,NeffP0,0, k2,NULL, P2,NULL,errP2,NULL,NeffP2,0, k4,NULL, P4,NULL,errP4,NULL,NeffP4,0,k11,NULL,k22,NULL,k33,NULL,B0,NULL,errB0,NULL,Bnoise,NULL,NeffB0,0, k0SGC,NULL, P0SGC,NULL,errP0SGC,NULL,NeffP0SGC,0, k2SGC,NULL, P2SGC,NULL,errP2SGC,NULL,NeffP2SGC,0, k4SGC,NULL, P4SGC,NULL,errP4SGC,NULL,NeffP4SGC,0,k11SGC,NULL,k22SGC,NULL,k33SGC,NULL,B0SGC,NULL,BnoiseSGC,NULL,NeffB0SGC,0, cov, covSGC, alpha_min,alpha_max, Sigma_def_type, Sigma_independent,  ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks,  path_output, identifier, do_plot,Nsteps,do_power_spectrum,do_bispectrum,Nalphas,Nsigmas_tot,Nsigmas_free, NULL,0, 0, 0, NULL, NULL,NULL,NULL,NULL, NULL, NULL,NULL,NULL,NULL,Sigma_smooth,factor_sampling_mask,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,NULL,type_of_analysis_BAO,type_of_analysis_FS);

}


}

