/*
 * BRASS: Bao and Rsd Algorithm for Spectroscopic Surveys
 * All rights reserved
 * Author: Hector Gil Marin
 * Date: 5th Mar 2021
 * email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu
 * */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "bao.h"
#include "functions.h"
#include "read.h"
#include "mcmc_bao.h"
#include "rsd.h"

int main(int argc, char *argv[])
{
FILE *f,*g;
int i;

char name_file_ini[2000];//name of inizialization file

char name_file[2000];
char spacing_dataNGC_bao[2000],spacing_dataSGC_bao[2000];
char spacing_dataNGC_rsd[2000],spacing_dataSGC_rsd[2000];
char spacing_theory[2000],spacing_theory2[2000];
char spacing_theory_rsd[2000];
char spacing_maskNGC[2000], spacing_maskSGC[2000];
char type_of_analysis_FS[20],type_of_analysis_BAO[20];
sprintf(type_of_analysis_FS,"no");
sprintf(type_of_analysis_BAO,"no");

char type_of_analysis[10];
char fit_BAO[1000],fit_RSD[1000];
double kminP0bao,kmaxP0bao;
double kminP2bao,kmaxP2bao;
double kminP4bao,kmaxP4bao;
double kminB0bao,kmaxB0bao;

double kminP0rsd,kmaxP0rsd;
double kminP2rsd,kmaxP2rsd;
double kminP4rsd,kmaxP4rsd;
double kminB0rsd,kmaxB0rsd;
double kmin,kmax;

int Nrealizations,Nrealizations_used_NGC,Nrealizations_used_SGC;
char path_to_data1_bao[2000];
char path_to_data2_bao[2000];
char path_to_mocks1_bao[2000];
char path_to_mocks2_bao[2000];
char *path_to_mocks1_ending;
char *path_to_mocks2_ending;
char path_to_data1_rsd[2000];
char path_to_data2_rsd[2000];
char path_to_mocks1_rsd[2000];
char path_to_mocks2_rsd[2000];
char covfile1[2000];
char covfile2[2000];
char path_to_mask1[2000];
char path_to_mask2[2000];
char mask_renormalization[20];
char mask_matrix[20];
char covariance_correction[100];//hartlap, sellentin-heavens
char covarianceFSa_option[200];//cov. for FSalpha varying or fixed for each real
double **MatrixFS_mask_NGC,**MatrixFS_mask_SGC;
double **MatrixBAO_mask_NGC,**MatrixBAO_mask_SGC;
int N_matrixFS_in,N_matrixBAO_in;
int N_matrixFS_in_SGC,N_matrixBAO_in_SGC;
double sumw_bao,I22_bao,sumw_rsd,I22_rsd;
char perturbation_theory_file[2000];
char ptmodel_ps[2000];
char rsdmodel_ps[2000];
char fogmodel_ps[2000];
char ptmodel_bs[2000];
char local_b2s2[2000];
char local_b3nl[2000];
char RSD_fit[2000];
char sigma8_free[2000];
char fog_free[2000];
char fog_bs[2000];

int Npolynomial;
int Nchunks;
char Sigma_def_type[2000];//effective para-perp
char Sigma_independent[2000];//yes no
double Sigma_nl_mean[6],Sigma_nl_stddev[6];
double Sigma_type[6];
double FSprior_mean[5],FSprior_stddev[5];
double FSprior_type[5];
char Sigma_readin[2000];
int noise_option;
double ffactor;

int nthreads;
char type_BAO_fit[10];
char path_to_Plin[2000];
char path_to_Olin[2000];
char use_prop_cov[2000];
char path_to_cov[2000];
long int Nsteps;
double step_size;
char do_bispectrum[10];
char bispectrum_BQ[2];

char do_compression[100];
char path_to_compression[2000];
double **matrix_linear_compression; matrix_linear_compression=NULL;
int input_elements_data_vector,output_elements_data_vector,bao,rsd,check,modeP0bao,modeP2bao,modeP4bao,modeB0bao,modeP0rsd,modeP2rsd,modeP4rsd,modeB0rsd;

char do_power_spectrum[10];
char path_to_data1_bis_bao[2000];
char path_to_data2_bis_bao[2000];
char path_to_mocks1_bis_bao[2000];
char path_to_mocks2_bis_bao[2000];
char path_to_data1_bis_rsd[2000];
char path_to_data2_bis_rsd[2000];
char path_to_mocks1_bis_rsd[2000];
char path_to_mocks2_bis_rsd[2000];
char path_output[2000];
char do_plot[10];
char identifier[1000];
int NeffP0bao,NeffP2bao,NeffP4bao,NeffB0bao,Nmask,NmaskSGC;//number of elements withing specified k-ranges
int  NeffP0rsd,NeffP2rsd,NeffP4rsd,NeffB0rsd;
int NeffP0baoSGC,NeffP2baoSGC,NeffP4baoSGC,NeffB0baoSGC;
int NeffP0rsdSGC,NeffP2rsdSGC,NeffP4rsdSGC,NeffB0rsdSGC; 
double *k0bao,*P0bao,*k2bao,*P2bao,*k4bao,*P4bao,*errP0bao,*errP2bao,*errP4bao;
double *k0rsd,*P0rsd,*k2rsd,*P2rsd,*k4rsd,*P4rsd,*errP0rsd,*errP2rsd,*errP4rsd;
double *kav0bao,*kav2bao,*kav4bao;
double *kav0rsd,*kav2rsd,*kav4rsd;
double *k0baoSGC,*P0baoSGC,*k2baoSGC,*P2baoSGC,*k4baoSGC,*P4baoSGC,*errP0baoSGC,*errP2baoSGC,*errP4baoSGC;
double *k0rsdSGC,*P0rsdSGC,*k2rsdSGC,*P2rsdSGC,*k4rsdSGC,*P4rsdSGC,*errP0rsdSGC,*errP2rsdSGC,*errP4rsdSGC;
double *kav0baoSGC,*kav2baoSGC,*kav4baoSGC;
double *kav0rsdSGC,*kav2rsdSGC,*kav4rsdSGC;

double Pnoise_bao,Pnoise_baoSGC;
double Pnoise_rsd,Pnoise_rsdSGC;

double *cov,*covSGC;

int Ncov,NcovP;
int N_inputs;

double *posAV,*pos,*W0,*W2,*W4,*W6,*W8;
double *posAVSGC,*posSGC,*W0SGC,*W2SGC,*W4SGC,*W6SGC,*W8SGC;

double *k11bao,*k22bao,*k33bao,*B0bao,*Bnoise_bao,*errB0bao;
double *k11rsd,*k22rsd,*k33rsd,*B0rsd,*Bnoise_rsd,*errB0rsd;


double *k11baoSGC,*k22baoSGC,*k33baoSGC,*B0baoSGC,*Bnoise_baoSGC,*errB0baoSGC;
double *k11rsdSGC,*k22rsdSGC,*k33rsdSGC,*B0rsdSGC,*Bnoise_rsdSGC,*errB0rsdSGC;


double **Theory;
int Ntheory;
double *k_in_mask;

double params[10];
int params_int[2];
double params_ks[2];
double paramsBAO[3];

int  N_Plin,N_Olin;
double *k_Plin, *Plin, *k_Olin, *Olin;

double alpha_min,alpha_max;
double alpha_step;

double Sigma_smooth;
int fraction;
int bispectrum_needs=0;
int olin_exists=0;

int factor_for_sampling=5;//sampling factor for applying the mask.


NeffP0bao=0;
NeffP2bao=0;
NeffP4bao=0;
NeffB0bao=0;
NeffP0rsd=0;
NeffP2rsd=0;
NeffP4rsd=0;
NeffB0rsd=0;

Nmask=0;
NmaskSGC=0;
NeffP0baoSGC=0;
NeffP2baoSGC=0;
NeffP4baoSGC=0;
NeffB0baoSGC=0;
NeffP0rsdSGC=0;
NeffP2rsdSGC=0;
NeffP4rsdSGC=0;
NeffB0rsdSGC=0;


sprintf(name_file_ini,argv[1]);
f=fopen(name_file_ini,"r");
if(f==NULL){printf("File %s not found...\t exiting now\n",name_file_ini);return 0;}
else{printf("Reading Inizialization file: %s\n\n",name_file_ini);}

fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %s\n\n",type_of_analysis);
if(strcmp(type_of_analysis, "BAOANISO") != 0 && strcmp(type_of_analysis, "BAOISO") != 0 && strcmp(type_of_analysis, "FS") != 0 && strcmp(type_of_analysis, "FSBAOISO") != 0 && strcmp(type_of_analysis, "FSBAOANISO") != 0 && strcmp(type_of_analysis,"FSalphasrecon") !=0 ){printf("Error with type of computation: %s. Exiting now...\n",type_of_analysis);exit(0);}

printf("%s-type computation\n",type_of_analysis);

//Power Spectrum reading stuff
fscanf(f,"%*s %*s %*s %*s %s\n",do_power_spectrum);
if(strcmp(do_power_spectrum, "yes") != 0 && strcmp(do_power_spectrum, "no") != 0){printf("Error do power spectrum has to be either yes or no: %s\n. Exiting now...\n",do_power_spectrum);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",fit_BAO);
if(strcmp(fit_BAO, "P0") != 0 && strcmp(fit_BAO, "P2") != 0 && strcmp(fit_BAO, "P4") != 0  && strcmp(fit_BAO, "P02") != 0  && strcmp(fit_BAO, "P04") != 0  && strcmp(fit_BAO, "P24") != 0 && strcmp(fit_BAO, "P024") != 0){printf("Error with power spectrum multipoles to be fitted: %s. Exiting now...\n",fit_BAO);exit(0);}

if(strcmp(type_of_analysis,"FSalphasrecon") == 0){sprintf(fit_BAO,"P0");}//just P0 for FSalphasrecon, regardless of the input. 

fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&kminP0bao,&kmaxP0bao);
if(kminP0bao>=kmaxP0bao || kminP0bao<0 || kmaxP0bao<0){printf("Error with kmin or kmax values P0-BAO: %lf %lf. Exiting now...\n",kminP0bao,kmaxP0bao);exit(0);}

fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&kminP2bao,&kmaxP2bao);
if(kminP2bao>=kmaxP2bao || kminP2bao<0 || kmaxP2bao<0){printf("Error with kmin or kmax values P2-BAO: %lf %lf. Exiting now...\n",kminP2bao,kmaxP2bao);exit(0);}

fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&kminP4bao,&kmaxP4bao);
if(kminP4bao>=kmaxP4bao || kminP4bao<0 || kmaxP4bao<0){printf("Error with kmin or kmax values P4-BAO: %lf %lf. Exiting now...\n",kminP4bao,kmaxP4bao);exit(0);}

fscanf(f,"%*s %*s %*s %*s %s\n",fit_RSD);
if(strcmp(fit_RSD, "P0") != 0 && strcmp(fit_RSD, "P2") != 0 && strcmp(fit_RSD, "P4") != 0  && strcmp(fit_RSD, "P02") != 0  && strcmp(fit_RSD, "P04") != 0  && strcmp(fit_RSD, "P24") != 0 && strcmp(fit_RSD, "P024") != 0){printf("Error with power spectrum multipoles to be fitted: %s. Exiting now...\n",fit_RSD);exit(0);}

fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&kminP0rsd,&kmaxP0rsd);
if(kminP0rsd>=kmaxP0rsd || kminP0rsd<0 || kmaxP0rsd<0){printf("Error with kmin or kmax values P0-RSD: %lf %lf. Exiting now...\n",kminP0rsd,kmaxP0rsd);exit(0);}

fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&kminP2rsd,&kmaxP2rsd);
if(kminP2rsd>=kmaxP2rsd || kminP2rsd<0 || kmaxP2rsd<0){printf("Error with kmin or kmax values P2-RSD: %lf %lf. Exiting now...\n",kminP2rsd,kmaxP2rsd);exit(0);}

fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&kminP4rsd,&kmaxP4rsd);
if(kminP4rsd>=kmaxP4rsd || kminP4rsd<0 || kmaxP4rsd<0){printf("Error with kmin or kmax values P4-RSD: %lf %lf. Exiting now...\n",kminP4rsd,kmaxP4rsd);exit(0);}

fscanf(f,"%*s %*s %*s %d\n\n",&Nchunks);
if(Nchunks!=1 && Nchunks!=2){printf("Error with the number of chunks %d\n. Exiting now...\n",Nchunks);exit(0);}

fscanf(f,"%*s %*s %*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",path_to_data1_bao);
fscanf(f,"%*s %*s %*s %s\n",path_to_data2_bao);
fscanf(f,"%*s %*s %*s %s\n",path_to_data1_rsd);
fscanf(f,"%*s %*s %*s %s\n",path_to_data2_rsd);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks1_bao);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks2_bao);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks1_rsd);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks2_rsd);
fscanf(f,"%*s %*s %*s %d\n",&Nrealizations);
if(Nrealizations<0){printf("Error with number of realizations %d\n. Exiting now...\n",Nrealizations);exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %s\n",path_to_mask1);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %s\n",path_to_mask2);
fscanf(f,"%*s %*s %*s %s\n",mask_renormalization);
if(strcmp(mask_renormalization, "no") != 0 && strcmp(mask_renormalization, "yes") != 0){printf("Error with answer to Window Renormalization %s\n. Exiting now...\n",mask_renormalization);exit(0);}
//#Apply mask as a matrix (yes/no): no
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n",mask_matrix);
if(strcmp(mask_matrix, "no") != 0 && strcmp(mask_matrix, "yes") != 0){printf("Error with answer to Window Mask Matrix %s\n. Exiting now...\n",mask_matrix);exit(0);}
if(strcmp(mask_matrix, "yes") == 0 && strcmp(path_to_mask1,"none") == 0){printf("Error, Mask Matrix set to yes, but no RR counts provided. Exiting now ...\n");exit(0);}
if(Nchunks==2){if(strcmp(mask_matrix, "yes") == 0 && strcmp(path_to_mask2,"none") == 0){printf("Error, Mask Matrix set to yes, but no RR counts provided. Exiting now ...\n");exit(0);}}
fscanf(f,"%*s %*s %*s %s\n",covariance_correction);//covariance_correction
if(strcmp(covariance_correction, "none") != 0 && strcmp(covariance_correction, "Hartlap") != 0 && strcmp(covariance_correction, "Sellentin-Heavens") != 0){printf("Error with answer to covariance correction: %s\n. Exiting now...\n",covariance_correction);exit(0);}
if( Nrealizations==0 && strcmp(covariance_correction, "none") != 0){printf("Error, can't do Hartlap or Sellentin-Heavens correction with 0 realizations. Exiting now...\n");exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n\n",covarianceFSa_option);//covariance_correction
if( strcmp(type_of_analysis,"FSalphasrecon") ==0 )//only relevant for FSalphasrecon
{
if( strcmp(covarianceFSa_option,"fixed")!=0 && strcmp(covarianceFSa_option,"varying")!=0){printf("Error with answer to covariance option for FSalpha correction: %s\n. Exiting now...\n",covarianceFSa_option);exit(0);}
}
else
{
sprintf(covarianceFSa_option,"fixed");//if analysis is not FSalphasrecon != 0 then fixed
}

fscanf(f,"%*s %*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",do_bispectrum);
if(strcmp(do_bispectrum, "no") != 0 && strcmp(do_bispectrum, "yes") != 0){printf("Error with answer to using bispectrum %s\n. Exiting now...\n",do_bispectrum);exit(0);}
if(strcmp(do_bispectrum, "no") == 0 && strcmp(do_power_spectrum, "no") == 0){printf("Error. You must choose either power spectrum or bispectrum as inputs. Exiting now...\n");exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %s\n",bispectrum_BQ);
if( strcmp(do_bispectrum, "yes") != 0){
      if( strcmp(bispectrum_BQ,"B") !=0 && strcmp(bispectrum_BQ,"Q") != 0 ){  printf("Error %s, you must choose between bispectrum (B) or reduced bispectrum (Q). Exiting now...\n",bispectrum_BQ);exit(0);}
}

fscanf(f,"%*s %*s %*s %s\n",path_to_data1_bis_bao);
fscanf(f,"%*s %*s %*s %s\n",path_to_data2_bis_bao);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks1_bis_bao);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks2_bis_bao);
fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n\n",&kminB0bao,&kmaxB0bao);
if(kminB0bao>=kmaxB0bao || kminB0bao<0 || kmaxB0bao<0){printf("Error with kmin or kmax values B0: %lf %lf. Exiting now...\n",kminB0bao,kmaxB0bao);exit(0);}

fscanf(f,"%*s %*s %*s %s\n",path_to_data1_bis_rsd);
fscanf(f,"%*s %*s %*s %s\n",path_to_data2_bis_rsd);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks1_bis_rsd);
fscanf(f,"%*s %*s %*s %s\n",path_to_mocks2_bis_rsd);
fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n\n",&kminB0rsd,&kmaxB0rsd);
if(kminB0rsd>=kmaxB0rsd || kminB0rsd<0 || kmaxB0rsd<0){printf("Error with kmin or kmax values B0: %lf %lf. Exiting now...\n",kminB0rsd,kmaxB0rsd);exit(0);}


fscanf(f,"%*s %*s %*s %s\n",do_compression);
fscanf(f,"%*s %*s %*s %*s %s\n\n",path_to_compression);

if( strcmp(do_compression, "no") != 0 && strcmp(do_compression, "linear") != 0 && strcmp(do_compression, "nonlinear") != 0){printf("Error, options for Do Compression should be no/linear/nonlinear (%s). Exiting now...\n",do_compression);exit(0);}
if( strcmp(do_compression, "linear") == 0){
g=fopen(path_to_compression,"r");
if(g==NULL){printf("Error reading the compression file for linear compression %s. Exiting now...\n",path_to_compression);exit(0);}
else{fclose(g);}
}

if( strcmp(do_compression, "nonlinear") == 0){printf("Nonlinear compression not available at the moment!. Exiting now...\n");exit(0);}

///
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %d\n",&Npolynomial);
if(Npolynomial<0 || Npolynomial>10){printf("Error with the polynomial order %d\n. Exiting now...\n",Npolynomial);exit(0);}
fscanf(f,"%*s %*s %*s %s\n",Sigma_def_type);
if(strcmp(Sigma_def_type, "effective") != 0 && strcmp(Sigma_def_type, "para-perp") != 0){printf("Error Sigma_def type %s\n. Exiting now...\n",Sigma_def_type);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",Sigma_independent);
if(strcmp(Sigma_independent, "yes") != 0 && strcmp(Sigma_independent, "no") != 0){printf("Error Sigma_independent type %s\n. Exiting now...\n",Sigma_independent);exit(0);}
fscanf(f,"%*s %*s %*s %*s %lf\n\n",&ffactor);
if(ffactor<0){printf("Warning f-factor negative: %lf. Exiting  now...\n",ffactor);exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %lf\n",&Sigma_smooth);
if(Sigma_smooth<0){printf("Warning unpysical value for smoothing scale: %lf. Exiting  now...\n",Sigma_smooth);exit(0);}
if(Sigma_smooth==0){printf("Smoothing scale set to %lf (pre-recon analysis).\n",Sigma_smooth);}
else{printf("Smoothing scale set to %lf Mpc/h (post-recon analysis)\n",Sigma_smooth);}

fscanf(f,"%*s\n");
//Sigma0
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "fixed") != 0 && strcmp(Sigma_readin, "free") != 0){printf("Error Sigma_nl type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&Sigma_nl_mean[0],&Sigma_nl_stddev[0]);
if(Sigma_nl_mean[0]<0 || Sigma_nl_stddev[0]<0){printf("Error with Sigma_nl values: %lf %lf. Exiting now...\n",Sigma_nl_mean[0],Sigma_nl_stddev[0]);exit(0);}
if(strcmp(Sigma_readin, "fixed") == 0){Sigma_type[0]=0;}
if(strcmp(Sigma_readin, "gaussian") == 0){Sigma_type[0]=1;}
if(strcmp(Sigma_readin, "free") == 0){Sigma_type[0]=2;}
//Sigma2
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "fixed") != 0 && strcmp(Sigma_readin, "free") != 0){printf("Error Sigma_nl type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&Sigma_nl_mean[1],&Sigma_nl_stddev[1]);
if(Sigma_nl_mean[1]<0 || Sigma_nl_stddev[1]<0){printf("Error with Sigma_nl values: %lf %lf. Exiting now...\n",Sigma_nl_mean[1],Sigma_nl_stddev[1]);exit(0);}
if(strcmp(Sigma_readin, "fixed") == 0){Sigma_type[1]=0;}
if(strcmp(Sigma_readin, "gaussian") == 0){Sigma_type[1]=1;}
if(strcmp(Sigma_readin, "free") == 0){Sigma_type[1]=2;}
//Sigma4
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "fixed") != 0 && strcmp(Sigma_readin, "free") != 0){printf("Error Sigma_nl type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&Sigma_nl_mean[2],&Sigma_nl_stddev[2]);
if(Sigma_nl_mean[2]<0 || Sigma_nl_stddev[2]<0){printf("Error with Sigma_nl values: %lf %lf. Exiting now...\n",Sigma_nl_mean[2],Sigma_nl_stddev[2]);exit(0);}
if(strcmp(Sigma_readin, "fixed") == 0){Sigma_type[2]=0;}
if(strcmp(Sigma_readin, "gaussian") == 0){Sigma_type[2]=1;}
if(strcmp(Sigma_readin, "free") == 0){Sigma_type[2]=2;}
//SigmaB0
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "fixed") != 0 && strcmp(Sigma_readin, "free") != 0){printf("Error Sigma_nl type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&Sigma_nl_mean[5],&Sigma_nl_stddev[5]);
if(Sigma_nl_mean[5]<0 || Sigma_nl_stddev[5]<0){printf("Error with Sigma_nl values: %lf %lf. Exiting now...\n",Sigma_nl_mean[5],Sigma_nl_stddev[5]);exit(0);}
if(strcmp(Sigma_readin, "fixed") == 0){Sigma_type[5]=0;}
if(strcmp(Sigma_readin, "gaussian") == 0){Sigma_type[5]=1;}
if(strcmp(Sigma_readin, "free") == 0){Sigma_type[5]=2;}

//Sigma_para
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "fixed") != 0 && strcmp(Sigma_readin, "free") != 0){printf("Error Sigma_nl type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&Sigma_nl_mean[3],&Sigma_nl_stddev[3]);
if(Sigma_nl_mean[3]<0 || Sigma_nl_stddev[3]<0){printf("Error with Sigma_nl values: %lf %lf. Exiting now...\n",Sigma_nl_mean[3],Sigma_nl_stddev[3]);exit(0);}
if(strcmp(Sigma_readin, "fixed") == 0){Sigma_type[3]=0;}
if(strcmp(Sigma_readin, "gaussian") == 0){Sigma_type[3]=1;}
if(strcmp(Sigma_readin, "free") == 0){Sigma_type[3]=2;}
//Sigma_perp
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "fixed") != 0 && strcmp(Sigma_readin, "free") != 0){printf("Error Sigma_nl type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n\n",&Sigma_nl_mean[4],&Sigma_nl_stddev[4]);
if(Sigma_nl_mean[4]<0 || Sigma_nl_stddev[4]<0){printf("Error with Sigma_nl values: %lf %lf. Exiting now...\n",Sigma_nl_mean[4],Sigma_nl_stddev[4]);exit(0);}
if(strcmp(Sigma_readin, "fixed") == 0){Sigma_type[4]=0;}
if(strcmp(Sigma_readin, "gaussian") == 0){Sigma_type[4]=1;}
if(strcmp(Sigma_readin, "free") == 0){Sigma_type[4]=2;}

if(strcmp(Sigma_def_type, "para-perp") == 0)
{
Sigma_nl_mean[0]=Sigma_nl_mean[3];
Sigma_nl_mean[1]=Sigma_nl_mean[4];

Sigma_nl_stddev[0]=Sigma_nl_stddev[3];
Sigma_nl_stddev[1]=Sigma_nl_stddev[4];

Sigma_type[0]=Sigma_type[3];
Sigma_type[1]=Sigma_type[4];
}


if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "BAOISO") == 0 ||  strcmp(type_of_analysis, "FSBAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0)
{
//Condition here that all fixed or all free/gaussian
if(strcmp(Sigma_def_type, "para-perp") == 0)
{
if(Sigma_type[0]==0 && Sigma_type[1]!=0){printf("Error. Both Sigma_para and Sigma_perp must be in either free/gaussian or fixed mode. Exiting now...\n");exit(0);}
if(Sigma_type[0]!=0 && Sigma_type[1]==0){printf("Error. Both Sigma_para and Sigma_perp must be in either free/gaussian or fixed mode. Exiting now...\n");exit(0);}

if( strcmp(Sigma_independent, "yes") == 0 && Sigma_type[0]==0){printf("Sigma_para and Sigma_perp set to be fixed to %lf and %lf, respectively.\n",Sigma_nl_mean[0],Sigma_nl_mean[1]);}
if( strcmp(Sigma_independent, "no") == 0 && Sigma_type[0]==0){printf("Sigma_para and Sigma_perp set to be fixed to %lf and %lf (through f=%lf), respectively.\n",Sigma_nl_mean[0],Sigma_nl_mean[0]/(1+ffactor),ffactor);}

if(Sigma_type[0]==2){printf("Sigma_para set to be free\n");}
if(Sigma_type[0]==1){printf("Sigma_para set to be free+gaussian prior (mu=%lf, sigma=%lf)\n",Sigma_nl_mean[0],Sigma_nl_stddev[0]);}

if( strcmp(Sigma_independent, "no") == 0 && Sigma_type[0]!=0)
{    
printf("Sigma_perp set to be Sigma_para/(1+f), where f=%lf.\n",ffactor);
}
if( strcmp(Sigma_independent, "yes") == 0 && Sigma_type[0]!=0){

if(Sigma_type[1]==2){printf("Sigma_perp set to be free\n");}
if(Sigma_type[1]==1){printf("Sigma_perp set to be free+gaussian prior (mu=%lf, sigma=%lf)\n",Sigma_nl_mean[1],Sigma_nl_stddev[1]);}

}


}
if(strcmp(Sigma_def_type, "effective") == 0)
{
if(Sigma_type[0]==0 && Sigma_type[1]!=0){printf("Error. Both Sigma_eff0, Sigma_eff2 and Sigma_eff4 must be in either free/gaussian or fixed mode. Exiting now...\n");exit(0);}
if(Sigma_type[0]==0 && Sigma_type[2]!=0){printf("Error. Both Sigma_eff0, Sigma_eff2 and Sigma_eff4 must be in either free/gaussian or fixed mode. Exiting now...\n");exit(0);}
if(Sigma_type[0]!=0 && Sigma_type[1]==0){printf("Error. Both Sigma_eff0, Sigma_eff2 and Sigma_eff4 must be in either free/gaussian or fixed mode. Exiting now...\n");exit(0);}
if(Sigma_type[1]!=0 && Sigma_type[2]==0){printf("Error. Both Sigma_eff0, Sigma_eff2 and Sigma_eff4 must be in either free/gaussian or fixed mode. Exiting now...\n");exit(0);}

if(Sigma_type[0]==0 &&  strcmp(Sigma_independent, "yes") == 0){printf("Sigma_0, Sigma_2, Sigma_4 set to be fixed to %lf, %lf, %lf, respectively\n",Sigma_nl_mean[0],Sigma_nl_mean[1],Sigma_nl_mean[2]);}
if(Sigma_type[0]==0 &&  strcmp(Sigma_independent, "no") == 0){printf("Sigma_0, Sigma_2, Sigma_4 set to be all fixed to %lf.\n",Sigma_nl_mean[0]);}

if(Sigma_type[0]==2){printf("Sigma_0 set to be free\n");}
if(Sigma_type[0]==1){printf("Sigma_0 set to be free+gaussian prior (mu=%lf, sigma=%lf)\n",Sigma_nl_mean[0],Sigma_nl_stddev[0]);}

if( strcmp(Sigma_independent, "no") == 0 && Sigma_type[0]!=0)
{
printf("Sigma_2 and Sigma_4 set to be equal to Sigma_0.\n");
}
if( strcmp(Sigma_independent, "yes") == 0 && Sigma_type[0]!=0){

if(Sigma_type[1]==2){printf("Sigma_2 set to be free\n");}
if(Sigma_type[1]==1){printf("Sigma_2 set to be free+gaussian prior (mu=%lf, sigma=%lf)\n",Sigma_nl_mean[1],Sigma_nl_stddev[1]);}

if(Sigma_type[2]==2){printf("Sigma_4 set to be free\n");}
if(Sigma_type[2]==1){printf("Sigma_4 set to be free+gaussian prior (mu=%lf, sigma=%lf)\n",Sigma_nl_mean[2],Sigma_nl_stddev[2]);}

}

}

if(strcmp(do_bispectrum,"yes") == 0){
if(Sigma_type[5]==0){printf("SigmaB_0 set to be fixed to %lf.\n",Sigma_nl_mean[5]);}
if(Sigma_type[5]==2){printf("SigmaB_0 set to be free\n");}
if(Sigma_type[5]==1){printf("SigmaB_0 set to be free+gaussian prior (mu=%lf, sigma=%lf)\n",Sigma_nl_mean[5],Sigma_nl_stddev[5]);}
}

}

//prior on FS variables
fscanf(f,"%*s %*s\n");
//Anoise
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "model") != 0 && strcmp(Sigma_readin, "data") != 0){printf("Error Anoise option choise %s\n. Exiting now...\n",Sigma_readin);exit(0);}
if(strcmp(Sigma_readin, "data") == 0){noise_option=1;}
if(strcmp(Sigma_readin, "model") == 0){noise_option=0;}
//Anoise
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "flat") != 0){printf("Error Anoise prior type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&FSprior_mean[0],&FSprior_stddev[0]);
if(Sigma_nl_stddev[0]<0){printf("Error with Anoise prior stdev<0 value: %lf. Exiting now...\n",FSprior_stddev[0]);exit(0);}
if(strcmp(Sigma_readin, "gaussian") == 0){FSprior_type[0]=1;}
if(strcmp(Sigma_readin, "flat") == 0){FSprior_type[0]=2;}
//b2
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "flat") != 0){printf("Error b2 prior type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&FSprior_mean[1],&FSprior_stddev[1]);
if(Sigma_nl_stddev[1]<0){printf("Error with Anoise prior stdev<0 value: %lf. Exiting now...\n",FSprior_stddev[1]);exit(0);}
if(strcmp(Sigma_readin, "gaussian") == 0){FSprior_type[1]=1;}
if(strcmp(Sigma_readin, "flat") == 0){FSprior_type[1]=2;}
//b2s2
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "flat") != 0){printf("Error b2s2 prior type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&FSprior_mean[2],&FSprior_stddev[2]);
if(Sigma_nl_stddev[2]<0){printf("Error with Anoise prior stdev<0 value: %lf. Exiting now...\n",FSprior_stddev[2]);exit(0);}
if(strcmp(Sigma_readin, "gaussian") == 0){FSprior_type[2]=1;}
if(strcmp(Sigma_readin, "flat") == 0){FSprior_type[2]=2;}
//b3nl
fscanf(f,"%*s %*s %s\n",Sigma_readin);
if(strcmp(Sigma_readin, "gaussian") != 0 && strcmp(Sigma_readin, "flat") != 0){printf("Error Anoise b3nl type %s\n. Exiting now...\n",Sigma_readin);exit(0);}
fscanf(f,"%*s %*s %lf %lf\n",&FSprior_mean[3],&FSprior_stddev[3]);
if(Sigma_nl_stddev[3]<0){printf("Error with Anoise prior stdev<0 value: %lf. Exiting now...\n",FSprior_stddev[3]);exit(0);}
if(strcmp(Sigma_readin, "gaussian") == 0){FSprior_type[3]=1;}
if(strcmp(Sigma_readin, "flat") == 0){FSprior_type[3]=2;}


fscanf(f,"%*s %*s %*s %lf %lf\n",&alpha_min,&alpha_max);
if(alpha_min>=alpha_max){printf("Error, alpha_min > alpha_max, %lf, %lf. Exiting now...\n",alpha_min,alpha_max);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",path_to_Plin);
fscanf(f,"%*s %*s %*s %s\n",path_to_Olin);
fscanf(f,"%*s %*s %*s %*s %s\n",type_BAO_fit);
if(strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "BAOANISO") == 0 ||  strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 )
{
g=fopen(path_to_Plin,"r");
if(g==NULL){printf("File %s not found. Exiting now...\n",path_to_Plin);exit(0);}
fclose(g);
g=fopen(path_to_Olin,"r");
if(g==NULL){printf("File %s not found. Exiting now...\n",path_to_Olin);exit(0);}
fclose(g);
sprintf(type_of_analysis_BAO,"yes");
}
g=fopen(path_to_Olin,"r");
if(g!=NULL){fclose(g);}

if(strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 || bispectrum_needs==1){
olin_exists=1;
}

if(strcmp(type_BAO_fit, "mcmc") != 0 && strcmp(type_BAO_fit, "analytic") != 0){printf("Error type BAO fit type %s\n. Exiting now...\n",type_BAO_fit);exit(0);}
if(strcmp(fit_BAO, "P24") == 0 && strcmp(type_BAO_fit, "analytic") == 0 && strcmp(path_to_mask1, "none") != 0){printf("I'd like to recomend you to use the monopole signal if quadrupole is also used.");}
if(strcmp(fit_BAO, "P04") == 0 && strcmp(type_BAO_fit, "analytic") == 0  && strcmp(path_to_mask1, "none") != 0){printf("I'd like to recomend you to use the quadrupole signal if hexadecapole is also used.\n");}

fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %d\n",&nthreads);
if(nthreads != omp_get_max_threads()){printf("The specified number of threads does not match the number of available threads found by OpenMP. Please make sure that the environment variable OMP_NUM_THREADS is set to the same value as in the paramfile.\n");exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",use_prop_cov);
if(strcmp(use_prop_cov, "no") != 0 && strcmp(use_prop_cov, "yes") != 0){printf("Error with answer to using prop. cov. %s\n. Exiting now...\n",use_prop_cov);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",path_to_cov);
fscanf(f,"%*s %*s %*s %*s %*s %*s %ld\n\n",&Nsteps);
fscanf(f,"%*s %*s %*s %*s %*s %lf\n\n",&step_size);
if(strcmp(use_prop_cov, "no") == 0){fraction = 10;}
else if(strcmp(use_prop_cov, "yes")==0){fraction = 1;}
if (Nsteps%(fraction*nthreads) != 0){
Nsteps += fraction*nthreads - Nsteps%(fraction*nthreads);
printf("The Number of steps has been increased to %ld, so they can be equally distributed among the %d threads.\n", Nsteps, nthreads);
}
if(strcmp(type_BAO_fit, "mcmc") == 0 &&  step_size<=0){printf("Error, step size must be positive: %lf. Exiting now...\n",step_size);exit(0);}
if(strcmp(type_BAO_fit, "mcmc") == 0 && step_size!=1.9){printf("Warning step size recommended to be 1.9, now %lf.\n",step_size);}

fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %lf\n\n",&alpha_step);

if(strcmp(type_BAO_fit, "analytic") == 0){
if(alpha_step<=0 || alpha_step>alpha_max-alpha_min){printf("Error, strange value for the step of alpha: %lf. Exiting now...\n",alpha_step);exit(0);}
}

fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %s\n",perturbation_theory_file);
if(strcmp(type_of_analysis, "FS") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 || strcmp(type_of_analysis, "FSalphasrecon") == 0 )
{
g=fopen(perturbation_theory_file,"r");
if(g==NULL){printf("File %s not found. Exiting now...\n",perturbation_theory_file);exit(0);}
fclose(g);
sprintf(type_of_analysis_FS,"yes");
}

fscanf(f,"%*s %*s %*s %*s %s\n",ptmodel_ps);
if( strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(ptmodel_ps, "linear") != 0 && strcmp(ptmodel_ps, "1L-SPT") != 0 && strcmp(ptmodel_ps, "2L-SPT") != 0 && strcmp(ptmodel_ps, "1L-RPT") != 0 && strcmp(ptmodel_ps, "2L-RPT") != 0){printf("Unvalid option for PT-Model %s. Available options: 'linear', '1L-SPT', '2L-SPT','1L-RPT', '2L-RPT'. Exiting now...\n",ptmodel_ps);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",rsdmodel_ps);
if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(rsdmodel_ps, "Kaiser87") != 0 && strcmp(rsdmodel_ps, "Scoccimarro04") != 0 && strcmp(rsdmodel_ps, "TNS10") != 0){printf("Unvalid option for RSD-Model '%s'. Available options: 'Kaiser87', 'Scoccimarro04', 'TNS10','1L-RPT', '2L-RPT'. Exiting now...\n",rsdmodel_ps);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",fogmodel_ps);
if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(fogmodel_ps, "Exponential") != 0 &&  strcmp(fogmodel_ps, "Lorentzian") != 0 && strcmp(fogmodel_ps, "Exponential_avir") != 0){printf("Unvalid option for FoG model %s. Available options: 'Lorentzian', 'Exponential', 'Exponential_avir'. Exiting now...\n",fogmodel_ps);exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n",ptmodel_bs);
if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(ptmodel_bs, "tree-level") != 0 &&  strcmp(ptmodel_bs, "1L-SPT") != 0 && strcmp(ptmodel_bs, "GilMarin14") != 0){printf("Unvalid option for FoG model %s. Available options: 'tree-level', '1L-SPT', 'GilMarin14'. Exiting now...\n",ptmodel_bs);exit(0);}
if(strcmp(ptmodel_bs, "GilMarin14") == 0 && strcmp(do_bispectrum,"yes")==0 && strcmp(type_of_analysis_FS,"yes") == 0  ){bispectrum_needs=1;}
if( strcmp(ptmodel_bs, "1L-SPT") == 0){printf("Sorry, 1L-SPT not available for Bispectrum yet!. Exiting now...\n");exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %s\n",local_b2s2);
if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(local_b2s2, "yes") != 0 && strcmp(local_b2s2, "no") != 0 && strcmp(local_b2s2, "off") != 0){printf("Error local lagrangian b2s2 bias has to be either off, yes or no: %s\n. Exiting now...\n",local_b2s2);exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %s\n",local_b3nl);
if( strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(local_b3nl, "yes") != 0 && strcmp(local_b3nl, "no") != 0 && strcmp(local_b3nl, "off") != 0){printf("Error local lagrangian b3nl bias has to be either off, yes or no: %s\n. Exiting now...\n",local_b3nl);exit(0);}
if( strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(local_b3nl, "no") == 0 && strcmp(do_power_spectrum,"no") == 0  ){printf("Warning. In the current bispectrum model there is no b3nl explicit dependence, and therefore, it doesn't make sense to have it as a free parameter. Please set local lagrangian of b3nl to yes or off, to avoid having an unconstrained free parameter. Exiting now...\n");exit(0);}
fscanf(f,"%*s %*s %s\n",RSD_fit);
if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(RSD_fit, "yes") != 0 && strcmp(RSD_fit, "no") != 0 && strcmp(RSD_fit, "shape") != 0 && strcmp(RSD_fit, "shape2") != 0){printf("Error RSD-fit has to be either yes, shape, shape2 or no: %s\n. Exiting now...\n",RSD_fit);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",sigma8_free);
if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(sigma8_free, "yes") != 0 && strcmp(sigma8_free, "no") != 0){printf("Error sigma8-free has to be either yes or no: %s\n. Exiting now...\n",sigma8_free);exit(0);}
fscanf(f,"%*s %*s %*s %*s %s\n",fog_free);
if( strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(fog_free, "yes") != 0 && strcmp(fog_free, "no") != 0){printf("Error FoG as free parameter has to be either yes or no: %s\n. Exiting now...\n",fog_free);exit(0);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n\n",fog_bs);
if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(fog_bs, "yes") != 0 && strcmp(fog_bs, "no") != 0){printf("Error Same FoG for bispectrum has to be either yes or no: %s\n. Exiting now...\n",fog_bs);exit(0);}

if(strcmp(type_of_analysis_FS, "yes") == 0 && strcmp(fog_bs, "yes") ==0 && strcmp(do_power_spectrum,"no") == 0){printf("Error, can't set bispectrum FoG to be the same as power spectrum if no power spectrum is included. Either enanble power spectrum or disable same FoG for P and B. Exiting now...\n");exit(0);}

fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",path_output);
fscanf(f,"%*s %*s %*s %s\n",identifier);
fscanf(f,"%*s %*s %*s %s\n",do_plot);
fclose(f);
if(strcmp(do_plot, "no") != 0 && strcmp(do_plot, "yes") != 0){printf("Error with answer about plotting. %s\n. Exiting now...\n",do_plot);exit(0);}


if(strcmp(type_BAO_fit, "analytic") == 0 && strcmp(do_bispectrum, "yes") == 0){printf("Error, the BAO analytic solver is not implemented to the bispectrum. Exiting now...\n");exit(0);}

if( Sigma_type[0]!= 0 && strcmp(type_BAO_fit, "analytic") == 0){printf("Warning. We recomend to use the mcmc option instead when exploring non-fixed Sigma_nl values. Exiting now...\n");exit(0);}


//errors and incompatibilities
//
if( strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0)
{
if( strcmp(fit_BAO, "P0") == 0 && strcmp(fit_RSD, "P0") == 0){printf("Warning. Can't do an anisotropic type of fit (%s) with just one multipole. Exiting now...\n",type_of_analysis);exit(0);}

if( strcmp(fit_BAO, "P2") == 0 && strcmp(fit_RSD, "P2") == 0 && strcmp(do_bispectrum,"yes") == 0){printf("Warning. Can't do an anisotropic type of fit (%s) with just one multipole. Exiting now...\n",type_of_analysis);exit(0);}

if( strcmp(fit_BAO, "P4") == 0 && strcmp(fit_RSD, "P4") == 0 && strcmp(do_bispectrum,"yes") == 0){printf("Warning. Can't do an anisotropic type of fit (%s) with just one multipole. Exiting now...\n",type_of_analysis);exit(0);}

if( strcmp(type_BAO_fit, "analytic") == 0){printf("Warning. Can't do an analytic fit with BAOANISO or FSBAOANISO template. Exiting now...\n");exit(0);}
if( strcmp(Sigma_def_type, "effective") == 0){printf("Warning. Can't choose effective variables within with the  BAOANISO template. Exiting now...\n");exit(0);}
}
if( strcmp(type_of_analysis, "FSBAOISO") == 0)
{
if( strcmp(type_BAO_fit, "analytic") == 0){printf("Warning. Can't do an analytic fit with FSBAOISO template. Exiting now...\n");exit(0);}
}

if( strcmp(RSD_fit,"no") == 0)
{
if( strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 || strcmp(type_of_analysis, "FSalphasrecon") == 0){printf("Warning. Can't do a combined FS + BAO with RSD option set to no. Exitting now...\n");exit(0);}
}

if( strcmp(type_of_analysis, "FSalphasrecon") == 0 && strcmp(do_power_spectrum, "no") == 0)
{
printf("Error, FSalphasrecon option requires power spectrum enabled (otherwise isotropic fit!). Exitin now...\n");exit(0);//may change this if alpha's recon provide constrains on aiso, but we'll see...
}

//Warning about same output as input for prop. covariance when set to yes. 
if(strcmp(type_of_analysis, "BAOISO") == 0){sprintf(name_file,"%s/mcmcBAOISO_output_%s.txt",path_output,identifier);}
if(strcmp(type_of_analysis, "BAOANISO") == 0){sprintf(name_file,"%s/mcmcBAOANISO_output_%s.txt",path_output,identifier);}
if(strcmp(type_of_analysis, "FS") == 0){sprintf(name_file,"%s/mcmcFS_output_%s.txt",path_output,identifier);}
if(strcmp(type_of_analysis, "FSBAOISO") == 0){sprintf(name_file,"%s/mcmcFSBAOISO_output_%s.txt",path_output,identifier);}
if(strcmp(type_of_analysis, "FSBAOANISO") == 0){sprintf(name_file,"%s/mcmcFSBAOANISO_output_%s.txt",path_output,identifier);}
if(strcmp(type_of_analysis, "FSalphasrecon") == 0){sprintf(name_file,"%s/mcmcFSalphas_output_%s.txt",path_output,identifier);}

if(strcmp(path_to_cov, name_file) == 0 &&  strcmp(use_prop_cov, "yes") == 0){printf("Warning, output has the same name as the input file for the proposal covariance. Exiting now...\n");exit(0);}




if(strcmp(do_power_spectrum, "yes") == 0){

if(strcmp(type_of_analysis, "FSalphasrecon") == 0)
{
NeffP0bao=2;
NeffP2bao=0;
NeffP4bao=0;
P0bao = (double*) calloc(NeffP0bao, sizeof(double));
get_data_alphas(path_to_data1_bao,P0bao);
//printf("%lf %lf %s\n",P0bao[0],P0bao[1],type_of_analysis_BAO);exit(0);
}else{P0bao=NULL;}

if( strcmp(type_of_analysis_BAO ,"yes") ==0)
{

NeffP0bao=countPk(0,path_to_data1_bao,kminP0bao,kmaxP0bao);printf("Neff P0 BAO: %d\n",NeffP0bao);
NeffP2bao=countPk(0,path_to_data1_bao,kminP2bao,kmaxP2bao);printf("Neff P2 BAO: %d\n",NeffP2bao);
NeffP4bao=countPk(0,path_to_data1_bao,kminP4bao,kmaxP4bao);printf("Neff P4 BAO: %d\n",NeffP4bao);

                k0bao = (double*) calloc(NeffP0bao, sizeof(double));
                k2bao = (double*) calloc(NeffP2bao, sizeof(double));
                k4bao = (double*) calloc(NeffP4bao, sizeof(double));
                kav0bao = (double*) calloc(NeffP0bao, sizeof(double));
                kav2bao = (double*) calloc(NeffP2bao, sizeof(double));
                kav4bao = (double*) calloc(NeffP4bao, sizeof(double));
	        P0bao = (double*) calloc(NeffP0bao, sizeof(double));
                P2bao = (double*) calloc(NeffP2bao, sizeof(double));
                P4bao = (double*) calloc(NeffP4bao, sizeof(double));
                errP0bao = (double*) calloc(NeffP0bao, sizeof(double));
                errP2bao = (double*) calloc(NeffP2bao, sizeof(double));
                errP4bao = (double*) calloc(NeffP4bao, sizeof(double));


get_data(path_to_data1_bao,k0bao,kav0bao,P0bao,k2bao,kav2bao,P2bao,k4bao,kav4bao,P4bao,params,kminP0bao,kmaxP0bao,kminP2bao,kmaxP2bao,kminP4bao,kmaxP4bao,type_of_analysis,0);
Pnoise_bao=params[0];if(Pnoise_bao==0){printf("Warning, Pnoise for BAO is read to be 0.\n");}
sumw_bao=params[1];
I22_bao=params[2];
determine_spacing(spacing_dataNGC_bao,kav0bao,kav2bao,kav4bao,NeffP0bao,NeffP2bao,NeffP4bao);

}
if( strcmp(type_of_analysis_FS ,"yes") ==0)
{

NeffP0rsd=countPk(0,path_to_data1_rsd,kminP0rsd,kmaxP0rsd);printf("Neff P0 FS: %d\n",NeffP0rsd);
NeffP2rsd=countPk(0,path_to_data1_rsd,kminP2rsd,kmaxP2rsd);printf("Neff P2 FS: %d\n",NeffP2rsd);
NeffP4rsd=countPk(0,path_to_data1_rsd,kminP4rsd,kmaxP4rsd);printf("Neff P4 FS: %d\n",NeffP4rsd);

                k0rsd = (double*) calloc(NeffP0rsd, sizeof(double));
                k2rsd = (double*) calloc(NeffP2rsd, sizeof(double));
                k4rsd = (double*) calloc(NeffP4rsd, sizeof(double));
                kav0rsd = (double*) calloc(NeffP0rsd, sizeof(double));
                kav2rsd = (double*) calloc(NeffP2rsd, sizeof(double));
                kav4rsd = (double*) calloc(NeffP4rsd, sizeof(double));
                P0rsd = (double*) calloc(NeffP0rsd, sizeof(double));
                P2rsd = (double*) calloc(NeffP2rsd, sizeof(double));
                P4rsd = (double*) calloc(NeffP4rsd, sizeof(double));
                errP0rsd = (double*) calloc(NeffP0rsd, sizeof(double));
                errP2rsd = (double*) calloc(NeffP2rsd, sizeof(double));
                errP4rsd = (double*) calloc(NeffP4rsd, sizeof(double));


get_data(path_to_data1_rsd,k0rsd,kav0rsd,P0rsd,k2rsd,kav2rsd,P2rsd,k4rsd,kav4rsd,P4rsd,params,kminP0rsd,kmaxP0rsd,kminP2rsd,kmaxP2rsd,kminP4rsd,kmaxP4rsd,type_of_analysis,1);
Pnoise_rsd=params[0];if(Pnoise_rsd==0){printf("Warning, Pnoise for RSD is read to be 0.\n");}
sumw_rsd=params[1];
I22_rsd=params[2];
determine_spacing(spacing_dataNGC_rsd,kav0rsd,kav2rsd,kav4rsd,NeffP0rsd,NeffP2rsd,NeffP4rsd);

}

}
//mask has to be available for power spectrum=no
if(strcmp(path_to_mask1, "none") != 0){
Nmask=countlines(path_to_mask1);

                posAV = (double*) calloc(Nmask, sizeof(double));
                pos = (double*) calloc(Nmask, sizeof(double));
                W0 = (double*) calloc(Nmask, sizeof(double));
                W2 = (double*) calloc(Nmask, sizeof(double));
                W4 = (double*) calloc(Nmask, sizeof(double));
                W6 = (double*) calloc(Nmask, sizeof(double));
                W8 = (double*) calloc(Nmask, sizeof(double));

//if possible renormalize on RSD data
if( strcmp(type_of_analysis_FS ,"yes") ==0){
params[0]=Pnoise_rsd;
params[1]=sumw_rsd;
params[2]=I22_rsd;
}
if( strcmp(type_of_analysis_FS ,"no") ==0){
params[0]=Pnoise_bao;
params[1]=sumw_bao;
params[2]=I22_bao;
}

get_mask(path_to_mask1,posAV,pos,W0,W2,W4,W6,W8,Nmask,type_of_analysis,params,mask_renormalization);
Nmask=(int)(params[0]);//for(i=0;i<Nmask;i++){printf("%lf %lf %e %e %e %e %e\n",posAV[i],pos[i],W0[i],W2[i],W4[i],W6[i],W8[i]);}
determine_spacing_theo(spacing_maskNGC,posAV,Nmask,2);
}



if(strcmp(do_bispectrum, "yes") == 0){

if( strcmp(type_of_analysis_BAO ,"yes") ==0)
{
                NeffB0bao=countPk(1,path_to_data1_bis_bao,kminB0bao,kmaxB0bao);printf("Neff B0-BAO: %d\n",NeffB0bao);
                B0bao = (double*) calloc(NeffB0bao, sizeof(double));
                Bnoise_bao = (double*) calloc(NeffB0bao, sizeof(double));
                errB0bao = (double*) calloc(NeffB0bao, sizeof(double));
                k11bao = (double*) calloc(NeffB0bao, sizeof(double));
                k22bao = (double*) calloc(NeffB0bao, sizeof(double));
                k33bao = (double*) calloc(NeffB0bao, sizeof(double));

get_data_bis(path_to_data1_bis_bao,k11bao,k22bao,k33bao,B0bao,Bnoise_bao,kminB0bao,kmaxB0bao,bispectrum_BQ);
}

if( strcmp(type_of_analysis_FS ,"yes") ==0)
{
                NeffB0rsd=countPk(1,path_to_data1_bis_rsd,kminB0rsd,kmaxB0rsd);printf("Neff B0-FS: %d\n",NeffB0rsd);
                B0rsd = (double*) calloc(NeffB0rsd, sizeof(double));
                Bnoise_rsd = (double*) calloc(NeffB0rsd, sizeof(double));
                errB0rsd = (double*) calloc(NeffB0rsd, sizeof(double));
                k11rsd = (double*) calloc(NeffB0rsd, sizeof(double));
                k22rsd = (double*) calloc(NeffB0rsd, sizeof(double));
                k33rsd = (double*) calloc(NeffB0rsd, sizeof(double));

get_data_bis(path_to_data1_bis_rsd,k11rsd,k22rsd,k33rsd,B0rsd,Bnoise_rsd,kminB0rsd,kmaxB0rsd,bispectrum_BQ);
}

}

bao=0;rsd=0;
if( strcmp(type_of_analysis,"FSBAOISO") ==0 || strcmp(type_of_analysis,"FSBAOANISO") ==0 && strcmp(type_of_analysis,"FSalphasrecon") ==0){bao=1;rsd=1;}
if( strcmp(type_of_analysis,"BAOISO") ==0 || strcmp(type_of_analysis,"BAOANISO") ==0 ){bao=1;rsd=0;}
if( strcmp(type_of_analysis,"FS") ==0){rsd=1;bao=0;}

modeP0bao=0;modeP2bao=0;modeP4bao=0;modeB0bao=0;

if( strcmp(do_power_spectrum,"yes") == 0){
if( strcmp(fit_BAO,"P0") == 0 || strcmp(fit_BAO,"P02") == 0 || strcmp(fit_BAO,"P04") == 0 || strcmp(fit_BAO,"P024") == 0){modeP0bao=1;}
if( strcmp(fit_BAO,"P2") == 0 || strcmp(fit_BAO,"P02") == 0 || strcmp(fit_BAO,"P24") == 0 || strcmp(fit_BAO,"P024") == 0){modeP2bao=1;}
if( strcmp(fit_BAO,"P4") == 0 || strcmp(fit_BAO,"P04") == 0 || strcmp(fit_BAO,"P24") == 0 || strcmp(fit_BAO,"P024") == 0){modeP4bao=1;}
}
if( strcmp(do_bispectrum,"yes") == 0){
modeB0bao=1;
}

input_elements_data_vector=bao*(NeffP0bao*modeP0bao+NeffP2bao*modeP2bao+NeffP4bao*modeP4bao+NeffB0bao*modeB0bao);

modeP0rsd=0;modeP2rsd=0;modeP4rsd=0;modeB0rsd=0;
if( strcmp(do_power_spectrum,"yes") == 0){
if( strcmp(fit_RSD,"P0") == 0 || strcmp(fit_RSD,"P02") == 0 || strcmp(fit_RSD,"P04") == 0 || strcmp(fit_RSD,"P024") == 0){modeP0rsd=1;}
if( strcmp(fit_RSD,"P2") == 0 || strcmp(fit_RSD,"P02") == 0 || strcmp(fit_RSD,"P24") == 0 || strcmp(fit_RSD,"P024") == 0){modeP2rsd=1;}
if( strcmp(fit_RSD,"P4") == 0 || strcmp(fit_RSD,"P04") == 0 || strcmp(fit_RSD,"P24") == 0 || strcmp(fit_RSD,"P024") == 0){modeP4rsd=1;}
}
if( strcmp(do_bispectrum,"yes") == 0){
modeB0rsd=1;
}

//Linear Compression Matrix
if( strcmp(do_compression,"linear") == 0)
{

input_elements_data_vector=input_elements_data_vector+rsd*(NeffP0rsd*modeP0rsd+NeffP2rsd*modeP2rsd+NeffP4rsd*modeP4rsd+NeffB0rsd*modeB0rsd);

output_elements_data_vector=get_compression_file(path_to_compression,matrix_linear_compression,input_elements_data_vector);

                 matrix_linear_compression = (double **) calloc(output_elements_data_vector, sizeof(double*));//rows
                 for(i=0;i<output_elements_data_vector;i++){matrix_linear_compression[i]= (double *) calloc(input_elements_data_vector, sizeof(double));}//input
//Matrix_{rows,col} = Matrix_{output, input}

   check=get_compression_file(path_to_compression,matrix_linear_compression,input_elements_data_vector);
   if(check==0){printf("Linear compression matrix succesfully read/\n");}
   else{printf("Error reading %s matrix. Exiting now...\n",path_to_compression);exit(0);}

}


//Covariance

Ncov=0;NcovP=0;
if(strcmp(do_power_spectrum, "yes") == 0){

if( strcmp(type_of_analysis_BAO ,"yes") == 0)
{
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){Ncov=Ncov+NeffP0bao;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){Ncov=Ncov+NeffP2bao;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){Ncov=Ncov+NeffP4bao;}
NcovP=Ncov;
if(strcmp(do_bispectrum, "yes") == 0){Ncov=Ncov+NeffB0bao;}
}
if( strcmp(type_of_analysis_FS ,"yes") == 0)
{
if(strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0 ){Ncov=Ncov+NeffP0rsd;}
if(strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0){Ncov=Ncov+NeffP2rsd;}
if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0){Ncov=Ncov+NeffP4rsd;}
NcovP=Ncov;
if(strcmp(do_bispectrum, "yes") == 0){Ncov=Ncov+NeffB0rsd;}
}

if( strcmp(type_of_analysis ,"FSalphasrecon") == 0){NcovP=NcovP+NeffP0bao;Ncov=Ncov+NeffP0bao;}

}
else{//if power spectrum==no, then bispectrum must be yes (otherwise no and no, is not allowed, check)
if( strcmp(type_of_analysis_BAO ,"yes") == 0){Ncov=Ncov+NeffB0bao;}
if( strcmp(type_of_analysis_FS ,"yes") == 0){Ncov=Ncov+NeffB0rsd;}
//recall, for bispectrum, either BAO, or FS, but not both (same catalogue, can't be reconstructed bispectrum)
}

 
                cov= (double*) calloc( Ncov*Ncov, sizeof(double));
printf("Ncov=%d, NcovP=%d\n",Ncov,NcovP);

if( strcmp(type_of_analysis_BAO ,"yes") == 0){path_to_mocks1_ending = strrchr(path_to_mocks1_bao, '.');} //this function sets a pointer to the last occurence of a character in a string. Used to identify the filetype.
if( strcmp(type_of_analysis_FS ,"yes") == 0){path_to_mocks1_ending = strrchr(path_to_mocks1_rsd, '.');} //this function sets a pointer to the last occurence of a character in a string. Used to identify the filetype.


//new
if (path_to_mocks1_ending != NULL){
//printf("here\n");exit(0);
  if ((strcmp(path_to_mocks1_ending, ".cov") == 0)){

if( strcmp(type_of_analysis_BAO ,"yes") == 0 && strcmp(type_of_analysis_FS ,"yes") == 0){printf("Error, no file-covariance for BAO+FS type of analysis. Exiting now...\n");exit(0);}    //need to fix this

if( strcmp(type_of_analysis ,"FSalphasrecon") == 0){printf("Error, no file-covariance for alphasrec+FS type of analysis. Exiting now...\n");exit(0);}    //need to fix this

if( strcmp(type_of_analysis_BAO ,"yes") == 0){ get_cov_from_file(path_to_mocks1_bao,path_to_mocks1_bis_bao,cov,Ncov,Nrealizations, NeffP0bao,NeffP0rsd, NeffP2bao,NeffP2rsd, NeffP4bao,NeffP4rsd, NeffB0bao,NeffB0rsd,errP0bao,errP0rsd,errP2bao,errP2rsd,errP4bao,errP4rsd,errB0bao,errB0rsd,do_power_spectrum, do_bispectrum,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction);Nrealizations_used_NGC=Nrealizations;}

if( strcmp(type_of_analysis_FS ,"yes") == 0){ get_cov_from_file(path_to_mocks1_rsd,path_to_mocks1_bis_rsd,cov,Ncov,Nrealizations, NeffP0bao,NeffP0rsd, NeffP2bao,NeffP2rsd, NeffP4bao,NeffP4rsd, NeffB0bao,NeffB0rsd,errP0bao,errP0rsd,errP2bao,errP2rsd,errP4bao,errP4rsd,errB0bao,errB0rsd,do_power_spectrum, do_bispectrum,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction);Nrealizations_used_NGC=Nrealizations;}

  }
  else {
    sprintf(covfile1,"%s/%s_1",path_output,identifier);

if( strcmp(type_of_analysis ,"FSalphasrecon") == 0 && strcmp(covarianceFSa_option,"varying") == 0 ){getBAOcov(paramsBAO,path_to_data1_bao);}//get BAOpostcov
else{paramsBAO[0]=-1;paramsBAO[1]=-1;paramsBAO[2]=-1;}

    Nrealizations_used_NGC=get_cov_from_mocks(path_to_mocks1_bao,path_to_mocks1_rsd,path_to_mocks1_bis_bao,path_to_mocks1_bis_rsd,covfile1,cov,Ncov,Nrealizations, NeffP0bao,NeffP0rsd, NeffP2bao,NeffP2rsd, NeffP4bao,NeffP4rsd, NeffB0bao,NeffB0rsd,errP0bao,errP0rsd,errP2bao,errP2rsd,errP4bao,errP4rsd,errB0bao,errB0rsd,kminP0bao,kminP0rsd,kmaxP0bao,kmaxP0rsd,kminP2bao,kminP2rsd,kmaxP2bao,kmaxP2rsd,kminP4bao,kminP4rsd,kmaxP4bao,kmaxP4rsd,kminB0bao,kminB0rsd,kmaxB0bao,kmaxB0rsd,type_of_analysis,fit_BAO,fit_RSD,do_power_spectrum, do_bispectrum,bispectrum_BQ,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction,covarianceFSa_option,paramsBAO);
  }
}
else{
//printf("here2\n");exit(0);
  sprintf(covfile1,"%s/%s_1",path_output,identifier);
  printf("Covfile = %s\n",covfile1);

if( strcmp(type_of_analysis ,"FSalphasrecon") == 0 && strcmp(covarianceFSa_option,"varying") == 0 ){getBAOcov(paramsBAO,path_to_data1_bao);}//get BAOpostcov
else{paramsBAO[0]=-1;paramsBAO[1]=-1;paramsBAO[2]=-1;}

  Nrealizations_used_NGC=get_cov_from_mocks(path_to_mocks1_bao,path_to_mocks1_rsd,path_to_mocks1_bis_bao,path_to_mocks1_bis_rsd,covfile1,cov,Ncov,Nrealizations, NeffP0bao,NeffP0rsd, NeffP2bao,NeffP2rsd, NeffP4bao,NeffP4rsd, NeffB0bao,NeffB0rsd,errP0bao,errP0rsd,errP2bao,errP2rsd,errP4bao,errP4rsd,errB0bao,errB0rsd,kminP0bao,kminP0rsd,kmaxP0bao,kmaxP0rsd,kminP2bao,kminP2rsd,kmaxP2bao,kmaxP2rsd,kminP4bao,kminP4rsd,kmaxP4bao,kmaxP4rsd,kminB0bao,kminB0rsd,kmaxB0bao,kmaxB0rsd,type_of_analysis,fit_BAO,fit_RSD,do_power_spectrum, do_bispectrum,bispectrum_BQ,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction,covarianceFSa_option,paramsBAO);
 
}

//exit(0);

if(Nchunks==2)//both chunks considered independent
{

if(strcmp(do_power_spectrum, "yes") == 0){

if(strcmp(type_of_analysis, "FSalphasrecon") == 0)
{
NeffP0baoSGC=2;
NeffP2baoSGC=0;
NeffP4baoSGC=0;
P0baoSGC = (double*) calloc(NeffP0baoSGC, sizeof(double));
get_data_alphas(path_to_data2_bao,P0baoSGC);
}else{P0baoSGC=NULL;}

if( strcmp(type_of_analysis_BAO ,"yes") ==0)
{

//ensure NeffP0SGC is equal to Neff

NeffP0baoSGC=countPk(0,path_to_data2_bao,kminP0bao,kmaxP0bao);if(NeffP0bao!=NeffP0baoSGC){printf("Warning, Neff for NGC and SGC are different for P0-BAO.\n");}
NeffP2baoSGC=countPk(0,path_to_data2_bao,kminP2bao,kmaxP2bao);if(NeffP2bao!=NeffP2baoSGC){printf("Warning, Neff for NGC and SGC are different for P2-BAO.\n");}
NeffP4baoSGC=countPk(0,path_to_data2_bao,kminP4bao,kmaxP4bao);if(NeffP4bao!=NeffP4baoSGC){printf("Warning, Neff for NGC and SGC are different for P4-BAO.\n");}


                k0baoSGC = (double*) calloc(NeffP0baoSGC, sizeof(double));
                k2baoSGC = (double*) calloc(NeffP2baoSGC, sizeof(double));
                k4baoSGC = (double*) calloc(NeffP4baoSGC, sizeof(double));
                kav0baoSGC = (double*) calloc(NeffP0baoSGC, sizeof(double));
                kav2baoSGC = (double*) calloc(NeffP2baoSGC, sizeof(double));
                kav4baoSGC = (double*) calloc(NeffP4baoSGC, sizeof(double));
                P0baoSGC = (double*) calloc(NeffP0baoSGC, sizeof(double));
                P2baoSGC = (double*) calloc(NeffP2baoSGC, sizeof(double));
                P4baoSGC = (double*) calloc(NeffP4baoSGC, sizeof(double));
                errP0baoSGC = (double*) calloc(NeffP0baoSGC, sizeof(double));
                errP2baoSGC = (double*) calloc(NeffP2baoSGC, sizeof(double));
                errP4baoSGC = (double*) calloc(NeffP4baoSGC, sizeof(double));


get_data(path_to_data2_bao,k0baoSGC,kav0baoSGC,P0baoSGC,k2baoSGC,kav2baoSGC,P2baoSGC,k4baoSGC,kav4baoSGC,P4baoSGC,params,kminP0bao,kmaxP0bao,kminP2bao,kmaxP2bao,kminP4bao,kmaxP4bao,type_of_analysis,0);
Pnoise_baoSGC=params[0];if(Pnoise_baoSGC==0){printf("Warning, Pnoise for BAO-SGC is read to be 0.\n");}
sumw_bao=params[1];
I22_bao=params[2];

determine_spacing(spacing_dataSGC_bao,kav0baoSGC,kav2baoSGC,kav4baoSGC,NeffP0baoSGC,NeffP2baoSGC,NeffP4baoSGC);

if( strcmp(spacing_dataNGC_bao,spacing_dataSGC_bao) != 0){printf("Warning, data in NGC and SGC, spaced differently\n");}

}

if( strcmp(type_of_analysis_FS ,"yes") ==0)
{

NeffP0rsdSGC=countPk(0,path_to_data2_rsd,kminP0rsd,kmaxP0rsd);if(NeffP0rsd!=NeffP0rsdSGC){printf("Warning, Neff for NGC and SGC are different for P0-RSD.\n");}
NeffP2rsdSGC=countPk(0,path_to_data2_rsd,kminP2rsd,kmaxP2rsd);if(NeffP2rsd!=NeffP2rsdSGC){printf("Warning, Neff for NGC and SGC are different for P2-RSD.\n");}
NeffP4rsdSGC=countPk(0,path_to_data2_rsd,kminP4rsd,kmaxP4rsd);if(NeffP4rsd!=NeffP4rsdSGC){printf("Warning, Neff for NGC and SGC are different for P4-RSD.\n");}


                k0rsdSGC = (double*) calloc(NeffP0rsdSGC, sizeof(double));
                k2rsdSGC = (double*) calloc(NeffP2rsdSGC, sizeof(double));
                k4rsdSGC = (double*) calloc(NeffP4rsdSGC, sizeof(double));
                kav0rsdSGC = (double*) calloc(NeffP0rsdSGC, sizeof(double));
                kav2rsdSGC = (double*) calloc(NeffP2rsdSGC, sizeof(double));
                kav4rsdSGC = (double*) calloc(NeffP4rsdSGC, sizeof(double));
                P0rsdSGC = (double*) calloc(NeffP0rsdSGC, sizeof(double));
                P2rsdSGC = (double*) calloc(NeffP2rsdSGC, sizeof(double));
                P4rsdSGC = (double*) calloc(NeffP4rsdSGC, sizeof(double));
                errP0rsdSGC = (double*) calloc(NeffP0rsdSGC, sizeof(double));
                errP2rsdSGC = (double*) calloc(NeffP2rsdSGC, sizeof(double));
                errP4rsdSGC = (double*) calloc(NeffP4rsdSGC, sizeof(double));


get_data(path_to_data2_rsd,k0rsdSGC,kav0rsdSGC,P0rsdSGC,k2rsdSGC,kav2rsdSGC,P2rsdSGC,k4rsdSGC,kav4rsdSGC,P4rsdSGC,params,kminP0rsd,kmaxP0rsd,kminP2rsd,kmaxP2rsd,kminP4rsd,kmaxP4rsd,type_of_analysis,1);
Pnoise_rsdSGC=params[0];if(Pnoise_rsdSGC==0){printf("Warning, Pnoise for RSD-SGC is read to be 0.\n");}
sumw_rsd=params[1];
I22_rsd=params[2];

determine_spacing(spacing_dataSGC_rsd,kav0rsdSGC,kav2rsdSGC,kav4rsdSGC,NeffP0rsdSGC,NeffP2rsdSGC,NeffP4rsdSGC);

if( strcmp(spacing_dataNGC_rsd,spacing_dataSGC_rsd) != 0){printf("Warning, data in NGC and SGC, spaced differently\n");}


}

}
//mask has to be available for power spectrum=0 as well (except for normalization)
if(strcmp(path_to_mask2, "none") != 0){
NmaskSGC=countlines(path_to_mask2);

                posAVSGC = (double*) calloc(NmaskSGC, sizeof(double));
                posSGC = (double*) calloc(NmaskSGC, sizeof(double));
                W0SGC = (double*) calloc(NmaskSGC, sizeof(double));
                W2SGC = (double*) calloc(NmaskSGC, sizeof(double));
                W4SGC = (double*) calloc(NmaskSGC, sizeof(double));
                W6SGC = (double*) calloc(NmaskSGC, sizeof(double));
                W8SGC = (double*) calloc(NmaskSGC, sizeof(double));

if( strcmp(type_of_analysis_FS ,"yes") ==0){
params[0]=Pnoise_rsd;
params[1]=sumw_rsd;
params[2]=I22_rsd;
}
if( strcmp(type_of_analysis_FS ,"no") ==0){
params[0]=Pnoise_bao;
params[1]=sumw_bao;
params[2]=I22_bao;
}


get_mask(path_to_mask2,posAVSGC,posSGC,W0SGC,W2SGC,W4SGC,W6SGC,W8SGC,NmaskSGC,type_of_analysis,params,mask_renormalization);
NmaskSGC=(int)(params[0]);

determine_spacing_theo(spacing_maskSGC,posAVSGC,NmaskSGC,3);

}


if(strcmp(do_bispectrum, "yes") == 0){

if( strcmp(type_of_analysis_BAO ,"yes") == 0)
{

      NeffB0baoSGC=countPk(1,path_to_data2_bis_bao,kminB0bao,kmaxB0bao);if(NeffB0bao!=NeffB0baoSGC){printf("Warning, Neff for NGC and SGC are different for B0.\n");}

                B0baoSGC = (double*) calloc(NeffB0baoSGC, sizeof(double));
                Bnoise_baoSGC = (double*) calloc(NeffB0baoSGC, sizeof(double));
                errB0baoSGC = (double*) calloc(NeffB0baoSGC, sizeof(double));
                k11baoSGC = (double*) calloc(NeffB0baoSGC, sizeof(double));
                k22baoSGC = (double*) calloc(NeffB0baoSGC, sizeof(double));
                k33baoSGC = (double*) calloc(NeffB0baoSGC, sizeof(double));

get_data_bis(path_to_data2_bis_bao,k11baoSGC,k22baoSGC,k33baoSGC,B0baoSGC,Bnoise_baoSGC,kminB0bao,kmaxB0bao,bispectrum_BQ);
}
if( strcmp(type_of_analysis_FS ,"yes") == 0)
{
      NeffB0rsdSGC=countPk(1,path_to_data2_bis_rsd,kminB0rsd,kmaxB0rsd);if(NeffB0rsd!=NeffB0rsdSGC){printf("Warning, Neff for NGC and SGC are different for B0.\n");}

                B0rsdSGC = (double*) calloc(NeffB0rsdSGC, sizeof(double));
                Bnoise_rsdSGC = (double*) calloc(NeffB0rsdSGC, sizeof(double));
                errB0rsdSGC = (double*) calloc(NeffB0rsdSGC, sizeof(double));
                k11rsdSGC = (double*) calloc(NeffB0rsdSGC, sizeof(double));
                k22rsdSGC = (double*) calloc(NeffB0rsdSGC, sizeof(double));
                k33rsdSGC = (double*) calloc(NeffB0rsdSGC, sizeof(double));

get_data_bis(path_to_data2_bis_rsd,k11rsdSGC,k22rsdSGC,k33rsdSGC,B0rsdSGC,Bnoise_rsdSGC,kminB0rsd,kmaxB0rsd,bispectrum_BQ);


}
}
                 covSGC= (double*) calloc( Ncov*Ncov, sizeof(double));



if( strcmp(type_of_analysis_BAO ,"yes") == 0){path_to_mocks2_ending = strrchr(path_to_mocks2_bao, '.');} //this function sets a pointer to the last occurence of a character in a string. Used to identify the filetype.
if( strcmp(type_of_analysis_FS ,"yes") == 0){path_to_mocks2_ending = strrchr(path_to_mocks2_rsd, '.');} //this function sets a pointer to the last occurence of a character in a string. Used to identify the filetype.

// read the covariance from a file, if the extension is ".cov", otherwise "path_to_mocks_ending" is  NULL Pointer and the covariance is read from mocks.
if (path_to_mocks2_ending != NULL){
  if ((strcmp(path_to_mocks2_ending, ".cov") == 0)){

if( strcmp(type_of_analysis_BAO ,"yes") == 0 && strcmp(type_of_analysis_FS ,"yes") == 0){printf("Error, no file-covariance for BAO+FS type of analysis. Exiting now...\n");exit(0);}    //need to fix this

if( strcmp(type_of_analysis ,"FSalphasrecon") == 0){printf("Error, no file-covariance for alphasrec+FS type of analysis. Exiting now...\n");exit(0);}    //need to fix this

if( strcmp(type_of_analysis_BAO ,"yes") == 0){get_cov_from_file(path_to_mocks2_bao,path_to_mocks2_bis_bao,covSGC,Ncov,Nrealizations, NeffP0baoSGC,NeffP0rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, NeffB0baoSGC,NeffB0rsdSGC,errP0baoSGC,errP0rsdSGC,errP2baoSGC,errP2rsdSGC,errP4baoSGC,errP4rsdSGC,errB0baoSGC,errB0rsdSGC,do_power_spectrum, do_bispectrum,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction);Nrealizations_used_SGC=Nrealizations;}

if( strcmp(type_of_analysis_FS ,"yes") == 0){get_cov_from_file(path_to_mocks2_rsd,path_to_mocks2_bis_rsd,covSGC,Ncov,Nrealizations, NeffP0baoSGC,NeffP0rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, NeffB0baoSGC,NeffB0rsdSGC,errP0baoSGC,errP0rsdSGC,errP2baoSGC,errP2rsdSGC,errP4baoSGC,errP4rsdSGC,errB0baoSGC,errB0rsdSGC,do_power_spectrum, do_bispectrum,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction);Nrealizations_used_SGC=Nrealizations;}


  }  
  else {
    sprintf(covfile2,"%s/%s_2",path_output,identifier);

if( strcmp(type_of_analysis ,"FSalphasrecon") == 0 && strcmp(covarianceFSa_option,"varying") == 0 ){getBAOcov(paramsBAO,path_to_data2_bao);}//get BAOpostcov
else{paramsBAO[0]=-1;paramsBAO[1]=-1;paramsBAO[2]=-1;}

    Nrealizations_used_SGC=get_cov_from_mocks(path_to_mocks2_bao,path_to_mocks2_rsd,path_to_mocks2_bis_bao,path_to_mocks2_bis_rsd,covfile2,covSGC,Ncov,Nrealizations, NeffP0baoSGC,NeffP0rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, NeffB0baoSGC,NeffB0rsdSGC,errP0baoSGC,errP0rsdSGC,errP2baoSGC,errP2rsdSGC,errP4baoSGC,errP4rsdSGC,errB0baoSGC,errB0rsdSGC,kminP0bao,kminP0rsd,kmaxP0bao,kmaxP0rsd,kminP2bao,kminP2rsd,kmaxP2bao,kmaxP2rsd,kminP4bao,kminP4rsd,kmaxP4bao,kmaxP4rsd,kminB0bao,kminB0rsd,kmaxB0bao,kmaxB0rsd,type_of_analysis,fit_BAO,fit_RSD,do_power_spectrum, do_bispectrum,bispectrum_BQ,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction,covarianceFSa_option,paramsBAO);
  }
} 
else{
  sprintf(covfile2,"%s/%s_2",path_output,identifier);
  printf("Covfile = %s\n",covfile2);

if( strcmp(type_of_analysis ,"FSalphasrecon") == 0 && strcmp(covarianceFSa_option,"varying") == 0 ){getBAOcov(paramsBAO,path_to_data2_bao);}//get BAOpostcov
else{paramsBAO[0]=-1;paramsBAO[1]=-1;paramsBAO[2]=-1;}

  Nrealizations_used_SGC=get_cov_from_mocks(path_to_mocks2_bao,path_to_mocks2_rsd,path_to_mocks2_bis_bao,path_to_mocks2_bis_rsd,covfile2,covSGC,Ncov,Nrealizations, NeffP0baoSGC,NeffP0rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, NeffB0baoSGC,NeffB0rsdSGC,errP0baoSGC,errP0rsdSGC,errP2baoSGC,errP2rsdSGC,errP4baoSGC,errP4rsdSGC,errB0baoSGC,errB0rsdSGC,kminP0bao,kminP0rsd,kmaxP0bao,kmaxP0rsd,kminP2bao,kminP2rsd,kmaxP2bao,kmaxP2rsd,kminP4bao,kminP4rsd,kmaxP4bao,kmaxP4rsd,kminB0bao,kminB0rsd,kmaxB0bao,kmaxB0rsd,type_of_analysis,fit_BAO,fit_RSD,do_power_spectrum, do_bispectrum,bispectrum_BQ,matrix_linear_compression,input_elements_data_vector,output_elements_data_vector,covariance_correction,covarianceFSa_option,paramsBAO);

}


}//Nchunks=2


if( strcmp(do_compression,"linear") == 0)
{//Free matrix covariance
freeTokens(matrix_linear_compression,output_elements_data_vector);
}

//Read theory

if(strcmp(type_of_analysis, "FS") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 || strcmp(type_of_analysis, "FSalphasrecon") == 0){

//Read 2-loop stuff (re-write 2loop program!)

                 Ntheory=countlines(perturbation_theory_file);
                 N_inputs=42;
                 //Theory[Ntheory][N_inputs];
                 Theory = (double **) calloc(Ntheory, sizeof(double*));
                 for(i=0;i<Ntheory;i++){Theory[i]= (double *) calloc(N_inputs, sizeof(double));}
                 k_in_mask = (double *) calloc(N_inputs, sizeof(double));

                 get_Pk_theory(perturbation_theory_file,Theory,Ntheory,N_inputs);
                 for(i=0;i<N_inputs;i++){k_in_mask[i]=Theory[0][i];}

determine_spacing_theo2(spacing_theory_rsd,Theory,Ntheory);


               if(strcmp(do_power_spectrum, "yes") == 0){

  if( strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0)
      {
           if(Theory[0][0]>k0rsd[0] || Theory[Ntheory-1][0]<k0rsd[NeffP0rsd-1]){ printf("Error, power spectrum P0 NGC data out of range of Ptheory file. Exiting now...\n");exit(0);}
           if(Nchunks==2){
          if(Theory[0][0]>k0rsdSGC[0] || Theory[Ntheory-1][0]<k0rsdSGC[NeffP0rsdSGC-1]){ printf("Error, power spectrum P0 SGC data out of range of Ptheory file. Exiting now...\n");exit(0);}
                         }
      }

      if( strcmp(fit_RSD, "P2") == 0 ||  strcmp(fit_RSD, "P02") == 0 ||  strcmp(fit_RSD, "P24") == 0 ||  strcmp(fit_RSD, "P024") == 0)
      {
        if(Theory[0][0]>k2rsd[0] || Theory[Ntheory-1][0]<k2rsd[NeffP2rsd-1]){ printf("Error, power spectrum P2 NGC data out of range of Ptheory file: %lf>%lf, %lf>%lf. Exiting now...\n",Theory[0][0],k2rsd[0],Theory[Ntheory-1][0],k2rsd[NeffP2rsd-1]);exit(0);}      
        if(Nchunks==2){
        if(Theory[0][0]>k2rsdSGC[0] || Theory[Ntheory-1][0]<k2rsdSGC[NeffP2rsdSGC-1]){ printf("Error, power spectrum P2 SGC data out of range of Ptheory file. Exiting now...\n");exit(0);}
                       }

      }
      if( strcmp(fit_RSD, "P4") == 0 ||  strcmp(fit_RSD, "P04") == 0 ||  strcmp(fit_RSD, "P24") == 0 ||  strcmp(fit_RSD, "P024") == 0)
      {
        if(Theory[0][0]>k4rsd[0] || Theory[Ntheory-1][0]<k4rsd[NeffP4rsd-1]){ printf("Error, power spectrum P4 NGC data out of range of Ptheory file. Exiting now...\n");exit(0);}   
        if(Nchunks==2)
        {
        if(Theory[0][0]>k4rsdSGC[0] || Theory[Ntheory-1][0]<k4rsdSGC[NeffP4rsdSGC-1]){ printf("Error, power spectrum P4 SGC data out of range of Ptheory file. Exiting now...\n");exit(0);}
        }

      }


               }

                if(strcmp(do_bispectrum, "yes") == 0){

           if(Theory[0][0]>k11rsd[0] || Theory[0][0]>k22rsd[0] || Theory[0][0]>k33rsd[0] || Theory[Ntheory-1][0]<k11rsd[NeffB0rsd-1] || Theory[Ntheory-1][0]<k22rsd[NeffB0rsd-1] || Theory[Ntheory-1][0]<k33rsd[NeffB0rsd-1]){ printf("Error, bispectrum data out of range of Ptheory file. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(Theory[0][0]>k11rsdSGC[0] || Theory[0][0]>k22rsdSGC[0] || Theory[0][0]>k33rsdSGC[0] || Theory[Ntheory-1][0]<k11rsdSGC[NeffB0rsdSGC-1] || Theory[Ntheory-1][0]<k22rsdSGC[NeffB0rsdSGC-1] || Theory[Ntheory-1][0]<k33rsdSGC[NeffB0rsdSGC-1]){ printf("Error, bispectrum data out of range of Ptheory file. Exiting now...\n");exit(0);}
        }

}


}
//printf("%d %d\n",bispectrum_needs,olin_exists);
N_Plin=0;
k_Plin=NULL;
Plin=NULL;
if(strcmp(type_of_analysis, "BAOISO") == 0 || strcmp(type_of_analysis, "BAOANISO") == 0 || strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0 || bispectrum_needs==1){
//add do_bispectrum with GM14 method


                 N_Plin=countlines(path_to_Plin);
                 k_Plin= (double*) calloc( N_Plin, sizeof(double));
                 Plin= (double*) calloc( N_Plin, sizeof(double));
N_Olin=0;
k_Olin=NULL;
Olin=NULL;
if( olin_exists==1){//file N_Olin exists
                 N_Olin=countlines(path_to_Olin);
                 k_Olin= (double*) calloc( N_Olin, sizeof(double));
                 Olin= (double*) calloc( N_Olin, sizeof(double));

if(N_Plin != N_Olin){printf("Warning, files %s and %s should have the same k-sampling. Exiting now...\n",path_to_Plin,path_to_Olin);exit(0);}
}

get_Pk_bao(path_to_Plin,k_Plin, Plin);
determine_spacing_theo(spacing_theory,k_Plin,N_Plin,0);

if(olin_exists==1){
get_Pk_bao(path_to_Olin,k_Olin, Olin);
determine_spacing_theo(spacing_theory2,k_Olin,N_Olin,1);

if( strcmp(spacing_theory,spacing_theory2) !=0){printf("Warning, theory in Plin and Olin, spaced differently: Plin is spaced %s, whereas Olin is spaced %s. Exiting now...\n",spacing_theory,spacing_theory2);exit(0);}
}
//printf("%lf %lf\n",k_Plin[0],k_Plin[N_Plin-1]);
if(strcmp(do_power_spectrum, "yes") == 0 &&  olin_exists==1){

      if( strcmp(fit_BAO, "P024") == 0 && strcmp(type_of_analysis_BAO,"yes") == 0)
      {
           if(k_Plin[0]>k0bao[0] || k_Plin[0]>k2bao[0] || k_Plin[0]>k4bao[0] || k_Olin[0]>k0bao[0] || k_Olin[0]>k2bao[0] || k_Olin[0]>k4bao[0] || k_Plin[N_Plin-1]<k0bao[NeffP0bao-1] || k_Plin[N_Plin-1]<k2bao[NeffP2bao-1] || k_Plin[N_Plin-1]<k4bao[NeffP4bao-1] || k_Olin[N_Olin-1]<k0bao[NeffP0bao-1] || k_Olin[N_Olin-1]<k2bao[NeffP2bao-1] || k_Olin[N_Olin-1]<k4bao[NeffP4bao-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k0baoSGC[0] || k_Plin[0]>k2baoSGC[0] || k_Plin[0]>k4baoSGC[0] || k_Olin[0]>k0baoSGC[0] || k_Olin[0]>k2baoSGC[0] || k_Olin[0]>k4baoSGC[0] || k_Plin[N_Plin-1]<k0baoSGC[NeffP0baoSGC-1] || k_Plin[N_Plin-1]<k2baoSGC[NeffP2baoSGC-1] || k_Plin[N_Plin-1]<k4baoSGC[NeffP4baoSGC-1] || k_Olin[N_Olin-1]<k0baoSGC[NeffP0baoSGC-1] || k_Olin[N_Olin-1]<k2baoSGC[NeffP2baoSGC-1] || k_Olin[N_Olin-1]<k4baoSGC[NeffP4baoSGC-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        }

      } 

      if( strcmp(fit_BAO, "P02") == 0 && strcmp(type_of_analysis_BAO,"yes") == 0)
      {
           if(k_Plin[0]>k0bao[0] || k_Plin[0]>k2bao[0]  || k_Olin[0]>k0bao[0] || k_Olin[0]>k2bao[0]  || k_Plin[N_Plin-1]<k0bao[NeffP0bao-1] || k_Plin[N_Plin-1]<k2bao[NeffP2bao-1] || k_Olin[N_Olin-1]<k0bao[NeffP0bao-1] || k_Olin[N_Olin-1]<k2bao[NeffP2bao-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k0baoSGC[0] || k_Plin[0]>k2baoSGC[0]  || k_Olin[0]>k0baoSGC[0] || k_Olin[0]>k2baoSGC[0]  || k_Plin[N_Plin-1]<k0baoSGC[NeffP0baoSGC-1] || k_Plin[N_Plin-1]<k2baoSGC[NeffP2baoSGC-1] || k_Olin[N_Olin-1]<k0baoSGC[NeffP0baoSGC-1] || k_Olin[N_Olin-1]<k2baoSGC[NeffP2baoSGC-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        }


      }

      if( strcmp(fit_BAO, "P04") == 0 && strcmp(type_of_analysis_BAO,"yes") == 0)
      {
           if(k_Plin[0]>k0bao[0] || k_Plin[0]>k4bao[0] || k_Olin[0]>k0bao[0] || k_Olin[0]>k4bao[0] || k_Plin[N_Plin-1]<k0bao[NeffP0bao-1] || k_Plin[N_Plin-1]<k4bao[NeffP4bao-1] || k_Olin[N_Olin-1]<k0bao[NeffP0bao-1] || k_Olin[N_Olin-1]<k4bao[NeffP4bao-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k0baoSGC[0] || k_Plin[0]>k4baoSGC[0] || k_Olin[0]>k0baoSGC[0] || k_Olin[0]>k4baoSGC[0] || k_Plin[N_Plin-1]<k0baoSGC[NeffP0baoSGC-1] || k_Plin[N_Plin-1]<k4baoSGC[NeffP4baoSGC-1] || k_Olin[N_Olin-1]<k0baoSGC[NeffP0baoSGC-1] || k_Olin[N_Olin-1]<k4baoSGC[NeffP4baoSGC-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        }


      }

      if( strcmp(fit_BAO, "P24") == 0 && strcmp(type_of_analysis_BAO,"yes") == 0)
      {
           if(k_Plin[0]>k2bao[0] || k_Plin[0]>k4bao[0] ||  k_Olin[0]>k2bao[0] || k_Olin[0]>k4bao[0] || k_Plin[N_Plin-1]<k2bao[NeffP2bao-1] || k_Plin[N_Plin-1]<k4bao[NeffP4bao-1] || k_Olin[N_Olin-1]<k2bao[NeffP2bao-1] || k_Olin[N_Olin-1]<k4bao[NeffP4bao-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k2baoSGC[0] || k_Plin[0]>k4baoSGC[0] ||  k_Olin[0]>k2baoSGC[0] || k_Olin[0]>k4baoSGC[0] || k_Plin[N_Plin-1]<k2baoSGC[NeffP2baoSGC-1] || k_Plin[N_Plin-1]<k4baoSGC[NeffP4baoSGC-1] || k_Olin[N_Olin-1]<k2baoSGC[NeffP2baoSGC-1] || k_Olin[N_Olin-1]<k4baoSGC[NeffP4baoSGC-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        }


      }

      if( strcmp(fit_BAO, "P0") == 0 && strcmp(type_of_analysis_BAO,"yes") == 0)
      {
           if(k_Plin[0]>k0bao[0] || k_Olin[0]>k0bao[0] || k_Plin[N_Plin-1]<k0bao[NeffP0bao-1] || k_Olin[N_Olin-1]<k0bao[NeffP0bao-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k0baoSGC[0] || k_Olin[0]>k0baoSGC[0] || k_Plin[N_Plin-1]<k0baoSGC[NeffP0baoSGC-1] || k_Olin[N_Olin-1]<k0baoSGC[NeffP0baoSGC-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        }


      }

      if( strcmp(fit_BAO, "P2") == 0 && strcmp(type_of_analysis_BAO,"yes") == 0)
      {
           if(k_Plin[0]>k2bao[0] || k_Olin[0]>k2bao[0] || k_Plin[N_Plin-1]<k2bao[NeffP2bao-1] || k_Olin[N_Olin-1]<k2bao[NeffP2bao-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k2baoSGC[0] || k_Olin[0]>k2baoSGC[0] || k_Plin[N_Plin-1]<k2baoSGC[NeffP2baoSGC-1] || k_Olin[N_Olin-1]<k2baoSGC[NeffP2baoSGC-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        }


      }

      if( strcmp(fit_BAO, "P4") == 0 && strcmp(type_of_analysis_BAO,"yes") == 0)
      {
           if(k_Plin[0]>k4bao[0] || k_Olin[0]>k4bao[0] || k_Plin[N_Plin-1]<k4bao[NeffP4bao-1] || k_Olin[N_Olin-1]<k4bao[NeffP4bao-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k4baoSGC[0] || k_Olin[0]>k4baoSGC[0] || k_Plin[N_Plin-1]<k4baoSGC[NeffP4baoSGC-1] || k_Olin[N_Olin-1]<k4baoSGC[NeffP4baoSGC-1]){ printf("Error, power spectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        }


      }



}


if(strcmp(do_bispectrum, "yes") == 0 &&  olin_exists==1 && strcmp(type_of_analysis_BAO,"yes") == 0){

           if(k_Plin[0]>k11bao[0] || k_Plin[0]>k22bao[0] || k_Plin[0]>k33bao[0] || k_Olin[0]>k11bao[0] || k_Olin[0]>k22bao[0] || k_Olin[0]>k33bao[0] || k_Plin[N_Plin-1]<k11bao[NeffB0bao-1] || k_Plin[N_Plin-1]<k22bao[NeffB0bao-1] || k_Plin[N_Plin-1]<k33bao[NeffB0bao-1] || k_Olin[N_Olin-1]<k11bao[NeffB0bao-1] || k_Olin[N_Olin-1]<k22bao[NeffB0bao-1] || k_Olin[N_Olin-1]<k33bao[NeffB0bao-1]){ printf("Error, bispectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k11baoSGC[0] || k_Plin[0]>k22baoSGC[0] || k_Plin[0]>k33baoSGC[0] || k_Olin[0]>k11baoSGC[0] || k_Olin[0]>k22baoSGC[0] || k_Olin[0]>k33baoSGC[0] || k_Plin[N_Plin-1]<k11baoSGC[NeffB0baoSGC-1] || k_Plin[N_Plin-1]<k22baoSGC[NeffB0baoSGC-1] || k_Plin[N_Plin-1]<k33baoSGC[NeffB0baoSGC-1] || k_Olin[N_Olin-1]<k11baoSGC[NeffB0baoSGC-1] || k_Olin[N_Olin-1]<k22baoSGC[NeffB0baoSGC-1] || k_Olin[N_Olin-1]<k33baoSGC[NeffB0baoSGC-1]){ printf("Error, bispectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}
        }

}

if(strcmp(do_bispectrum, "yes") == 0 &&  olin_exists==1 && strcmp(type_of_analysis_FS,"yes") == 0){

           if(k_Plin[0]>k11rsd[0] || k_Plin[0]>k22rsd[0] || k_Plin[0]>k33rsd[0] || k_Olin[0]>k11rsd[0] || k_Olin[0]>k22rsd[0] || k_Olin[0]>k33rsd[0] || k_Plin[N_Plin-1]<k11rsd[NeffB0rsd-1] || k_Plin[N_Plin-1]<k22rsd[NeffB0rsd-1] || k_Plin[N_Plin-1]<k33rsd[NeffB0rsd-1] || k_Olin[N_Olin-1]<k11rsd[NeffB0rsd-1] || k_Olin[N_Olin-1]<k22rsd[NeffB0rsd-1] || k_Olin[N_Olin-1]<k33rsd[NeffB0rsd-1]){ printf("Error, bispectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k11rsdSGC[0] || k_Plin[0]>k22rsdSGC[0] || k_Plin[0]>k33rsdSGC[0] || k_Olin[0]>k11rsdSGC[0] || k_Olin[0]>k22rsdSGC[0] || k_Olin[0]>k33rsdSGC[0] || k_Plin[N_Plin-1]<k11rsdSGC[NeffB0rsdSGC-1] || k_Plin[N_Plin-1]<k22rsdSGC[NeffB0rsdSGC-1] || k_Plin[N_Plin-1]<k33rsdSGC[NeffB0rsdSGC-1] || k_Olin[N_Olin-1]<k11rsdSGC[NeffB0rsdSGC-1] || k_Olin[N_Olin-1]<k22rsdSGC[NeffB0rsdSGC-1] || k_Olin[N_Olin-1]<k33rsdSGC[NeffB0rsdSGC-1]){ printf("Error, bispectrum data out of range of Plin or Olin ranges. Exiting now...\n");exit(0);}
        }

}


if(strcmp(do_bispectrum, "yes") == 0 && olin_exists==0 && strcmp(type_of_analysis_FS,"yes") == 0){

           if(k_Plin[0]>k11rsd[0] || k_Plin[0]>k22rsd[0] || k_Plin[0]>k33rsd[0] || k_Plin[N_Plin-1]<k11rsd[NeffB0rsd-1] || k_Plin[N_Plin-1]<k22rsd[NeffB0rsd-1] || k_Plin[N_Plin-1]<k33rsd[NeffB0rsd-1]){ printf("Error, bispectrum data out of range of Plin ranges. Exiting now...\n");exit(0);}

        if(Nchunks==2)
        {
           if(k_Plin[0]>k11rsdSGC[0] || k_Plin[0]>k22rsdSGC[0] || k_Plin[0]>k33rsdSGC[0] || k_Plin[N_Plin-1]<k11rsdSGC[NeffB0rsdSGC-1] || k_Plin[N_Plin-1]<k22rsdSGC[NeffB0rsdSGC-1] || k_Plin[N_Plin-1]<k33rsdSGC[NeffB0rsdSGC-1]){ printf("Error, bispectrum data out of range of Plin ranges. Exiting now...\n");exit(0);}
        }

}

} 


//mask as a matrix (for P-RSD and P-BAOANISO only)
//check what happens if bispectrum = yes


MatrixFS_mask_NGC=NULL;N_matrixFS_in=0;
MatrixFS_mask_SGC=NULL;N_matrixFS_in_SGC=0;

MatrixBAO_mask_NGC=NULL;N_matrixBAO_in=0;
MatrixBAO_mask_SGC=NULL;N_matrixBAO_in_SGC=0;

if( strcmp(path_to_mask1, "none") != 0 && strcmp(mask_matrix,"yes") == 0 )
{

if( strcmp(do_power_spectrum,"yes") == 0 )
{
printf("Generating mask-matrix for NGC...\n");
if( strcmp(type_of_analysis_FS ,"yes") == 0)//imply FS, FSBAOISO, FSBAOANISO, FSalphasrecon
{

get_kmin_kmax(params_ks,modeP0rsd,modeP2rsd,modeP4rsd,k0rsd[0],k0rsd[NeffP0rsd-1],k2rsd[0],k2rsd[NeffP2rsd-1],k4rsd[0],k4rsd[NeffP4rsd-1]);
kmax=params_ks[1];
kmin=params_ks[0];

get_NeffP(params_int,NeffP0rsd,NeffP2rsd,NeffP4rsd,spacing_dataNGC_rsd,kmin,kmax,k0rsd[0],k0rsd[NeffP0rsd-1],k2rsd[0],k2rsd[NeffP2rsd-1],k4rsd[0],k4rsd[NeffP4rsd-1], Ntheory );
if(strcmp(spacing_dataNGC_rsd,"irregular") == 0){N_matrixFS_in=params_int[0];}
else{N_matrixFS_in=params_int[0]*factor_for_sampling;}

N_matrixFS_in=N_matrixFS_in*(modeP0rsd+modeP2rsd+modeP4rsd);

MatrixFS_mask_NGC = (double **) calloc(N_matrixFS_in, sizeof(double*));
                 for(i=0;i<N_matrixFS_in;i++){MatrixFS_mask_NGC[i]= (double *) calloc(NcovP, sizeof(double));}

do_window_matrix_for_P(MatrixFS_mask_NGC, type_of_analysis, modeP0rsd, modeP2rsd, modeP4rsd, k_in_mask, Ntheory, NeffP0rsd, NeffP2rsd, NeffP4rsd, factor_for_sampling, pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC, k0rsd, k2rsd, k4rsd, k0rsd[0], k0rsd[NeffP0rsd-1], k2rsd[0], k2rsd[NeffP2rsd-1], k4rsd[0], k4rsd[NeffP4rsd-1],  spacing_dataNGC_rsd, NcovP, path_output, identifier, 1, 0);


}


if( strcmp(type_of_analysis_BAO ,"yes") == 0)
{

get_kmin_kmax(params_ks,modeP0bao,modeP2bao,modeP4bao,k0bao[0],k0bao[NeffP0bao-1],k2bao[0],k2bao[NeffP2bao-1],k4bao[0],k4bao[NeffP4bao-1]);
kmax=params_ks[1];
kmin=params_ks[0];

get_NeffP(params_int,NeffP0bao,NeffP2bao,NeffP4bao,spacing_dataNGC_bao,kmin,kmax,k0bao[0],k0bao[NeffP0bao-1],k2bao[0],k2bao[NeffP2bao-1],k4bao[0],k4bao[NeffP4bao-1], N_Plin );
if(strcmp(spacing_dataNGC_bao,"irregular") == 0){N_matrixBAO_in=params_int[0];}
else{N_matrixBAO_in=params_int[0]*factor_for_sampling;}

N_matrixBAO_in=N_matrixBAO_in*(modeP0bao+modeP2bao+modeP4bao);

MatrixBAO_mask_NGC = (double **) calloc(N_matrixBAO_in, sizeof(double*));
                 for(i=0;i<N_matrixBAO_in;i++){MatrixBAO_mask_NGC[i]= (double *) calloc(NcovP, sizeof(double));}

do_window_matrix_for_P(MatrixBAO_mask_NGC, type_of_analysis, modeP0bao, modeP2bao, modeP4bao, k_Plin, N_Plin, NeffP0bao, NeffP2bao, NeffP4bao, factor_for_sampling, pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC, k0bao, k2bao, k4bao, k0bao[0], k0bao[NeffP0bao-1], k2bao[0], k2bao[NeffP2bao-1], k4bao[0], k4bao[NeffP4bao-1],  spacing_dataNGC_bao, NcovP, path_output, identifier, 1, 1);



if( strcmp(type_of_analysis_BAO ,"yes") == 0 && strcmp(type_of_analysis ,"BAOISO") == 0 &&  strcmp(type_BAO_fit, "analytic") == 0)//check if needed
{
printf("Warning No mask-matrix implemented for BAOISO-analytic analysis. Exiting now ...\n");exit(0);
}

}

printf("done\n");
}//power spectrum yes
}//mask_ngc = yes

if(Nchunks==2){

if( strcmp(path_to_mask1, "none") != 0 && strcmp(mask_matrix,"yes") == 0 )
{

if( strcmp(do_power_spectrum,"yes") == 0 )
{

printf("Generating mask-matrix for SGC...\n");

if( strcmp(type_of_analysis_FS ,"yes") == 0)//imply FS, FSBAOISO, FSBAOANISO, FSalphasrecon
{


get_kmin_kmax(params_ks,modeP0rsd,modeP2rsd,modeP4rsd,k0rsdSGC[0],k0rsdSGC[NeffP0rsdSGC-1],k2rsdSGC[0],k2rsdSGC[NeffP2rsdSGC-1],k4rsdSGC[0],k4rsdSGC[NeffP4rsdSGC-1]);
kmax=params_ks[1];
kmin=params_ks[0];

get_NeffP(params_int,NeffP0rsdSGC,NeffP2rsdSGC,NeffP4rsdSGC,spacing_dataSGC_rsd,kmin,kmax,k0rsdSGC[0],k0rsdSGC[NeffP0rsdSGC-1],k2rsdSGC[0],k2rsdSGC[NeffP2rsdSGC-1],k4rsdSGC[0],k4rsdSGC[NeffP4rsdSGC-1], Ntheory );
if(strcmp(spacing_dataSGC_rsd,"irregular") == 0){N_matrixFS_in_SGC=params_int[0];}
else{N_matrixFS_in_SGC=params_int[0]*factor_for_sampling;}

N_matrixFS_in_SGC=N_matrixFS_in_SGC*(modeP0rsd+modeP2rsd+modeP4rsd);

MatrixFS_mask_SGC = (double **) calloc(N_matrixFS_in_SGC, sizeof(double*));
                 for(i=0;i<N_matrixFS_in_SGC;i++){MatrixFS_mask_SGC[i]= (double *) calloc(NcovP, sizeof(double));}

do_window_matrix_for_P(MatrixFS_mask_SGC, type_of_analysis, modeP0rsd, modeP2rsd, modeP4rsd, k_in_mask, Ntheory, NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, factor_for_sampling, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC, k0rsdSGC, k2rsdSGC, k4rsdSGC, k0rsdSGC[0], k0rsdSGC[NeffP0rsdSGC-1], k2rsdSGC[0], k2rsdSGC[NeffP2rsdSGC-1], k4rsdSGC[0], k4rsdSGC[NeffP4rsdSGC-1],  spacing_dataSGC_rsd, NcovP, path_output, identifier, 2, 0);


}

if( strcmp(type_of_analysis_BAO ,"yes") == 0)
{

get_kmin_kmax(params_ks,modeP0bao,modeP2bao,modeP4bao,k0baoSGC[0],k0baoSGC[NeffP0baoSGC-1],k2baoSGC[0],k2baoSGC[NeffP2baoSGC-1],k4baoSGC[0],k4baoSGC[NeffP4baoSGC-1]);
kmax=params_ks[1];
kmin=params_ks[0];

get_NeffP(params_int,NeffP0baoSGC,NeffP2baoSGC,NeffP4baoSGC,spacing_dataSGC_bao,kmin,kmax,k0baoSGC[0],k0baoSGC[NeffP0baoSGC-1],k2baoSGC[0],k2baoSGC[NeffP2baoSGC-1],k4baoSGC[0],k4baoSGC[NeffP4baoSGC-1], N_Plin );
if(strcmp(spacing_dataSGC_bao,"irregular") == 0){N_matrixBAO_in_SGC=params_int[0];}
else{N_matrixBAO_in_SGC=params_int[0]*factor_for_sampling;}

N_matrixBAO_in_SGC=N_matrixBAO_in_SGC*(modeP0bao+modeP2bao+modeP4bao);

MatrixBAO_mask_SGC = (double **) calloc(N_matrixBAO_in_SGC, sizeof(double*));
                 for(i=0;i<N_matrixBAO_in_SGC;i++){MatrixBAO_mask_SGC[i]= (double *) calloc(NcovP, sizeof(double));}

do_window_matrix_for_P(MatrixBAO_mask_SGC, type_of_analysis, modeP0bao, modeP2bao, modeP4bao, k_Plin, N_Plin, NeffP0baoSGC, NeffP2baoSGC, NeffP4baoSGC, factor_for_sampling, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC, k0baoSGC, k2baoSGC, k4baoSGC, k0baoSGC[0], k0baoSGC[NeffP0baoSGC-1], k2baoSGC[0], k2baoSGC[NeffP2baoSGC-1], k4baoSGC[0], k4baoSGC[NeffP4baoSGC-1],  spacing_dataSGC_bao, NcovP, path_output, identifier, 2, 1);



if( strcmp(type_of_analysis_BAO ,"yes") == 0 && strcmp(type_of_analysis ,"BAOISO") == 0 &&  strcmp(type_BAO_fit, "analytic") == 0)//check if needed
{
printf("Warning No mask-matrix implemented for BAOISO-analytic analysis. Exiting now ...\n");exit(0);
}

}

printf("done\n");
}//power spectrum yes
}//mask_sgc = yes 


}//Nchucks=2



//Do BAO & RSD


if( strcmp(type_of_analysis_FS ,"yes") == 0 && strcmp(type_of_analysis_BAO ,"no") == 0){

 do_rsd_mcmc(nthreads,type_of_analysis,fit_RSD,Theory,Ntheory,k_Plin,Plin,N_Plin,k_Olin, Olin, N_Olin, pos, W0, W2, W4,W6, W8,Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0rsd, P0rsd, errP0rsd,Pnoise_rsd, NeffP0rsd, k2rsd, P2rsd, errP2rsd, NeffP2rsd, k4rsd, P4rsd, errP4rsd, NeffP4rsd, k11rsd, k22rsd, k33rsd, B0rsd, errB0rsd, Bnoise_rsd, NeffB0rsd, k0rsdSGC, P0rsdSGC, errP0rsdSGC,Pnoise_rsdSGC,NeffP0rsdSGC, k2rsdSGC, P2rsdSGC, errP2rsdSGC,NeffP2rsdSGC, k4rsdSGC, P4rsdSGC, errP4rsdSGC,NeffP4rsdSGC, k11rsdSGC, k22rsdSGC, k33rsdSGC,B0rsdSGC,errB0rsdSGC, Bnoise_rsdSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max,ptmodel_ps,rsdmodel_ps,fogmodel_ps,ptmodel_bs,local_b2s2,local_b3nl,RSD_fit,sigma8_free,fog_free,fog_bs ,Nchunks, path_output, identifier, do_plot, use_prop_cov, path_to_cov, Nsteps, do_power_spectrum, do_bispectrum,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS,bispectrum_BQ, factor_for_sampling, mask_matrix,MatrixFS_mask_NGC, MatrixFS_mask_SGC,FSprior_type,FSprior_mean,FSprior_stddev,noise_option,step_size,covariance_correction, Nrealizations_used_NGC,Nrealizations_used_SGC,P0bao,P0baoSGC );

}

if(strcmp(type_of_analysis_BAO ,"yes") == 0 && strcmp(type_of_analysis_FS ,"no") == 0){

//analytic way
if(strcmp(type_BAO_fit, "analytic") == 0 && strcmp(type_of_analysis, "BAOISO") == 0){

//bispectrum not available here
do_bao_analytic(type_BAO_fit,type_of_analysis,fit_BAO,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, pos, W0,W2,W4,W6,W8,Nmask,path_to_mask1, spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao, P0bao, errP0bao,NeffP0bao, k2bao, P2bao, errP2bao,NeffP2bao, k4bao, P4bao, errP4bao,NeffP4bao,k0baoSGC,P0baoSGC,errP0baoSGC, NeffP0baoSGC,k2baoSGC,P2baoSGC,errP2baoSGC, NeffP2baoSGC,k4baoSGC,P4baoSGC,errP4baoSGC, NeffP4baoSGC, cov, covSGC, alpha_min,alpha_max,alpha_step, Sigma_def_type, Sigma_independent, ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev,    Npolynomial, Nchunks, path_output,identifier, do_plot, do_power_spectrum , do_bispectrum,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory,factor_for_sampling, covariance_correction, Nrealizations_used_NGC,Nrealizations_used_SGC);
}

//mcmc way
if(strcmp(type_BAO_fit, "mcmc") == 0 && strcmp(type_of_analysis, "BAOISO") == 0){

//here it goes the bispectrum as BAOISO

 do_bao_mcmc(nthreads,type_BAO_fit,type_of_analysis,fit_BAO,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4,W6, W8,Nmask, path_to_mask1, spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao, P0bao, errP0bao, NeffP0bao, k2bao, P2bao, errP2bao, NeffP2bao, k4bao, P4bao, errP4bao, NeffP4bao, k11bao, k22bao, k33bao, B0bao, errB0bao, Bnoise_bao, NeffB0bao, k0baoSGC, P0baoSGC, errP0baoSGC,NeffP0baoSGC, k2baoSGC, P2baoSGC, errP2baoSGC,NeffP2baoSGC, k4baoSGC, P4baoSGC, errP4baoSGC,NeffP4baoSGC, k11baoSGC, k22baoSGC, k33baoSGC,B0baoSGC,errB0baoSGC, Bnoise_baoSGC,NeffB0baoSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent, ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev,  Npolynomial,  Nchunks, path_output, identifier, do_plot, use_prop_cov, path_to_cov, Nsteps, do_power_spectrum, do_bispectrum,0,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS,bispectrum_BQ,factor_for_sampling, mask_matrix,MatrixBAO_mask_NGC, MatrixBAO_mask_SGC,step_size,covariance_correction, Nrealizations_used_NGC,Nrealizations_used_SGC);
}

if(strcmp(type_BAO_fit, "mcmc") == 0 && strcmp(type_of_analysis, "BAOANISO") == 0){

//we allow as well to do bispectrum under baoaniso (for P only)
 do_bao_mcmc(nthreads,type_BAO_fit,type_of_analysis,fit_BAO,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4,W6, W8,Nmask, path_to_mask1, spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao, P0bao, errP0bao, NeffP0bao, k2bao, P2bao, errP2bao, NeffP2bao, k4bao, P4bao, errP4bao, NeffP4bao, k11bao, k22bao, k33bao, B0bao, errB0bao, Bnoise_bao, NeffB0bao, k0baoSGC, P0baoSGC, errP0baoSGC,NeffP0baoSGC, k2baoSGC, P2baoSGC, errP2baoSGC,NeffP2baoSGC, k4baoSGC, P4baoSGC, errP4baoSGC,NeffP4baoSGC, k11baoSGC, k22baoSGC, k33baoSGC,B0baoSGC,errB0baoSGC, Bnoise_baoSGC,NeffB0baoSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent, ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev,  Npolynomial,  Nchunks, path_output, identifier, do_plot, use_prop_cov, path_to_cov, Nsteps, do_power_spectrum, do_bispectrum,Sigma_smooth,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS,bispectrum_BQ,factor_for_sampling,mask_matrix,MatrixBAO_mask_NGC, MatrixBAO_mask_SGC,step_size,covariance_correction, Nrealizations_used_NGC,Nrealizations_used_SGC);
}

}

if( strcmp(type_of_analysis_FS ,"yes") == 0 && strcmp(type_of_analysis_BAO ,"yes") == 0){

//to be done, include FS+BAO  for P+B
if(strcmp(do_bispectrum,"yes")==0){printf("FS+BAO not available for bispectrum so far. Exiting...\n");exit(0);}

do_rsd_bao_mcmc(nthreads,type_BAO_fit,type_of_analysis,fit_BAO,fit_RSD,Theory,Ntheory,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4,W6, W8,Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao,k0rsd, P0bao,P0rsd, errP0bao,errP0rsd,Pnoise_rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, errP2bao,errP2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, errP4bao,errP4rsd, NeffP4bao,NeffP4rsd, k11bao,k11rsd, k22bao,k22rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, errP0baoSGC,errP0rsdSGC,Pnoise_rsdSGC,NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, errP2baoSGC,errP2rsdSGC,NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, errP4baoSGC,errP4rsdSGC,NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC,B0baoSGC,B0rsdSGC,errB0baoSGC,errB0rsdSGC, Bnoise_baoSGC,Bnoise_rsdSGC,NeffB0baoSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent, ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev,  Npolynomial, ptmodel_ps,rsdmodel_ps,fogmodel_ps,ptmodel_bs,local_b2s2,local_b3nl,RSD_fit,sigma8_free,fog_free,fog_bs, Nchunks, path_output, identifier, do_plot, use_prop_cov, path_to_cov, Nsteps, do_power_spectrum, do_bispectrum,Sigma_smooth,spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS,bispectrum_BQ,factor_for_sampling,mask_matrix,MatrixBAO_mask_NGC, MatrixBAO_mask_SGC,MatrixFS_mask_NGC, MatrixFS_mask_SGC,FSprior_type,FSprior_mean,FSprior_stddev,noise_option,step_size,covariance_correction, Nrealizations_used_NGC,Nrealizations_used_SGC);

}

if( strcmp(type_of_analysis_FS ,"yes") == 0){
//free Theory
freeTokens(Theory,Ntheory);
free(k_in_mask);
if( strcmp(mask_matrix,"yes") == 0)
{
freeTokens(MatrixFS_mask_NGC,N_matrixFS_in);
if(Nchunks==2){freeTokens(MatrixFS_mask_SGC,N_matrixFS_in_SGC);}
}

}


if (strcmp(type_of_analysis_BAO ,"yes") == 0 && strcmp(mask_matrix,"yes") == 0)
{

if( strcmp(mask_matrix,"yes") == 0)
{
freeTokens(MatrixBAO_mask_NGC,N_matrixBAO_in);
if(Nchunks==2){freeTokens(MatrixBAO_mask_SGC,N_matrixBAO_in_SGC);}
}

}

if (strcmp(type_of_analysis_BAO ,"yes") == 0 || bispectrum_needs==1){
free(k_Plin);
free(Plin);
}

if (olin_exists==1){
free(k_Olin);
free(Olin);
}

//Free
if(strcmp(do_power_spectrum, "yes") == 0){

if( strcmp(type_of_analysis_BAO ,"yes") ==0)
{
free(P0bao);
free(P2bao);
free(P4bao);
free(errP0bao);
free(errP2bao);
free(errP4bao);
free(k0bao);
free(k2bao);
free(k4bao);
free(kav0bao);
free(kav2bao);
free(kav4bao);
}
if( strcmp(type_of_analysis ,"FSalphasrecon") ==0)
{
free(P0bao);
}

if( strcmp(type_of_analysis_FS ,"yes") ==0)
{
free(P0rsd);
free(P2rsd);
free(P4rsd);
free(errP0rsd);
free(errP2rsd);
free(errP4rsd);
free(k0rsd);
free(k2rsd);
free(k4rsd);
free(kav0rsd);
free(kav2rsd);
free(kav4rsd);
}

if(strcmp(path_to_mask1, "none") != 0){
free(pos); 
free(posAV);
free(W0);
free(W2);
free(W4);
free(W6);
free(W8);
}
}
if(strcmp(do_bispectrum, "yes") == 0){

if( strcmp(type_of_analysis_BAO ,"yes") ==0)
{
free(B0bao);
free(errB0bao);
free(k11bao);
free(k22bao);
free(k33bao);
free(Bnoise_bao);
}

if( strcmp(type_of_analysis_FS ,"yes") ==0)
{
free(B0rsd);
free(errB0rsd);
free(k11rsd);
free(k22rsd);
free(k33rsd);
free(Bnoise_rsd);
}

}
free(cov);

if(Nchunks==2)//both chunks considered independent
{

if(strcmp(do_power_spectrum, "yes") == 0){
if( strcmp(type_of_analysis_BAO ,"yes") ==0)
{
free(P0baoSGC);
free(P2baoSGC);
free(P4baoSGC);
free(errP0baoSGC);
free(errP2baoSGC);
free(errP4baoSGC);
free(k0baoSGC);
free(k2baoSGC);
free(k4baoSGC);
free(kav0baoSGC);
free(kav2baoSGC);
free(kav4baoSGC);
}
if( strcmp(type_of_analysis ,"FSalphasrecon") ==0)
{
free(P0baoSGC);
}

if( strcmp(type_of_analysis_FS ,"yes") ==0)
{
free(P0rsdSGC);
free(P2rsdSGC);
free(P4rsdSGC);
free(errP0rsdSGC);
free(errP2rsdSGC);
free(errP4rsdSGC);
free(k0rsdSGC);
free(k2rsdSGC);
free(k4rsdSGC);
free(kav0rsdSGC);
free(kav2rsdSGC);
free(kav4rsdSGC);
}

if(strcmp(path_to_mask2, "none") != 0){
free(posAVSGC);
free(posSGC);
free(W0SGC);
free(W2SGC);
free(W4SGC);
free(W6SGC);
free(W8SGC);
}
}
if(strcmp(do_bispectrum, "yes") == 0){

if( strcmp(type_of_analysis_BAO ,"yes") ==0)
{
free(B0baoSGC);
free(errB0baoSGC);
free(k11baoSGC);
free(k22baoSGC);
free(k33baoSGC);
free(Bnoise_baoSGC);
}

if( strcmp(type_of_analysis_FS ,"yes") ==0)
{
free(B0rsdSGC);
free(errB0rsdSGC);
free(k11rsdSGC);
free(k22rsdSGC);
free(k33rsdSGC);
free(Bnoise_rsdSGC);
}


}
free(covSGC);


}

return 0;
}
