#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#define Pi (4.*atan(1.))


void get_Pk_theory(char *perturbation_theory_file,double **Theory,int N_Plin,int N_inputs)
{
int i,j;
FILE *f;
f=fopen(perturbation_theory_file,"r");
for(i=0;i<N_Plin;i++)
{
     for(j=0;j<N_inputs;j++)
     {
        if(j==N_inputs-1){fscanf(f,"%lf\n",&Theory[i][j]);}
        else{fscanf(f,"%lf ",&Theory[i][j]);}//Theory[N_Plin][N_inputs];
     }

}

fclose(f);
}


int countPk(int mode, char *path ,double kmin,double kmax)
{
int Neff,N,i;
double k;
double k1,k2,k3;
FILE *f;
//N=countlines(path)-27;
N=countlines(path);
f=fopen(path,"r");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*f %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f %*s\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*f\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*s\n");
Neff=0;
for(i=0;i<N;i++)
{

if(mode==0){
fscanf(f,"%*f %lf %*f %*f %*f %*d %*f\n",&k);
if(k>kmin && k<kmax){Neff++;}
}

if(mode==1){
fscanf(f,"%*f %lf %*f %lf %*f %lf %*f %*f %*f %*f %*f\n",&k1,&k2,&k3);
if(k1>kmin && k1<kmax && k2>kmin && k2<kmax && k3>kmin && k3<kmax){Neff++;}
}

}
fclose(f);

return Neff;
}

void get_Pk_bao( char *path, double k[], double P[])
{
int N,i;
FILE *f;
N=countlines(path);
f=fopen(path,"r");
for(i=0;i<N;i++)
{
fscanf(f,"%lf %lf\n",&k[i],&P[i]);
}
fclose(f);
}

void get_data(char *path, double k0[], double kav0[], double P0[], double k2[], double kav2[], double P2[], double k4[], double kav4[], double P4[], double parameter_value[], double kminP0,double kmaxP0, double kminP2, double kmaxP2, double kminP4, double kmaxP4, char *type_BAORSD, int baorsd, int bao)
{
int i,N;
int i0,i2,i4;
FILE *f;
double Pnoise,sumw,I22;
double k,p0,p2,p4,kav;
//N=countlines(path)-27;
N=countlines(path);
f=fopen(path,"r");

fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %lf\n",&sumw);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*f %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f %*s\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %lf\n",&I22);
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*s\n");
printf("I22=%lf, Sumw=%lf\n",I22,sumw);

   i0=0;i2=0;i4=0;
   for(i=0;i<N;i++)
   {
     fscanf(f,"%lf %lf %lf %lf %lf %*d %lf\n",&kav,&k,&p0,&p2,&p4,&Pnoise);
   
        if(k>kminP0 && k<kmaxP0)
        {
           k0[i0]=k;
           kav0[i0]=kav;
           P0[i0]=p0;
           i0++;
        }
        if(strcmp(type_BAORSD, "BAOISO") == 0 || strcmp(type_BAORSD, "FSBAOISO") == 0){

           if(k>kminP2 && k<kmaxP2)
           {
             k2[i2]=k;
             kav2[i2]=kav;
             if(baorsd==0){P2[i2]=p0+2./5.*p2;}
             else{P2[i2]=p2;}
             i2++;
           }
           if(k>kminP4 && k<kmaxP4)
           {
             k4[i4]=k;
             kav4[i4]=kav;
             if(bao==0){P4[i4]=p0+4./7.*p2+8./63.*p4;}
             else{P4[i4]=p4;}
             i4++;
           }
        }
        if(strcmp(type_BAORSD, "FS") == 0 || strcmp(type_BAORSD, "BAOANISO") == 0 ||  strcmp(type_BAORSD, "FSBAOANISO") == 0){

           if(k>kminP2 && k<kmaxP2)
           {
             k2[i2]=k;
             kav2[i2]=kav;
             P2[i2]=p2;
             i2++;
           }
           if(k>kminP4 && k<kmaxP4)
           {
             k4[i4]=k;
             kav4[i4]=kav;
             P4[i4]=p4;
             i4++;
           }
        }


   } 
fclose(f);
parameter_value[0]=Pnoise;
parameter_value[1]=sumw;
parameter_value[2]=I22;
printf("%lf %lf %lf\n",parameter_value[0],parameter_value[1],parameter_value[2]);
//exit(0);
}

void get_data_bis(char *path, double k11[], double k22[], double k33[], double B0[], double Bnoise[], double kminB0, double kmaxB0)
{

int i,N;
int i0;
FILE *f;
double k1,k2,k3,b0,noise;
N=countlines(path);
f=fopen(path,"r");

   i0=0;
   for(i=0;i<N;i++)
   {
     fscanf(f,"%*f %lf %*f %lf %*f %lf %lf %lf %*f %*f %*f\n",&k1,&k2,&k3,&b0,&noise);

        if(k1>kminB0 && k1<kmaxB0 && k2>kminB0 && k2<kmaxB0 && k3>kminB0 && k3<kmaxB0)
        {
           k11[i0]=k1;
           k22[i0]=k2;
           k33[i0]=k3;
           B0[i0]=b0;
           Bnoise[i0]=noise;
           i0++;
        }
   }
   fclose(f);

}

//void get_cov(char *path_to_mocks, char *path_to_mocks_bis, double cov[], int Ncov,int Nrealizations, int NeffP0, int NeffP2, int NeffP4, int NeffB0, double errP0[],double errP2[], double errP4[], double errB0[], double kminP0,double kmaxP0,  double kminP2,double kmaxP2,double kminP4,double kmaxP4,double kminB0,double kmaxB0, char *type_BAORSD, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum)
void get_cov(char *path_to_mocks_bao, char *path_to_mocks_rsd, char *path_to_mocks_bis_bao, char *path_to_mocks_bis_rsd, double cov[], int Ncov,int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao,  int NeffP2rsd, int NeffP4bao, int NeffP4rsd , int NeffB0bao, int NeffB0rsd, double errP0bao[], double errP0rsd[],double errP2bao[], double errP2rsd[], double errP4bao[],  double errP4rsd[] , double errB0bao[], double errB0rsd[], double kminP0bao,  double kminP0rsd,  double kmaxP0bao,  double kmaxP0rsd,  double kminP2bao, double kminP2rsd ,double kmaxP2bao, double kmaxP2rsd,double kminP4bao, double kminP4rsd,double kmaxP4bao, double kmaxP4rsd,double kminB0bao, double kminB0rsd,double kmaxB0bao, double kmaxB0rsd, char *type_BAORSD, char *fit_BAO, char *fit_RSD, char *do_power_spectrum, char *do_bispectrum)
{
int bao,rsd;
int iteration,iteration_ini,iteration_fin;
double scaling_factor;//number of realizations when fitting the mean
int i,j,Nlines,Ncounter,Ncounterbis,i0,i2,i4,i0bis,i_real;
int l0rsd,l0bao,l2rsd,l2bao,l4rsd,l4bao;
int i0rsd,i0bao,i2rsd,i2bao,i4rsd,i4bao;
//int Ncounterbao,Ncounterrsd;
int l,l0,l2,l4;
int s,fail;
double jelement,ielement;
FILE *f1bao,*f2bao,*f3bao,*f4bao;
FILE *f1rsd,*f2rsd,*f3rsd,*f4rsd;
FILE *f,*g;
double kminP0,kminP2,kminP4,kminB0;
double kmaxP0,kmaxP2,kmaxP4,kmaxB0;
int modersd=-1;
int modebao=-1;
double k,k1,k2,k3,p0,p2,p4,b0;
int baoshift,rsdshift;

char full_path_bao[2000];
char full_path2_bao[2000];
char full_path3_bao[2000];
char full_path4_bao[2000];

char full_path_rsd[2000];
char full_path2_rsd[2000];
char full_path3_rsd[2000];
char full_path4_rsd[2000];


double *P0baoav,*P2baoav, *P4baoav, *B0baoav, *Pvectorav;
double **P0bao, **P2bao, **P4bao, **B0bao, **Pvector;
double *P0rsdav,*P2rsdav, *P4rsdav, *B0rsdav;
double **P0rsd, **P2rsd, **P4rsd, **B0rsd;


double *Variance;
//exit(0);
l0bao=0;
l0rsd=0;
l2bao=0;
l2rsd=0;
l4bao=0;
l4rsd=0;
bao=0;rsd=0;baoshift=0;rsdshift=0;
if( strcmp(type_BAORSD,"FSBAOISO") ==0 || strcmp(type_BAORSD,"FSBAOANISO") ==0 ){iteration_ini=1;iteration_fin=2;bao=1;rsd=1;}
if( strcmp(type_BAORSD,"BAOISO") ==0 || strcmp(type_BAORSD,"BAOANISO") ==0 ){iteration_ini=1;iteration_fin=1;bao=1;rsd=0;}
if( strcmp(type_BAORSD,"FS") ==0){iteration_ini=2;iteration_fin=2;rsd=1;bao=0;}
if(bao==1 && rsd==1)
{
if( strcmp(fit_BAO,"P0") ==0 || strcmp(fit_BAO,"P02") ==0 ||  strcmp(fit_BAO,"P04") ==0 || strcmp(fit_BAO,"P024")==0 ){baoshift=baoshift+NeffP0bao;}
if( strcmp(fit_BAO,"P2") ==0 || strcmp(fit_BAO,"P02") ==0 ||  strcmp(fit_BAO,"P24") ==0 || strcmp(fit_BAO,"P024")==0 ){baoshift=baoshift+NeffP2bao;}
if( strcmp(fit_BAO,"P4") ==0 || strcmp(fit_BAO,"P24") ==0 ||  strcmp(fit_BAO,"P24") ==0 || strcmp(fit_BAO,"P024")==0 ){baoshift=baoshift+NeffP4bao;}

if( strcmp(fit_RSD,"P0") ==0 || strcmp(fit_RSD,"P02") ==0 ||  strcmp(fit_RSD,"P04") ==0 || strcmp(fit_RSD,"P024")==0 ){rsdshift=rsdshift+NeffP0rsd;}
if( strcmp(fit_RSD,"P2") ==0 || strcmp(fit_RSD,"P02") ==0 ||  strcmp(fit_RSD,"P24") ==0 || strcmp(fit_RSD,"P024")==0 ){rsdshift=rsdshift+NeffP2rsd;}
if( strcmp(fit_RSD,"P4") ==0 || strcmp(fit_RSD,"P24") ==0 ||  strcmp(fit_RSD,"P24") ==0 || strcmp(fit_RSD,"P024")==0 ){rsdshift=rsdshift+NeffP4rsd;}

}

scaling_factor=1.;

//l=0;l0=2;l2=2;l4=0;

Variance = (double*) calloc( Ncov*Ncov, sizeof(double));

if( strcmp(do_power_spectrum, "yes") == 0)
{

if(bao==1){
P0baoav=  (double*) calloc( NeffP0bao, sizeof(double));
P2baoav=  (double*) calloc( NeffP2bao, sizeof(double));
P4baoav=  (double*) calloc( NeffP4bao, sizeof(double));}

if(rsd==1){
P0rsdav=  (double*) calloc( NeffP0rsd, sizeof(double));
P2rsdav=  (double*) calloc( NeffP2rsd, sizeof(double));
P4rsdav=  (double*) calloc( NeffP4rsd, sizeof(double));}


Pvectorav = (double*) calloc( Ncov, sizeof(double));

if(bao==1){
P0bao = (double**) calloc(Nrealizations,sizeof(double*));
P2bao = (double**) calloc(Nrealizations,sizeof(double*));
P4bao = (double**) calloc(Nrealizations,sizeof(double*));}

if(rsd==1){
P0rsd = (double**) calloc(Nrealizations,sizeof(double*));
P2rsd = (double**) calloc(Nrealizations,sizeof(double*));
P4rsd = (double**) calloc(Nrealizations,sizeof(double*));}
//exit(0);

Pvector = (double**) calloc(Nrealizations,sizeof(double*));

for(j=0;j<Nrealizations;j++)
{

if(bao==1){
   P0bao[j] = (double*) calloc(NeffP0bao,sizeof(double));
   P2bao[j] = (double*) calloc(NeffP2bao,sizeof(double));
   P4bao[j] = (double*) calloc(NeffP4bao,sizeof(double));}

if(rsd==1){
   P0rsd[j] = (double*) calloc(NeffP0rsd,sizeof(double));
   P2rsd[j] = (double*) calloc(NeffP2rsd,sizeof(double));
   P4rsd[j] = (double*) calloc(NeffP4rsd,sizeof(double));}

   Pvector[j] = (double*) calloc(Ncov,sizeof(double));


}
//exit(0);

   Ncounter=0;
   for(j=1;j<=Nrealizations;j++)
   {

if(bao==1){
        sprintf(full_path_bao,"%s%.4d.txt",path_to_mocks_bao,j);
        sprintf(full_path2_bao,"%s%d.txt",path_to_mocks_bao,j);
        sprintf(full_path3_bao,"%s%.4d.dat",path_to_mocks_bao,j);
        sprintf(full_path4_bao,"%s%d.dat",path_to_mocks_bao,j);}

if(rsd==1){
        sprintf(full_path_rsd,"%s%.4d.txt",path_to_mocks_rsd,j);
        sprintf(full_path2_rsd,"%s%d.txt",path_to_mocks_rsd,j);
        sprintf(full_path3_rsd,"%s%.4d.dat",path_to_mocks_rsd,j);
        sprintf(full_path4_rsd,"%s%d.dat",path_to_mocks_rsd,j);}

if(bao==1){
        f1bao=fopen(full_path_bao,"r");
        f2bao=fopen(full_path2_bao,"r");
        f3bao=fopen(full_path3_bao,"r");
        f4bao=fopen(full_path4_bao,"r");}

if(rsd==1){
        f1rsd=fopen(full_path_rsd,"r");
        f2rsd=fopen(full_path2_rsd,"r");
        f3rsd=fopen(full_path3_rsd,"r");
        f4rsd=fopen(full_path4_rsd,"r");}

        fail=0;
        if(f1bao==NULL && f2bao==NULL && f3bao==NULL && f4bao==NULL && bao==1 && rsd==0 )
        {
          printf("Warning, realization %d missing at %s.\n",j,path_to_mocks_bao);fail=1;
        }
        if(f1rsd==NULL && f2rsd==NULL && f3rsd==NULL && f4rsd==NULL && bao==0 && rsd==1 )
        {
          printf("Warning, realization %d missing at %s.\n",j,path_to_mocks_rsd);fail=1;
        }

        if( f1bao==NULL && f2bao==NULL && f3bao==NULL && f4bao==NULL && f1rsd==NULL && f2rsd==NULL && f3rsd==NULL && f4rsd==NULL && bao==1 && rsd==1 )
        {
          printf("Warning, realization %d missing at %s or %s.\n",j,path_to_mocks_bao, path_to_mocks_rsd);fail=1;
        }

if(f1bao!=NULL && bao==1){fclose(f1bao);modebao=1;}
if(f2bao!=NULL && bao==1){fclose(f2bao);modebao=2;}
if(f3bao!=NULL && bao==1){fclose(f3bao);modebao=3;}
if(f4bao!=NULL && bao==1){fclose(f4bao);modebao=4;}
if(f1rsd!=NULL && rsd==1){fclose(f1rsd);modersd=1;}
if(f2rsd!=NULL && rsd==1){fclose(f2rsd);modersd=2;}
if(f3rsd!=NULL && rsd==1){fclose(f3rsd);modersd=3;}
if(f4rsd!=NULL && rsd==1){fclose(f4rsd);modersd=4;}


        if(fail==0)
        {

//iterate one or two times
for(iteration=iteration_ini;iteration<=iteration_fin;iteration++){
  
           Nlines=0;
if(iteration==1){
           if(modebao==1){Nlines=countlines(full_path_bao);}//if(j==1508){printf("1: %d\n",Nlines);Nlines=0;}
           if(modebao==2){Nlines=countlines(full_path2_bao);}//if(j==1508){printf("2: %d\n",Nlines);Nlines=0;}
           if(modebao==3){Nlines=countlines(full_path3_bao);}//if(j==1508){printf("3: %d\n",Nlines);Nlines=0;}
           if(modebao==4){Nlines=countlines(full_path4_bao);}//if(j==1508){printf("4: %d\n",Nlines);Nlines=0;}
}

if(iteration==2){
           if(modersd==1){Nlines=countlines(full_path_rsd);}
           if(modersd==2){Nlines=countlines(full_path2_rsd);}
           if(modersd==3){Nlines=countlines(full_path3_rsd);}
           if(modersd==4){Nlines=countlines(full_path4_rsd);}}

//printf("%d %d (%d)\n",j,Nlines,modersd);
        

           if(Nlines!=0)
           {  
              if(bao==1 && rsd==1){
              if(iteration==1){Ncounter++;}}//count them only once
              else{Ncounter++;}

              i0=0;i2=0;i4=0;l0=0;l2=0;l4=0;
if(iteration==1){
if(modebao==1){f=fopen(full_path_bao,"r");}
if(modebao==2){f=fopen(full_path2_bao,"r");}
if(modebao==3){f=fopen(full_path3_bao,"r");}
if(modebao==4){f=fopen(full_path4_bao,"r");}
}
if(iteration==2){
if(modersd==1){f=fopen(full_path_rsd,"r");}
if(modersd==2){f=fopen(full_path2_rsd,"r");}
if(modersd==3){f=fopen(full_path3_rsd,"r");}
if(modersd==4){f=fopen(full_path4_rsd,"r");}
}


fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s %*f %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*f %*s\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*f\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*f\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*d\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*s\n");

if(iteration==1){kminP0=kminP0bao;kminP2=kminP2bao;kminP4=kminP4bao;kmaxP0=kmaxP0bao;kmaxP2=kmaxP2bao;kmaxP4=kmaxP4bao;}
if(iteration==2){kminP0=kminP0rsd;kminP2=kminP2rsd;kminP4=kminP4rsd;kmaxP0=kmaxP0rsd;kmaxP2=kmaxP2rsd;kmaxP4=kmaxP4rsd;}


              for(i=0;i<Nlines;i++)
              {
                  fscanf(f,"%*f %lf %lf %lf %lf %*d %*f\n",&k,&p0,&p2,&p4);


                  if(k>kminP0 && k<kmaxP0)
                  {
                      if(iteration==1){P0bao[Ncounter-1][i0]=p0;}
                      if(iteration==2){P0rsd[Ncounter-1][i0]=p0;}

                      if(iteration==1){
                      if( strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0){Pvector[Ncounter-1][l0]=p0;}
                      i0++;
                      if( strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0){l0++;}
                      }
                      if(iteration==2){
                      if( strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0){Pvector[Ncounter-1][l0+baoshift]=p0;}
                      i0++;
                      if( strcmp(fit_RSD, "P0") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0){l0++;}
                      }

                      
                  }                 

                  if( strcmp(type_BAORSD, "BAOISO") == 0 || strcmp(type_BAORSD, "FSBAOISO") == 0)
                  {
                       
                      if(k>kminP2 && k<kmaxP2)
                      {
                         if(iteration==1){P2bao[Ncounter-1][i2]=p0+2./5.*p2;}
                         if(iteration==2){P2rsd[Ncounter-1][i2]=p2;}

                         if(iteration==1){
                         if( strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P24") == 0){Pvector[Ncounter-1][l2]=p0+2./5.*p2;}
                         if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P024") == 0){Pvector[Ncounter-1][l2+NeffP0bao]=p0+2./5.*p2;}
                         i2++;
                         if( strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){l2++;}
                         }
                         if(iteration==2){
                         if( strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P24") == 0){Pvector[Ncounter-1][l2+baoshift]=p0+2./5.*p2;}
                         if( strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P024") == 0){Pvector[Ncounter-1][l2+NeffP0rsd+baoshift]=p2;}
                         i2++;
                         if( strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0){l2++;}

                         }

                      }
                      if(k>kminP4 && k<kmaxP4)
                      {
                         if(iteration==1){P4bao[Ncounter-1][i4]=p0+4./7.*p2+8./63.*p4;}
                         if(iteration==2){P4rsd[Ncounter-1][i4]=p4;}

                         if(iteration==1){
                         if( strcmp(fit_BAO, "P4") == 0){Pvector[Ncounter-1][l4]=p0+4./7.*p2+8./63.*p4;}
                         if( strcmp(fit_BAO, "P04") == 0){Pvector[Ncounter-1][l4+NeffP0bao]=p0+4./7.*p2+8./63.*p4;}
                         if( strcmp(fit_BAO, "P24") == 0){Pvector[Ncounter-1][l4+NeffP2bao]=p0+4./7.*p2+8./63.*p4;}
                         if( strcmp(fit_BAO, "P024") == 0){Pvector[Ncounter-1][l4+NeffP0bao+NeffP2bao]=p0+4./7.*p2+8./63.*p4;}
                         i4++;
                         if( strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){l4++;}
                         }
                         if(iteration==2){
                         if( strcmp(fit_RSD, "P4") == 0){Pvector[Ncounter-1][l4+baoshift]=p4;}
                         if( strcmp(fit_RSD, "P04") == 0){Pvector[Ncounter-1][l4+NeffP0rsd+baoshift]=p4;}
                         if( strcmp(fit_RSD, "P24") == 0){Pvector[Ncounter-1][l4+NeffP2rsd+baoshift]=p4;}
                         if( strcmp(fit_RSD, "P024") == 0){Pvector[Ncounter-1][l4+NeffP0rsd+NeffP2rsd+baoshift]=p4;}
                         i4++;
                         if( strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0){l4++;}

                         }


                      }

                  }

                  if( strcmp(type_BAORSD, "FS") == 0 || strcmp(type_BAORSD, "BAOANISO") == 0 || strcmp(type_BAORSD, "FSBAOANISO") == 0)
                  {

                      if(k>kminP2 && k<kmaxP2)
                      {
                         if(iteration==1){P2bao[Ncounter-1][i2]=p2;}
                         if(iteration==2){P2rsd[Ncounter-1][i2]=p2;}

                         if(iteration==1){
                         if( strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P24") == 0){Pvector[Ncounter-1][l2]=p2;}
                         if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P024") == 0){Pvector[Ncounter-1][l2+NeffP0bao]=p2;}

                         i2++;
                         if( strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){l2++;}
                         }
                         if(iteration==2){
                         if( strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P24") == 0){Pvector[Ncounter-1][l2+baoshift]=p2;}
                         if( strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P024") == 0){Pvector[Ncounter-1][l2+NeffP0rsd+baoshift]=p2;}

                         i2++;
                         if( strcmp(fit_RSD, "P2") == 0 || strcmp(fit_RSD, "P02") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0){l2++;}

                         }

                      }
                      if(k>kminP4 && k<kmaxP4)
                      {
                         if(iteration==1){P4bao[Ncounter-1][i4]=p4;}
                         if(iteration==2){P4rsd[Ncounter-1][i4]=p4;}

                         if(iteration==1){
	                 if( strcmp(fit_BAO, "P4") == 0){Pvector[Ncounter-1][l4]=p4;}
                         if( strcmp(fit_BAO, "P04") == 0){Pvector[Ncounter-1][l4+NeffP0bao]=p4;}
                         if( strcmp(fit_BAO, "P24") == 0){Pvector[Ncounter-1][l4+NeffP2bao]=p4;}
                         if( strcmp(fit_BAO, "P024") == 0){Pvector[Ncounter-1][l4+NeffP0bao+NeffP2bao]=p4;}
                         i4++;
                         if( strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0){l4++;}
                         }
                         if(iteration==2)
                         {
                         if( strcmp(fit_RSD, "P4") == 0){Pvector[Ncounter-1][l4+baoshift]=p4;}
                         if( strcmp(fit_RSD, "P04") == 0){Pvector[Ncounter-1][l4+NeffP0rsd+baoshift]=p4;}
                         if( strcmp(fit_RSD, "P24") == 0){Pvector[Ncounter-1][l4+NeffP2rsd+baoshift]=p4;}
                         if( strcmp(fit_RSD, "P024") == 0){Pvector[Ncounter-1][l4+NeffP0rsd+NeffP2rsd+baoshift]=p4;}
                         i4++;
                         if( strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P024") == 0){l4++;}
                         }

                      }
                  }

              }//end of i-for
              //if(mode==1){fclose(f1);}
              //if(mode==2){fclose(f2);}
              //if(mode==3){fclose(f3);}
              //if(mode==4){fclose(f4);}
                fclose(f);              
//exit(0);
if(iteration==1){l0bao=l0;l2bao=l2;l4bao=l4;i0bao=i0;i2bao=i2;i4bao=i4;}
if(iteration==2){l0rsd=l0;l2rsd=l2;l4rsd=l4;i0rsd=i0;i2rsd=i2;i4rsd=i4;}


           }//lines !=0
           else
           {
               printf("Warning, realization %d with no data at %s or %s\n",j,path_to_mocks_bao,path_to_mocks_rsd);
           }
          
}//iteration=1,2 loop

        }//if fail=0

   }//end of j-for


//Get the mean of P0,P2 and P4
//exit(0);
//printf("%d, %d %d %d, %d, %d, %d, %d, %d, %d\n",Ncounter,i0rsd,i2rsd,i4rsd,l0bao,l2bao,l4bao,l0rsd,l2rsd,l4rsd);
for(j=0;j<Ncounter;j++)
{

if(bao==1){
    for(i=0;i<i0bao;i++)
    {
        P0baoav[i]=P0baoav[i]+P0bao[j][i]/Ncounter*1.;
    }
    for(i=0;i<i2bao;i++)
    {
        P2baoav[i]=P2baoav[i]+P2bao[j][i]/Ncounter*1.;

    }
    for(i=0;i<i4bao;i++)
    {
        P4baoav[i]=P4baoav[i]+P4bao[j][i]/Ncounter*1.;

    }
}
if(rsd==1){
    for(i=0;i<i0rsd;i++)
    {
        P0rsdav[i]=P0rsdav[i]+P0rsd[j][i]/Ncounter*1.;
    }
    for(i=0;i<i2rsd;i++)
    {
        P2rsdav[i]=P2rsdav[i]+P2rsd[j][i]/Ncounter*1.;

    }
    for(i=0;i<i4rsd;i++)
    {
        P4rsdav[i]=P4rsdav[i]+P4rsd[j][i]/Ncounter*1.;

    }
}
//exit(0);
    for(i=0;i<l0bao+l2bao+l4bao+l0rsd+l2rsd+l4rsd;i++)
    {
        Pvectorav[i]=Pvectorav[i]+Pvector[j][i]/Ncounter*1.;
    }

//exit(0);
}
//exit(0);
//for(i=0;i<i0;i++){printf("P[%d]=%lf,%lf,%lf\n",i,P0baoav[i],P2baoav[i],P4baoav[i]);}

}//end of do power spectrum if

//printf("%d\n",Nrealizations);
//exit(0);
if(Ncounter!=Nrealizations){printf("Warning, %d/%d realizations used for Power Spectrum.\n",Ncounter,Nrealizations);}
//exit(0);
/*
if( strcmp(do_bispectrum, "yes") == 0)
{

if(bao==1){

B0baoav=  (double*) calloc( NeffB0bao, sizeof(double));

B0bao = (double**) calloc(Nrealizations,sizeof(double*));

for(j=0;j<Nrealizations;j++)
{
   B0bao[j] = (double*) calloc(NeffB0bao,sizeof(double));
}

}

if(rsd==1){

B0rsdav=  (double*) calloc( NeffB0rsd, sizeof(double));

B0rsd = (double**) calloc(Nrealizations,sizeof(double*));

for(j=0;j<Nrealizations;j++)
{
   B0rsd[j] = (double*) calloc(NeffB0rsd,sizeof(double));
}

}


Ncounterbis=0;
   for(j=1;j<=Nrealizations;j++)
   {
        sprintf(full_path,"%s%.4d.txt",path_to_mocks_bis,j);
        sprintf(full_path2,"%s%d.txt",path_to_mocks_bis,j);
        sprintf(full_path3,"%s%.4d.dat",path_to_mocks_bis,j);
        sprintf(full_path4,"%s%d.dat",path_to_mocks_bis,j);

        f1=fopen(full_path,"r");
        f2=fopen(full_path2,"r");
        f3=fopen(full_path3,"r");
        f4=fopen(full_path4,"r");

        if(f1==NULL && f2==NULL && f3==NULL && f4==NULL)
        {
          //do nothing
        }
        else
        {
           Nlines=0;
           if(f1!=NULL){mode=1;Nlines=countlines(full_path);}
           if(f2!=NULL){mode=2;Nlines=countlines(full_path2);}
           if(f3!=NULL){mode=3;Nlines=countlines(full_path3);}
           if(f4!=NULL){mode=4;Nlines=countlines(full_path4);}

           if(Nlines!=0)
           {
              Ncounterbis++;
              i0bis=0;l=0;
              for(i=0;i<Nlines;i++)
              {
                  if(mode==1){fscanf(f1,"%*f %lf %*f %lf %*f %lf %lf %*f %*f %*f %*f\n",&k1,&k2,&k3,&b0);}
                  if(mode==2){fscanf(f2,"%*f %lf %*f %lf %*f %lf %lf %*f %*f %*f %*f\n",&k1,&k2,&k3,&b0);}
                  if(mode==3){fscanf(f3,"%*f %lf %*f %lf %*f %lf %lf %*f %*f %*f %*f\n",&k1,&k2,&k3,&b0);}
                  if(mode==4){fscanf(f4,"%*f %lf %*f %lf %*f %lf %lf %*f %*f %*f %*f\n",&k1,&k2,&k3,&b0);}

                  if(k1>kminB0 && k1<kmaxB0 && k2>kminB0 && k2<kmaxB0 && k3>kminB0 && k3<kmaxB0)
                  {
                      B0[Ncounterbis-1][i0bis]=b0;
                      if( strcmp(do_power_spectrum, "no") == 0){Pvector[Ncounterbis-1][l]=b0;}
                      else
                      {
                             if( strcmp(fit_BAO, "P0") == 0){Pvector[Ncounterbis-1][l+NeffP0]=b0;}
                             if( strcmp(fit_BAO, "P2") == 0){Pvector[Ncounterbis-1][l+NeffP2]=b0;}
                             if( strcmp(fit_BAO, "P4") == 0){Pvector[Ncounterbis-1][l+NeffP4]=b0;}
                             if( strcmp(fit_BAO, "P02") == 0){Pvector[Ncounterbis-1][l+NeffP0+NeffP2]=b0;}
                             if( strcmp(fit_BAO, "P04") == 0){Pvector[Ncounterbis-1][l+NeffP0+NeffP4]=b0;}
                             if( strcmp(fit_BAO, "P24") == 0){Pvector[Ncounterbis-1][l+NeffP2+NeffP4]=b0;}
                             if( strcmp(fit_BAO, "P024") == 0){Pvector[Ncounterbis-1][l+NeffP0+NeffP2+NeffP4]=b0;}

                      }
                      i0bis++;
                      l++;
                  }

              }//end of i-for
              if(mode==1){fclose(f1);}
              if(mode==2){fclose(f2);}
              if(mode==3){fclose(f3);}
              if(mode==4){fclose(f4);}

           }//lines !=0
           else
           {
            
           }

        }//else f=open not null

   }//end of j-for

//Get mean of B0

for(j=0;j<Ncounterbis;j++)
{
    for(i=0;i<i0bis;i++)
    {
        B0av[i]=B0av[i]+B0[j][i]/Ncounterbis*1.;
    }

    for(i=0;i<l;i++)
    {

           if( strcmp(do_power_spectrum, "no") == 0){Pvectorav[i]=Pvectorav[i]+B0[j][i]/Ncounterbis*1.;}
           else
           {
               if( strcmp(fit_BAO, "P0") == 0){Pvectorav[i+NeffP0]=Pvectorav[i+NeffP0]+B0[j][i+NeffP0]/Ncounterbis*1.;}
               if( strcmp(fit_BAO, "P2") == 0){Pvectorav[i+NeffP2]=Pvectorav[i+NeffP2]+B0[j][i+NeffP2]/Ncounterbis*1.;}
               if( strcmp(fit_BAO, "P4") == 0){Pvectorav[i+NeffP4]=Pvectorav[i+NeffP4]+B0[j][i+NeffP4]/Ncounterbis*1.;}
               if( strcmp(fit_BAO, "P02") == 0){Pvectorav[i+NeffP0+NeffP2]=Pvectorav[i+NeffP0+NeffP2]+B0[j][i+NeffP0+NeffP2]/Ncounterbis*1.;}
               if( strcmp(fit_BAO, "P04") == 0){Pvectorav[i+NeffP0+NeffP4]=Pvectorav[i+NeffP0+NeffP4]+B0[j][i+NeffP0+NeffP4]/Ncounterbis*1.;}
               if( strcmp(fit_BAO, "P24") == 0){Pvectorav[i+NeffP2+NeffP4]=Pvectorav[i+NeffP2+NeffP4]+B0[j][i+NeffP2+NeffP4]/Ncounterbis*1.;}
               if( strcmp(fit_BAO, "P024") == 0){Pvectorav[i+NeffP0+NeffP2+NeffP4]=Pvectorav[i+NeffP0+NeffP2+NeffP4]+B0[j][i+NeffP0+NeffP2+NeffP4]/Ncounterbis*1.;}

            }


    }
   
   
}

if(Ncounterbis!=Nrealizations){printf("Warning, %d/%d realizations used for Bispectrum.\n",Ncounter,Nrealizations);}



}//do bispectrum if


if( strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 && Ncounter!=Ncounterbis){printf("Error, number of used realizations is different for the power spectrum and bispectrum: %d != %d. Exiting now...\n",Ncounter,Ncounterbis);exit(0);} 

if( strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0){Ncounter=Ncounterbis;}
*/
//Build errors

//printf("%d %d\n",bao,rsd);
//exit(0);
if( strcmp(do_power_spectrum, "yes") == 0)
{

if( bao==1){

  for(i_real=0;i_real<Ncounter;i_real++)
  {
    for(i=0;i<i0bao;i++)
    {
        errP0bao[i]=errP0bao[i]+pow( P0baoav[i]-P0bao[i_real][i],2);
    }

    for(i=0;i<i2bao;i++)
    {
        errP2bao[i]=errP2bao[i]+pow( P2baoav[i]-P2bao[i_real][i],2);
    }

    for(i=0;i<i4bao;i++)
    {
        errP4bao[i]=errP4bao[i]+pow( P4baoav[i]-P4bao[i_real][i],2);
    }

  }

    for(i=0;i<i0bao;i++)
    {
        errP0bao[i]=sqrt(errP0bao[i]/(Ncounter-1.));
    }
    for(i=0;i<i2bao;i++)
    {
        errP2bao[i]=sqrt(errP2bao[i]/(Ncounter-1.));

    }
    for(i=0;i<i4bao;i++)
    {
        errP4bao[i]=sqrt(errP4bao[i]/(Ncounter-1.));
    }
}

//printf("%d %d %d\n",i0rsd,i2rsd,i4rsd);
//exit(0);
if( rsd==1){

  for(i_real=0;i_real<Ncounter;i_real++)
  {
    for(i=0;i<i0rsd;i++)
    {
        errP0rsd[i]=errP0rsd[i]+pow( P0rsdav[i]-P0rsd[i_real][i],2);
    }

    for(i=0;i<i2rsd;i++)
    {
        errP2rsd[i]=errP2rsd[i]+pow( P2rsdav[i]-P2rsd[i_real][i],2);
    }

    for(i=0;i<i4rsd;i++)
    {
        errP4rsd[i]=errP4rsd[i]+pow( P4rsdav[i]-P4rsd[i_real][i],2);
    }

  }

    for(i=0;i<i0rsd;i++)
    {
        errP0rsd[i]=sqrt(errP0rsd[i]/(Ncounter-1.));
    }
    for(i=0;i<i2rsd;i++)
    {
        errP2rsd[i]=sqrt(errP2rsd[i]/(Ncounter-1.));

    }
    for(i=0;i<i4rsd;i++)
    {
        errP4rsd[i]=sqrt(errP4rsd[i]/(Ncounter-1.));
    }


}

}
//exit(0);
/*
if( strcmp(do_bispectrum, "yes") == 0)
{

  for(i_real=0;i_real<Ncounter;i_real++)
  {
    for(i=0;i<i0bis;i++)
    {
        errB0[i]=errB0[i]+pow( B0av[i]-B0[i_real][i],2);
    }

  }

    for(i=0;i<i0bis;i++)
    {
        errB0[i]=sqrt(errB0[i]/(Ncounter-1.));
    }

}
*/
//Build covariance

for(i_real=0;i_real<Ncounter;i_real++)
{
for(j=0;j<Ncov;j++)
{
jelement=Pvectorav[j]-Pvector[i_real][j];
for(i=0;i<Ncov;i++)
{
ielement=Pvectorav[i]-Pvector[i_real][i];
Variance[j+i*Ncov]=Variance[j+i*Ncov]+ielement*jelement/(Ncounter*1.-1.);
}
}
}

//write somewhere the covariance used
/*
for(j=0;j<Ncov;j++)
{
for(i=0;i<Ncov;i++)
{
if(i!=Ncov-1){printf("%lf\t",Variance[j+i*Ncov]/sqrt(Variance[i+i*Ncov]*Variance[j+j*Ncov]));}
else{printf("%lf\n",Variance[j+i*Ncov]/sqrt(Variance[i+i*Ncov]*Variance[j+j*Ncov]));}
}
}
*/
//exit(0);
       //invert covariance

        gsl_matrix * m = gsl_matrix_alloc (Ncov, Ncov);
        gsl_matrix * inverse = gsl_matrix_alloc (Ncov, Ncov);
        gsl_matrix * identity = gsl_matrix_alloc (Ncov, Ncov);

        gsl_permutation * perm = gsl_permutation_alloc (Ncov);
         for(i=0;i<Ncov;i++)
         {
         for(j=0;j<Ncov;j++)
         {
                if(Variance[i+i*Ncov]==0 || Variance[j+j*Ncov]==0){printf("Error 0-element(s) in the diagonal of covarince: Exiting now...\n");exit(0);}
                gsl_matrix_set (m, i, j, Variance[j+i*Ncov]);
         }
         }
           gsl_linalg_LU_decomp (m, perm, &s);
           gsl_linalg_LU_invert (m, perm, inverse);

         for(i=0;i<Ncov;i++)
         {
         for(j=0;j<Ncov;j++)
         {
           cov[j+Ncov*i]=1./((1.-(Ncov+1.)/(Ncounter*1.-1.))*gsl_matrix_get (inverse, i, j))*(1./scaling_factor);
         }
         }
         printf("Hartlap correction of covariance is %lf, %d x %d elements ( %d for P0-BAO, %d for P2-BAO, %d for P4-BAO, %d for P0-RSD, %d for P2-RSD, %d for P4-RSD; %d realizations)\n",1./((1.-(Ncov+1.)/(Ncounter*1.-1.))),Ncov,Ncov,l0bao,l2bao,l4bao,l0rsd,l2rsd,l4rsd,Ncounter);
         if(Ncov>=Ncounter){printf("Error, not enough realizations to perform a reliable estimate of the covariance. Please reduce the dimension of your covariance or add more mock realizations. Exiting now...\n");exit(0);}

//free stuff!
free(Variance);
free(Pvectorav);
freeTokens(Pvector,Nrealizations);

gsl_matrix_free(m);
gsl_matrix_free(inverse);
gsl_matrix_free(identity);
gsl_permutation_free(perm);

if( strcmp(do_power_spectrum, "yes") == 0)
{
if(bao==1){
free(P0baoav);
free(P2baoav);
free(P4baoav);
freeTokens(P0bao,Nrealizations);
freeTokens(P2bao,Nrealizations);
freeTokens(P4bao,Nrealizations);}

if(rsd==1){
free(P0rsdav);
free(P2rsdav);
free(P4rsdav);
freeTokens(P0rsd,Nrealizations);
freeTokens(P2rsd,Nrealizations);
freeTokens(P4rsd,Nrealizations);}


}/*
if( strcmp(do_bispectrum, "yes") == 0)
{
free(B0av);
freeTokens(B0,Nrealizations);
}
*/
}//end of get_cov

void get_mask(char *path, double posAV[], double pos[],double W0[], double W2[], double W4[], double W6[], double W8[], int Nmask, char *type_BAORSD, double params[],char *renormalize_window)
{
int i;
FILE *f;
double p,pav,w0,w2,w4,w6,w8;
double Norm,sumw_ran,sumw_dat;
double correction,I22;
sumw_dat=params[1];
I22=params[2];
f=fopen(path,"r");
fscanf(f,"%*s %lf %*s %lf\n",&Norm,&sumw_ran);
if( strcmp(renormalize_window,"yes") == 0){
correction=Norm*pow(sumw_dat/sumw_ran,2)/(2.*Pi*I22);
printf("Correcting %s window by factor %lf\n",path,correction);
if(correction<=0.01){printf("Error with correction: %lf. Norm=%lf, sumw_dat=%lf, sumw_ran=%lf, I22=%lf. Exiting now...\n",correction,Norm,sumw_dat,sumw_ran,I22);exit(0);}
if(correction>2){printf("Warning. Large correction applied: %lf. Norm=%lf, sumw_dat=%lf, sumw_ran=%lf, I22=%lf. Exiting now...\n",correction,Norm,sumw_dat,sumw_ran,I22);exit(0);}
}
if( strcmp(renormalize_window,"no") == 0){correction=1;}
for(i=0;i<Nmask;i++)
{
fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %*f\n",&pav,&p,&w0,&w2,&w4,&w6,&w8);

if( isnan(w0)!=0 || isnan(w2)!=0 || isnan(w4)!=0 || isnan(w6)!=0 || isnan(w8)!=0 || isnan(p)!=0)
{
break;
}
else
{
posAV[i]=pav;
pos[i]=p;
W0[i]=w0*correction;
W2[i]=w2*correction;
W4[i]=w4*correction;
W6[i]=w6*correction;
W8[i]=w8*correction;
}

}
fclose(f);
if(i!=Nmask){printf("Warning, only %d/%d first rows used from file %s\n",i,Nmask,path);}
params[0]=i;
}

