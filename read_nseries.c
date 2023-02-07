#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"
#include "get_line.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#define Pi (4.*atan(1.))
#define SCALING_COVARIANCE (1.12*84)
int get_compression_file(char *path_to_compression, double **matrix_linear_compression, int input_elements_data_vector)
{
int input,output,columns;
FILE *f;
int N,Nall,i,j,ch;
char dummy[200];
int dummy2,read_header1,read_header2;
N=countlines(path_to_compression);//lines with entries (output elements)
Nall=count_all_lines(path_to_compression);//lines with entries + header

f=fopen(path_to_compression,"r");

for(i=0;i<Nall-N;i++){

read_header1=0;
read_header2=0;
do{

fscanf(f,"%s",dummy);
ch = fgetc(f);

dummy2=atoi(dummy);

if(read_header1==1 && dummy2>0){   input=dummy2; read_header1=0;}
if(strcmp(dummy,"#Inputs/columns") == 0){read_header1++;}

if(read_header2==1 && dummy2>0){   output=dummy2; read_header2=0;}
if(strcmp(dummy,"#Outputs/rows") == 0){read_header2++;}

}while(ch != '\n');

}//skip header line

if(input!=input_elements_data_vector){printf("Warning, the inputed elemens in the file header (%d), do not match with the actual data-vector number of elements (%d). Exiting now...\n",input,input_elements_data_vector);exit(0);}

if(output!=N){printf("Warning, the number of rows of the file (%d), do not match with the output elements in the file header (%d). Exiting now...\n",N,output);exit(0);}

fclose(f);

f=fopen(path_to_compression,"r");

for(i=0;i<Nall-N;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line

//check number of columns that correspond to  input value
columns=0;
do{fscanf(f,"%s",dummy);ch = fgetc(f);columns++;}while(ch != '\n');
fclose(f);

if(columns!=input){printf("Warning, the input elements in the covariance (number of columns, %d) do not correspond to the stated number in the header (%d). Exiting now...\n",columns,input);exit(0);}

if(input==output){printf("Warning, no compression applied. Same number of elements after 'compressing'.\n"); }
if(input>output){printf("Compressing from %d ---> %d, (%.1lf %)...\n",input,output,100.*input/output*1.); }
if(input<output){printf("Warning, the data-vector has more elements after compressing (so not a compression). Check your files. Exiting now...\n");exit(0); }

if(matrix_linear_compression==NULL){return output;}//matrix has not been yet allocated. Return output value to do so
else{//matrix has been allocated

f=fopen(path_to_compression,"r");

for(i=0;i<Nall-N;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line

for(i=0;i<output;i++)
{

     for(j=0;j<input;j++)
     {
         if(j==input-1){fscanf(f,"%lf\n",&matrix_linear_compression[i][j]);}
         else{fscanf(f,"%lf\n",&matrix_linear_compression[i][j]);}  
     }
//read compressed covariance

}

fclose(f);
return 0;
}

}


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
int Neff,N,i,j;
char dummy[200];
int ch;
double k;
double k1,k2,k3;
FILE *f;
double params[10];
int Nall;//N=countlines(path)-27;
N=countlines(path);//lines with entries
Nall=count_all_lines(path);//lines with entries + header
f=fopen(path,"r");

for(i=0;i<Nall-N;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line



Neff=0;
for(i=0;i<N;i++)
{

if(mode==0){// power spectrum file structure for data

get_lineP(f, params,0);
k=params[1];
if(k>kmin && k<kmax){/*printf("k=%lf\n",k);*/Neff++;}
}

if(mode==1){//bispectrum file structure for data
get_lineB(f, params,0);
k1=params[1];
k2=params[3];
k3=params[5];
if(k1>kmin && k1<kmax && k2>kmin && k2<kmax && k3>kmin && k3<kmax){/*printf("%d) %lf %lf %lf\n",i,k1,k2,k3);*/Neff++;}
}

}
fclose(f);

return Neff;
}

void get_Pk_bao( char *path, double k[], double P[])
{
int N,Nall,i,ch;
FILE *f;
N=countlines(path);//lines with entries
Nall=count_all_lines(path);//lines with entries + header
f=fopen(path,"r");

for(i=0;i<Nall-N;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line

for(i=0;i<N;i++)
{
fscanf(f,"%lf %lf\n",&k[i],&P[i]);
}
fclose(f);
}

void get_data_alphas(char *path, double *P0bao)
{
FILE *f;
int ch;
double params[10];
f=fopen(path,"r");
do{ch = fgetc(f);}while(ch != '\n');//skip 1st line (header)
  get_linealphas(f,params,0);
  P0bao[0]=params[2];
  get_linealphas(f,params,0);
  P0bao[1]=params[2];
fclose(f);
}


void get_data(char *path, double k0[], double kav0[], double P0[], double k2[], double kav2[], double P2[], double k4[], double kav4[], double P4[], double parameter_value[], double kminP0,double kmaxP0, double kminP2, double kmaxP2, double kminP4, double kmaxP4, char *type_BAORSD, int baorsd, int bao)
{
int read_header1,read_header2;
char dummy[200];
double dummy2;
double params[10];
int i,N,Nall,ch;
int i0,i2,i4;
FILE *f;
double Pnoise,Pnoise_prev,sumw,I22;sumw=0;I22=0;
double k,p0,p2,p4,kav;
//N=countlines(path)-27;
N=countlines(path);
Nall=count_all_lines(path);//lines with entries + header
f=fopen(path,"r");

for(i=0;i<Nall-N;i++){

read_header1=0;
read_header2=0;
do{

fscanf(f,"%s",dummy); 
ch = fgetc(f);

dummy2=atof(dummy);

if(read_header1==1 && dummy2>0){   I22=dummy2; read_header1=0;}
if(strcmp(dummy,"#Normalization") == 0){read_header1++;}

if(read_header2==8 && dummy2>0){ sumw=dummy2; read_header2=0;}
if(strcmp(dummy,"#Effective") == 0){read_header2++;}
if(strcmp(dummy,"number") == 0){read_header2++;}
if(strcmp(dummy,"of") == 0){read_header2++;}
if(strcmp(dummy,"data") == 0){read_header2++;}
if(strcmp(dummy,"elements") == 0){read_header2++;}
if(strcmp(dummy,"weighted") == 0){read_header2++;}
if(strcmp(dummy,"by") == 0){read_header2++;}
if(strcmp(dummy,"wcol*wsys*wfkp:") == 0){read_header2++;}

}while(ch != '\n');


}//skip header line
if(I22>0 && sumw>0){printf("I22=%lf, Sumw=%lf\n",I22,sumw);}

fclose(f);

f=fopen(path,"r");

for(i=0;i<Nall-N;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line  
 
   i0=0;i2=0;i4=0;
   for(i=0;i<N;i++)
   {
  get_lineP(f,params,0);
  kav=params[0];
  k=params[1];
  p0=params[2];
  p2=params[3];
  p4=params[4];
  Pnoise=params[5];

  if(i>0 && Pnoise!=Pnoise_prev){printf("Warning, shot noise value changes across the file! Pnoise[%d]=%lf; Pnoise[%d]=%lf. Exiting now...\n",i,Pnoise,i-1,Pnoise_prev);exit(0);}
  Pnoise_prev=Pnoise;

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
        if(strcmp(type_BAORSD, "FS") == 0 || strcmp(type_BAORSD, "BAOANISO") == 0 ||  strcmp(type_BAORSD, "FSBAOANISO") == 0 || strcmp(type_BAORSD, "FSalphasrecon") == 0){

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
}

void get_data_bis(char *path, double k11[], double k22[], double k33[], double B0[], double Bnoise[], double kminB0, double kmaxB0,char *bispectrumBQ)
{
double params[10];
int i,N,Nall,ch;
int i0;
FILE *f;
double k1,k2,k3,b0,noise,q0,noiseq;
N=countlines(path);//lines with entries
Nall=count_all_lines(path);//lines with entries + header
f=fopen(path,"r");

for(i=0;i<Nall-N;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line

   i0=0;
   for(i=0;i<N;i++)
   {
       get_lineB(f,params,0);
       k1=params[1];
       k2=params[3];
       k3=params[5];
       b0=params[6];
       noise=params[7];
       q0=params[8];
       noiseq=params[9];
        if(k1>kminB0 && k1<kmaxB0 && k2>kminB0 && k2<kmaxB0 && k3>kminB0 && k3<kmaxB0)
        {
           k11[i0]=k1;
           k22[i0]=k2;
           k33[i0]=k3;
         if(strcmp(bispectrumBQ,"B") == 0){
           B0[i0]=b0;
           Bnoise[i0]=noise;}
if(strcmp(bispectrumBQ,"Q") == 0){
//          printf("error not Q available. Exiting now...\n");exit(0);//comment this line if Q actually available
          B0[i0]=q0;
          Bnoise[i0]=noiseq;
        }
           i0++;
        }
   }
   fclose(f);

}


void get_cov_from_file(char *path_to_cov, char *path_to_cov_bis, double cov[], int Ncov,int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao,  int NeffP2rsd, int NeffP4bao, int NeffP4rsd , int NeffB0bao, int NeffB0rsd, double errP0bao[], double errP0rsd[],double errP2bao[], double errP2rsd[], double errP4bao[],  double errP4rsd[], double errB0bao[], double errB0rsd[], char *do_power_spectrum, char *do_bispectrum, double **matrix_compression, int input, int output, char *covariance_correction)
{

FILE *f;
double scaling_factor=SCALING_COVARIANCE;//number of realizations when fitting the mean & extra scaling
double hartlap_factor;
double k,k1,k2,k3,p0,p2,p4,b0;
int Nlines, i, j, s;

double *Variance;

Nlines = countlines(path_to_cov);

if (Nlines != Ncov)
{
  printf("Error: The dimension of your input covariance of %d does not match the expected number %d corresponding to the chosen fittype. \n. Exiting now...\n",Nlines,Ncov);exit(0);
}

f=fopen(path_to_cov,"r");
skipheader(f);

Variance = (double*) calloc( Ncov*Ncov, sizeof(double));
for(j=0;j<Ncov;j++)
{
for(i=0;i<Ncov;i++)
{
fscanf(f, "%lf", &Variance[j+i*Ncov]);
}
}

if( strcmp(covariance_correction,"none") == 0 || strcmp(covariance_correction,"Sellentin-Heavens") == 0)
{
hartlap_factor = 1.;
}
else{//Hartlap

if (Nrealizations<2){
  hartlap_factor = 1.;
}
else {
  hartlap_factor = 1./(1.-(Ncov+1.)/(Nrealizations*1.-1.));
}

}

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
           cov[j+Ncov*i]=hartlap_factor/(gsl_matrix_get (inverse, i, j))*(1./scaling_factor);
         }
         }

gsl_matrix_free(m);
gsl_matrix_free(inverse);
gsl_matrix_free(identity);
gsl_permutation_free(perm);

i=0;
if(NeffP0bao>0)
{
for(j=0;j<NeffP0bao;j++){errP0bao[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

if(NeffP2bao>0)
{
for(j=0;j<NeffP2bao;j++){errP2bao[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

if(NeffP4bao>0)
{
for(j=0;j<NeffP4bao;j++){errP4bao[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

if(NeffP0rsd>0)
{
for(j=0;j<NeffP0rsd;j++){errP0rsd[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

if(NeffP2rsd>0)
{
for(j=0;j<NeffP2rsd;j++){errP2rsd[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

if(NeffP4rsd>0)
{
for(j=0;j<NeffP4rsd;j++){errP4rsd[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

if(NeffB0bao>0)
{
for(j=0;j<NeffB0bao;j++){errB0bao[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

if(NeffB0rsd>0)
{
for(j=0;j<NeffB0rsd;j++){errB0rsd[j]=sqrt(Variance[i+Ncov*i])/sqrt(scaling_factor);i++;}
}

free(Variance);

}//end of get_cov_from_file


void write_cov_file(char *covfile, double Covariance[], int Ncov, int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao, int NeffP2rsd, int NeffP4bao, int NeffP4rsd, int NeffB0bao, int NeffB0rsd)
{
FILE *f,*g;
int i,j;
char rfile[2000];
sprintf(rfile,"%s.rterms",covfile);
f=fopen(covfile,"w");
g=fopen(rfile,"w");
if(f==NULL){printf("Error. Can't write covariance in %s. Exiting now...\n",covfile);exit(0);}

fprintf(f,"# Nrealizations %d, \t NeffP0bao %d, \t NeffP2bao %d, \t NeffP4bao %d,\t NeffP0rsd %d, \t NeffP2rsd %d, \t NeffP4rsd %d, \t NeffB0bao %d, \t NeffB0rsd %d\n", Nrealizations, NeffP0bao, NeffP2bao, NeffP4bao, NeffP0rsd, NeffP2rsd, NeffP4rsd, NeffB0bao, NeffB0rsd);

for(j=0;j<Ncov;j++)//Write Covariance
{
for(i=0;i<Ncov;i++)
{
if(i!=Ncov-1){fprintf(f,"%e\t",Covariance[j+i*Ncov]);fprintf(g,"%e\t",Covariance[j+i*Ncov]/sqrt(Covariance[i+i*Ncov]*Covariance[j+j*Ncov]));}
else{fprintf(f,"%e\n",Covariance[j+i*Ncov]);fprintf(g,"%e\n",Covariance[j+i*Ncov]/sqrt(Covariance[i+i*Ncov]*Covariance[j+j*Ncov]));}
}
}
fclose(f);
fclose(g);
}

void get_mask(char *path, double posAV[], double pos[],double W0[], double W2[], double W4[], double W6[], double W8[], int Nmask, char *type_BAORSD, double params[],char *renormalize_window)
{
char dummy[200];
double dummy2;
int ch;
int read_header1,read_header2;
int i;
FILE *f;
double paramsmask[10];
double p,pav,w0,w2,w4,w6,w8;
double Norm,sumw_ran,sumw_dat;
double correction,I22;
int N,Nall;
sumw_dat=params[1];
I22=params[2];
f=fopen(path,"r");

Norm=0;sumw_ran=0;

N=countlines(path);//lines with entries
Nall=count_all_lines(path);//lines with entries + header
f=fopen(path,"r");


for(i=0;i<Nall-N;i++){

read_header1=0;
read_header2=0;
do{

fscanf(f,"%s",dummy);
ch = fgetc(f);

dummy2=atof(dummy);

if(read_header1==1 && dummy2>0){   Norm=dummy2; read_header1=0;}
if(strcmp(dummy,"#N=") == 0){read_header1++;}

if(read_header2==1 && dummy2>0){   sumw_ran=dummy2; read_header2=0;}
if(strcmp(dummy,"Sumwtot=") == 0){read_header2++;}

//printf("%s %lf; %d %d; %lf %lf \n",dummy,dummy2,read_header1,read_header2,Norm,sumw_ran);

}while(ch != '\n');


}//skip header line
if(Norm>0 && sumw_ran>0){printf("N_win=%lf, Sumw_ran_win=%lf\n",Norm,sumw_ran);}

//exit(0);

if( strcmp(renormalize_window,"yes") == 0){
correction=Norm*pow(sumw_dat/sumw_ran,2)/(2.*Pi*I22);
printf("Correcting %s window by factor %lf\n",path,correction);

if(correction<=0.01){printf("Error with correction: %lf. Norm=%lf, sumw_dat=%lf, sumw_ran=%lf, I22=%lf. Exiting now...\n",correction,Norm,sumw_dat,sumw_ran,I22);exit(0);}

if(correction>2){printf("Warning. Large correction applied: %lf. Norm=%lf, sumw_dat=%lf, sumw_ran=%lf, I22=%lf. Exiting now...\n",correction,Norm,sumw_dat,sumw_ran,I22);exit(0);}
}

if( strcmp(renormalize_window,"no") == 0){correction=1;}
for(i=0;i<N;i++)
{
get_lineMask(f,paramsmask);
pav=paramsmask[0];
p=paramsmask[1];
w0=paramsmask[2];
w2=paramsmask[3];
w4=paramsmask[4];
w6=paramsmask[5];
w8=paramsmask[6];

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


int get_cov_from_mocks(char *path_to_mocks_bao, char *path_to_mocks_rsd, char *path_to_mocks_bis_bao, char *path_to_mocks_bis_rsd,  char *covfile, double cov[], int Ncov,int Nrealizations, int NeffP0bao, int NeffP0rsd, int NeffP2bao,  int NeffP2rsd, int NeffP4bao, int NeffP4rsd , int NeffB0bao, int NeffB0rsd, double errP0bao[], double errP0rsd[],double errP2bao[], double errP2rsd[], double errP4bao[],  double errP4rsd[] , double errB0bao[], double errB0rsd[], double kminP0bao,  double kminP0rsd,  double kmaxP0bao,  double kmaxP0rsd,  double kminP2bao, double kminP2rsd ,double kmaxP2bao, double kmaxP2rsd,double kminP4bao, double kminP4rsd,double kmaxP4bao, double kmaxP4rsd,double kminB0bao, double kminB0rsd,double kmaxB0bao, double kmaxB0rsd, char *type_BAORSD, char *fit_BAO, char *fit_RSD, char *do_power_spectrum, char *do_bispectrum, char *bispectrumBQ,double  **matrix_compression, int input, int output, char *covariance_correction,char *covarianceFSa_option, double *paramsBAO)
{
int out;
double params[10];
int bao,rsd;
int iteration,iteration_ini,iteration_fin;
double scaling_factor;//number of realizations when fitting the mean
int i,j,Nlines,ch,Nall,Ncounter,Ncounterbis,i0,i2,i4,i0bis,i_real;
int l0rsd,l0bao,l2rsd,l2bao,l4rsd,l4bao;
int i0rsd,i0bao,i2rsd,i2bao,i4rsd,i4bao,i0bisbao,i0bisrsd;
int l,l0,l2,l4;
int s,fail;
int ii,jj;
double jelement,ielement;
FILE *f1bao,*f2bao,*f3bao,*f4bao,*f5bao,*f6bao;
FILE *f1rsd,*f2rsd,*f3rsd,*f4rsd;
FILE *f,*g;
double kminP0,kminP2,kminP4,kminB0;
double kmaxP0,kmaxP2,kmaxP4,kmaxB0;
int modersd=-1;
int modebao=-1;
double k,k1,k2,k3,p0,p2,p4,b0,q0,bispectrum;
int baoshift,rsdshift,Pshift,Bbaoshift;
double diag_apara,diag_aperp;
char full_path_bao[2000];
char full_path2_bao[2000];
char full_path3_bao[2000];
char full_path4_bao[2000];
char full_path5_bao[2000];
char full_path6_bao[2000];

char full_path_rsd[2000];
char full_path2_rsd[2000];
char full_path3_rsd[2000];
char full_path4_rsd[2000];



int nP0bao,nP2bao,nP4bao,nB0bao,nP0rsd, nP2rsd, nP4rsd,nB0rsd;
double *P0baoav,*P2baoav, *P4baoav, *B0baoav, *Pvectorav;
double **P0bao, **P2bao, **P4bao, **B0bao, **Pvector;
double *P0rsdav,*P2rsdav, *P4rsdav, *B0rsdav;
double **P0rsd, **P2rsd, **P4rsd, **B0rsd;
char covfile_in[2000];

double *Variance,*Variance_comp;
double hartlap_factor;

gsl_matrix *m;
gsl_matrix *inverse;
gsl_matrix *identity;
gsl_permutation *perm;

//exit(0);
l0bao=0;
l0rsd=0;
l2bao=0;
l2rsd=0;
l4bao=0;
l4rsd=0;
i0bisbao=0;
i0bisrsd=0;
bao=0;rsd=0;baoshift=0;rsdshift=0;Pshift=0;
if( strcmp(type_BAORSD,"FSBAOISO") ==0 || strcmp(type_BAORSD,"FSBAOANISO") ==0 || strcmp(type_BAORSD,"FSalphasrecon") ==0){iteration_ini=1;iteration_fin=2;bao=1;rsd=1;}
if( strcmp(type_BAORSD,"BAOISO") ==0 || strcmp(type_BAORSD,"BAOANISO") ==0 ){iteration_ini=1;iteration_fin=1;bao=1;rsd=0;}
if( strcmp(type_BAORSD,"FS") ==0){iteration_ini=2;iteration_fin=2;rsd=1;bao=0;}

nP0bao=0;if( strcmp(fit_BAO,"P0") ==0 || strcmp(fit_BAO,"P02") ==0 ||  strcmp(fit_BAO,"P04") ==0 || strcmp(fit_BAO,"P024")==0 ){nP0bao=NeffP0bao;}
nP2bao=0;if( strcmp(fit_BAO,"P2") ==0 || strcmp(fit_BAO,"P02") ==0 ||  strcmp(fit_BAO,"P24") ==0 || strcmp(fit_BAO,"P024")==0 ){nP2bao=NeffP2bao;}
nP4bao=0;if( strcmp(fit_BAO,"P4") ==0 || strcmp(fit_BAO,"P24") ==0 ||  strcmp(fit_BAO,"P24") ==0 || strcmp(fit_BAO,"P024")==0 ){nP4bao=NeffP4bao;}
nB0bao=0;if( strcmp(do_bispectrum,"yes") == 0 && bao == 1){nB0bao=NeffB0bao;}
nP0rsd=0;if( strcmp(fit_RSD,"P0") ==0 || strcmp(fit_RSD,"P02") ==0 ||  strcmp(fit_RSD,"P04") ==0 || strcmp(fit_RSD,"P024")==0 ){nP0rsd=NeffP0rsd;}
nP2rsd=0;if( strcmp(fit_RSD,"P2") ==0 || strcmp(fit_RSD,"P02") ==0 ||  strcmp(fit_RSD,"P24") ==0 || strcmp(fit_RSD,"P024")==0 ){nP2rsd=NeffP2rsd;}
nP4rsd=0;if( strcmp(fit_RSD,"P4") ==0 || strcmp(fit_RSD,"P24") ==0 ||  strcmp(fit_RSD,"P24") ==0 || strcmp(fit_RSD,"P024")==0 ){nP4rsd=NeffP4rsd;}
nB0rsd=0;if( strcmp(do_bispectrum,"yes") == 0  && rsd == 1){nB0rsd=NeffB0rsd;}

if( strcmp(type_BAORSD,"FSalphasrecon") == 0 ){nP0bao=2;nP2bao=0;nP4bao=0;nB0bao=0;}

if(bao==1)
{
if( strcmp(fit_BAO,"P0") ==0 || strcmp(fit_BAO,"P02") ==0 ||  strcmp(fit_BAO,"P04") ==0 || strcmp(fit_BAO,"P024")==0 ){baoshift=baoshift+NeffP0bao;}
if( strcmp(fit_BAO,"P2") ==0 || strcmp(fit_BAO,"P02") ==0 ||  strcmp(fit_BAO,"P24") ==0 || strcmp(fit_BAO,"P024")==0 ){baoshift=baoshift+NeffP2bao;}
if( strcmp(fit_BAO,"P4") ==0 || strcmp(fit_BAO,"P24") ==0 ||  strcmp(fit_BAO,"P24") ==0 || strcmp(fit_BAO,"P024")==0 ){baoshift=baoshift+NeffP4bao;}
}
if(rsd==1)
{
if( strcmp(fit_RSD,"P0") ==0 || strcmp(fit_RSD,"P02") ==0 ||  strcmp(fit_RSD,"P04") ==0 || strcmp(fit_RSD,"P024")==0 ){rsdshift=rsdshift+NeffP0rsd;}
if( strcmp(fit_RSD,"P2") ==0 || strcmp(fit_RSD,"P02") ==0 ||  strcmp(fit_RSD,"P24") ==0 || strcmp(fit_RSD,"P024")==0 ){rsdshift=rsdshift+NeffP2rsd;}
if( strcmp(fit_RSD,"P4") ==0 || strcmp(fit_RSD,"P24") ==0 ||  strcmp(fit_RSD,"P24") ==0 || strcmp(fit_RSD,"P024")==0 ){rsdshift=rsdshift+NeffP4rsd;}
}

//printf("%d %d %d, %d %d %d\n",NeffP0bao,NeffP2bao,NeffP4bao, NeffP0rsd,NeffP2rsd,NeffP4rsd);
//printf("%d %d\n",baoshift,rsdshift);
//printf("%d %d %d, %d %d %d\n",nP0bao,nP2bao,nP4bao,nP0rsd,nP2rsd,nP4rsd);
//printf("Ncov=%d\n",Ncov);
//exit(0);
if(strcmp(do_power_spectrum,"yes") == 0){Pshift=baoshift+rsdshift;}//printf("Pshift=%d (%d,%d)\n",Pshift,baoshift,rsdshift);
scaling_factor=SCALING_COVARIANCE;

Variance = (double*) calloc( Ncov*Ncov, sizeof(double));
Pvectorav = (double*) calloc( Ncov, sizeof(double));
Pvector = (double**) calloc(Nrealizations,sizeof(double*));

for(j=0;j<Nrealizations;j++){Pvector[j] = (double*) calloc(Ncov,sizeof(double));}

if( strcmp(do_power_spectrum, "yes") == 0)
{

if(bao==1){

if(NeffP0bao>0){P0baoav=  (double*) calloc( NeffP0bao, sizeof(double));}
if(NeffP2bao>0){P2baoav=  (double*) calloc( NeffP2bao, sizeof(double));}
if(NeffP4bao>0){P4baoav=  (double*) calloc( NeffP4bao, sizeof(double));}

}

if(rsd==1){
if(NeffP0rsd>0){P0rsdav=  (double*) calloc( NeffP0rsd, sizeof(double));}
if(NeffP2rsd>0){P2rsdav=  (double*) calloc( NeffP2rsd, sizeof(double));}
if(NeffP4rsd>0){P4rsdav=  (double*) calloc( NeffP4rsd, sizeof(double));}

}


if(bao==1){
P0bao = (double**) calloc(Nrealizations,sizeof(double*));
P2bao = (double**) calloc(Nrealizations,sizeof(double*));
P4bao = (double**) calloc(Nrealizations,sizeof(double*));}

if(rsd==1){
P0rsd = (double**) calloc(Nrealizations,sizeof(double*));
P2rsd = (double**) calloc(Nrealizations,sizeof(double*));
P4rsd = (double**) calloc(Nrealizations,sizeof(double*));}


for(j=0;j<Nrealizations;j++)
{

if(bao==1){
if(NeffP0bao>0){   P0bao[j] = (double*) calloc(NeffP0bao,sizeof(double));}
if(NeffP2bao>0){   P2bao[j] = (double*) calloc(NeffP2bao,sizeof(double));}
if(NeffP4bao>0){   P4bao[j] = (double*) calloc(NeffP4bao,sizeof(double));}

}

if(rsd==1){
if(NeffP0rsd>0){   P0rsd[j] = (double*) calloc(NeffP0rsd,sizeof(double));}
if(NeffP2rsd>0){   P2rsd[j] = (double*) calloc(NeffP2rsd,sizeof(double));}
if(NeffP4rsd>0){   P4rsd[j] = (double*) calloc(NeffP4rsd,sizeof(double));}

}



}
//exit(0);

   Ncounter=0;
   for(j=1;j<=Nrealizations;j++)
   {

if(bao==1){
        sprintf(full_path_bao,"%s%.4d.txt",path_to_mocks_bao,j);
        sprintf(full_path2_bao,"%s%d.txt",path_to_mocks_bao,j);
        sprintf(full_path3_bao,"%s%.4d.dat",path_to_mocks_bao,j);
        sprintf(full_path4_bao,"%s%d.dat",path_to_mocks_bao,j);
        sprintf(full_path5_bao,"%s%.4d.cov",path_to_mocks_bao,j);
        sprintf(full_path6_bao,"%s%d.cov",path_to_mocks_bao,j);}

if(rsd==1){
        sprintf(full_path_rsd,"%s%.4d.txt",path_to_mocks_rsd,j);
        sprintf(full_path2_rsd,"%s%d.txt",path_to_mocks_rsd,j);
        sprintf(full_path3_rsd,"%s%.4d.dat",path_to_mocks_rsd,j);
        sprintf(full_path4_rsd,"%s%d.dat",path_to_mocks_rsd,j);}

f1bao=NULL;
f2bao=NULL;
f3bao=NULL;
f4bao=NULL;
f5bao=NULL;
f6bao=NULL;
if(bao==1){
        f1bao=fopen(full_path_bao,"r");
        f2bao=fopen(full_path2_bao,"r");
        f3bao=fopen(full_path3_bao,"r");
        f4bao=fopen(full_path4_bao,"r");
        f5bao=fopen(full_path5_bao,"r");
        f6bao=fopen(full_path6_bao,"r");}

f1rsd=NULL;
f2rsd=NULL;
f3rsd=NULL;
f4rsd=NULL;
if(rsd==1){
        f1rsd=fopen(full_path_rsd,"r");
        f2rsd=fopen(full_path2_rsd,"r");
        f3rsd=fopen(full_path3_rsd,"r");
        f4rsd=fopen(full_path4_rsd,"r");}

        fail=0;
        if(f1bao==NULL && f2bao==NULL && f3bao==NULL && f4bao==NULL && f5bao==NULL && f6bao==NULL && bao==1 && rsd==0 )
        {
          printf("Warning, realization %d missing at %s.\n",j,path_to_mocks_bao);fail=1;
        }
        if(f1rsd==NULL && f2rsd==NULL && f3rsd==NULL && f4rsd==NULL && bao==0 && rsd==1 )
        {
          printf("Warning, realization %d missing at %s.\n",j,path_to_mocks_rsd);fail=1;
        }

        if( f1bao==NULL && f2bao==NULL && f3bao==NULL && f4bao==NULL && f5bao==NULL && f6bao==NULL && f1rsd==NULL && f2rsd==NULL && f3rsd==NULL && f4rsd==NULL && bao==1 && rsd==1 )
        {
          printf("Warning, realization %d missing at %s or %s.\n",j,path_to_mocks_bao, path_to_mocks_rsd);fail=1;
        }

if(f1bao!=NULL && bao==1){fclose(f1bao);modebao=1;}
if(f2bao!=NULL && bao==1){fclose(f2bao);modebao=2;}
if(f3bao!=NULL && bao==1){fclose(f3bao);modebao=3;}
if(f4bao!=NULL && bao==1){fclose(f4bao);modebao=4;}
if(f5bao!=NULL && bao==1){fclose(f5bao);modebao=5;}
if(f6bao!=NULL && bao==1){fclose(f6bao);modebao=6;}
if(f1rsd!=NULL && rsd==1){fclose(f1rsd);modersd=1;}
if(f2rsd!=NULL && rsd==1){fclose(f2rsd);modersd=2;}
if(f3rsd!=NULL && rsd==1){fclose(f3rsd);modersd=3;}
if(f4rsd!=NULL && rsd==1){fclose(f4rsd);modersd=4;}

//printf("%d %d\n",modebao,modersd);
//printf("iterini=%d, iterfin=%d\n",iteration_ini,iteration_fin);
//exit(0);
        if(fail==0)
        {

//iterate one or two times
for(iteration=iteration_ini;iteration<=iteration_fin;iteration++){
  
           Nlines=0;
if(iteration==1){
           if(modebao==1){Nlines=countlines(full_path_bao);Nall=count_all_lines(full_path_bao);}
           if(modebao==2){Nlines=countlines(full_path2_bao);Nall=count_all_lines(full_path2_bao);}
           if(modebao==3){Nlines=countlines(full_path3_bao);Nall=count_all_lines(full_path3_bao);}
           if(modebao==4){Nlines=countlines(full_path4_bao);Nall=count_all_lines(full_path4_bao);}
           if(modebao==5){Nlines=countlines(full_path5_bao);Nall=count_all_lines(full_path5_bao);}
           if(modebao==6){Nlines=countlines(full_path6_bao);Nall=count_all_lines(full_path6_bao);}

if( strcmp(type_BAORSD,"FSalphasrecon") ==0 ){Nlines=2;Nall=3;}


}

if(iteration==2){
           if(modersd==1){Nlines=countlines(full_path_rsd);Nall=count_all_lines(full_path_rsd);}
           if(modersd==2){Nlines=countlines(full_path2_rsd);Nall=count_all_lines(full_path2_rsd);}
           if(modersd==3){Nlines=countlines(full_path3_rsd);Nall=count_all_lines(full_path3_rsd);}
           if(modersd==4){Nlines=countlines(full_path4_rsd);Nall=count_all_lines(full_path4_rsd);}}
        

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
if(modebao==5){f=fopen(full_path5_bao,"r");}
if(modebao==6){f=fopen(full_path6_bao,"r");}
}
if(iteration==2){
if(modersd==1){f=fopen(full_path_rsd,"r");}
if(modersd==2){f=fopen(full_path2_rsd,"r");}
if(modersd==3){f=fopen(full_path3_rsd,"r");}
if(modersd==4){f=fopen(full_path4_rsd,"r");}
}


for(i=0;i<Nall-Nlines;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line


if(iteration==1){kminP0=kminP0bao;kminP2=kminP2bao;kminP4=kminP4bao;kmaxP0=kmaxP0bao;kmaxP2=kmaxP2bao;kmaxP4=kmaxP4bao;}
if(iteration==2){kminP0=kminP0rsd;kminP2=kminP2rsd;kminP4=kminP4rsd;kmaxP0=kmaxP0rsd;kmaxP2=kmaxP2rsd;kmaxP4=kmaxP4rsd;}


              for(i=0;i<Nlines;i++)
              {

                  if(iteration==1 && strcmp(type_BAORSD,"FSalphasrecon") ==0 )
                  {
                  get_linealphas(f,params,1);
                  k=(kmaxP0+kminP0)*0.5;//always will pass kmin.kmax condition
                  p0=params[2];//printf("(iter=%d) Nreal=%d; Nline=%d alpha=%lf k=%lf (%lf,%lf)\n",iteration,j,i,p0,k,kminP0,kmaxP0);
                  }else{
                  get_lineP(f,params,1);
                  k=params[1];
                  p0=params[2];
                  p2=params[3];
                  p4=params[4];}

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

                  if( strcmp(type_BAORSD, "FS") == 0 || strcmp(type_BAORSD, "BAOANISO") == 0 || strcmp(type_BAORSD, "FSBAOANISO") == 0 || strcmp(type_BAORSD, "FSalphasrecon") == 0)
                  {

                      if(k>kminP2 && k<kmaxP2)
                      {
                         if(iteration==1 && strcmp(type_BAORSD, "FSalphasrecon") != 0){P2bao[Ncounter-1][i2]=p2;}
                         if(iteration==2){P2rsd[Ncounter-1][i2]=p2;}

                         if(iteration==1){
                         if( strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P24") == 0){Pvector[Ncounter-1][l2]=p2;}
                         if( strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P024") == 0){Pvector[Ncounter-1][l2+NeffP0bao]=p2;}

                         if(strcmp(type_BAORSD, "FSalphasrecon") != 0){i2++;}
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
                         if(iteration==1 && strcmp(type_BAORSD, "FSalphasrecon") != 0){P4bao[Ncounter-1][i4]=p4;}
                         if(iteration==2){P4rsd[Ncounter-1][i4]=p4;}

                         if(iteration==1){
	                 if( strcmp(fit_BAO, "P4") == 0){Pvector[Ncounter-1][l4]=p4;}
                         if( strcmp(fit_BAO, "P04") == 0){Pvector[Ncounter-1][l4+NeffP0bao]=p4;}
                         if( strcmp(fit_BAO, "P24") == 0){Pvector[Ncounter-1][l4+NeffP2bao]=p4;}
                         if( strcmp(fit_BAO, "P024") == 0){Pvector[Ncounter-1][l4+NeffP0bao+NeffP2bao]=p4;}
                         if(strcmp(type_BAORSD, "FSalphasrecon") != 0){i4++;}

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
                fclose(f);              
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

//printf("%d %d %d\n",i0bao,i2bao,i4bao);
//exit(0);
//Get the mean of P0,P2 and P4
for(j=0;j<Ncounter;j++)
{

if(bao==1){
    if(i0bao>0){
    for(i=0;i<i0bao;i++)
    {
        P0baoav[i]=P0baoav[i]+P0bao[j][i]/Ncounter*1.;
    }}
    if(i2bao>0){
    for(i=0;i<i2bao;i++)
    {
        P2baoav[i]=P2baoav[i]+P2bao[j][i]/Ncounter*1.;

    }}
    if(i4bao>0){
    for(i=0;i<i4bao;i++)
    {
        P4baoav[i]=P4baoav[i]+P4bao[j][i]/Ncounter*1.;

    }}
}
//exit(0);
if(rsd==1){
    if(i0rsd>0){
    for(i=0;i<i0rsd;i++)
    {
        P0rsdav[i]=P0rsdav[i]+P0rsd[j][i]/Ncounter*1.;
    }}
    if(i2rsd>0){
    for(i=0;i<i2rsd;i++)
    {
        P2rsdav[i]=P2rsdav[i]+P2rsd[j][i]/Ncounter*1.;

    }}
    if(i4rsd>0){
    for(i=0;i<i4rsd;i++)
    {
        P4rsdav[i]=P4rsdav[i]+P4rsd[j][i]/Ncounter*1.;

    }}
}

    for(i=0;i<l0bao+l2bao+l4bao+l0rsd+l2rsd+l4rsd;i++)
    {
        Pvectorav[i]=Pvectorav[i]+Pvector[j][i]/Ncounter*1.;
    }


}

if(Ncounter!=Nrealizations){printf("Warning, %d/%d realizations used for Power Spectrum.\n",Ncounter,Nrealizations);}

}//end of do power spectrum if

//exit(0);
if( strcmp(do_bispectrum, "yes") == 0)
{
modersd=-1;
modebao=-1;
Bbaoshift=0;if(bao==1 && rsd==1){Bbaoshift=NeffB0bao;}

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

if(bao==1){
        sprintf(full_path_bao,"%s%.4d.txt",path_to_mocks_bis_bao,j);
        sprintf(full_path2_bao,"%s%d.txt",path_to_mocks_bis_bao,j);
        sprintf(full_path3_bao,"%s%.4d.dat",path_to_mocks_bis_bao,j);
        sprintf(full_path4_bao,"%s%d.dat",path_to_mocks_bis_bao,j);
        sprintf(full_path5_bao,"%s%.4d.cov",path_to_mocks_bis_bao,j);
        sprintf(full_path6_bao,"%s%d.cov",path_to_mocks_bis_bao,j);
        }

if(rsd==1){
        sprintf(full_path_rsd,"%s%.4d.txt",path_to_mocks_bis_rsd,j);
        sprintf(full_path2_rsd,"%s%d.txt",path_to_mocks_bis_rsd,j);
        sprintf(full_path3_rsd,"%s%.4d.dat",path_to_mocks_bis_rsd,j);
        sprintf(full_path4_rsd,"%s%d.dat",path_to_mocks_bis_rsd,j);}

//avoids error on uninitialised variables
f1bao=NULL;
f2bao=NULL;
f3bao=NULL;        
f4bao=NULL;
f5bao=NULL;
f6bao=NULL;
if(bao==1){
        f1bao=fopen(full_path_bao,"r");
        f2bao=fopen(full_path2_bao,"r");
        f3bao=fopen(full_path3_bao,"r");
        f4bao=fopen(full_path4_bao,"r");
        f5bao=fopen(full_path5_bao,"r");        
        f6bao=fopen(full_path6_bao,"r");
        }

        f1rsd=NULL;
        f2rsd=NULL;
        f3rsd=NULL;
        f4rsd=NULL;
if(rsd==1){
        
        f1rsd=fopen(full_path_rsd,"r");
        f2rsd=fopen(full_path2_rsd,"r");
        f3rsd=fopen(full_path3_rsd,"r");
        f4rsd=fopen(full_path4_rsd,"r");}


    fail=0;
        if(f1bao==NULL && f2bao==NULL && f3bao==NULL && f4bao==NULL && f5bao==NULL && f6bao==NULL && bao==1 && rsd==0 )
        {
          printf("Warning, realization %d missing at %s.\n",j,path_to_mocks_bis_bao);fail=1;
        }
        if(f1rsd==NULL && f2rsd==NULL && f3rsd==NULL && f4rsd==NULL && bao==0 && rsd==1 )
        {
          printf("Warning, realization %d missing at %s.\n",j,path_to_mocks_bis_rsd);fail=1;
        }

        if( f1bao==NULL && f2bao==NULL && f3bao==NULL && f4bao==NULL && f5bao==NULL && f6bao==NULL && f1rsd==NULL && f2rsd==NULL && f3rsd==NULL && f4rsd==NULL && bao==1 && rsd==1 )
        {
          printf("Warning, realization %d missing at %s or %s.\n",j,path_to_mocks_bis_bao, path_to_mocks_bis_rsd);fail=1;
        }

if(f1bao!=NULL && bao==1){fclose(f1bao);modebao=1;}
if(f2bao!=NULL && bao==1){fclose(f2bao);modebao=2;}
if(f3bao!=NULL && bao==1){fclose(f3bao);modebao=3;}
if(f4bao!=NULL && bao==1){fclose(f4bao);modebao=4;}
if(f5bao!=NULL && bao==1){fclose(f5bao);modebao=4;}
if(f6bao!=NULL && bao==1){fclose(f6bao);modebao=4;}

if(f1rsd!=NULL && rsd==1){fclose(f1rsd);modersd=1;}
if(f2rsd!=NULL && rsd==1){fclose(f2rsd);modersd=2;}
if(f3rsd!=NULL && rsd==1){fclose(f3rsd);modersd=3;}
if(f4rsd!=NULL && rsd==1){fclose(f4rsd);modersd=4;}


if(f1bao==NULL && f2bao==NULL && f3bao==NULL && f4bao==NULL && f5bao==NULL && f6bao==NULL && bao==1){printf("Error with realization %d, exiting now...\n",j+1);exit(0);}
if(f1rsd==NULL && f2rsd==NULL && f3rsd==NULL && f4rsd==NULL && rsd==1){printf("Error with realization %d, exiting now...\n",j+1);exit(0);}


       if(fail==0)
       {

for(iteration=iteration_ini;iteration<=iteration_fin;iteration++){

           Nlines=0;
if(iteration==1){
           if(modebao==1){Nlines=countlines(full_path_bao);Nall=count_all_lines(full_path_bao);}
           if(modebao==2){Nlines=countlines(full_path2_bao);Nall=count_all_lines(full_path2_bao);}
           if(modebao==3){Nlines=countlines(full_path3_bao);Nall=count_all_lines(full_path3_bao);}
           if(modebao==4){Nlines=countlines(full_path4_bao);Nall=count_all_lines(full_path4_bao);}
           if(modebao==5){Nlines=countlines(full_path5_bao);Nall=count_all_lines(full_path5_bao);}
           if(modebao==6){Nlines=countlines(full_path6_bao);Nall=count_all_lines(full_path6_bao);}}

if(iteration==2){
           if(modersd==1){Nlines=countlines(full_path_rsd);Nall=count_all_lines(full_path_rsd);}
           if(modersd==2){Nlines=countlines(full_path2_rsd);Nall=count_all_lines(full_path2_rsd);}
           if(modersd==3){Nlines=countlines(full_path3_rsd);Nall=count_all_lines(full_path3_rsd);}
           if(modersd==4){Nlines=countlines(full_path4_rsd);Nall=count_all_lines(full_path4_rsd);}}


           if(Nlines!=0)
           {
              if(bao==1 && rsd==1){
              if(iteration==1){Ncounterbis++;}}//count them only once
              else{Ncounterbis++;}
if(Ncounterbis>Nrealizations){printf("Error Ncounter>Nreal for bispectrum. Exiting now...\n");exit(0);}

              i0bis=0;
if(iteration==1){
if(modebao==1){f=NULL;f=fopen(full_path_bao,"r");}
if(modebao==2){f=NULL;f=fopen(full_path2_bao,"r");}
if(modebao==3){f=NULL;f=fopen(full_path3_bao,"r");}
if(modebao==4){f=NULL;f=fopen(full_path4_bao,"r");}
if(modebao==5){f=NULL;f=fopen(full_path5_bao,"r");}
if(modebao==6){f=NULL;f=fopen(full_path6_bao,"r");}
}
if(iteration==2){
if(modersd==1){f=NULL;f=fopen(full_path_rsd,"r");}
if(modersd==2){f=NULL;f=fopen(full_path2_rsd,"r");}
if(modersd==3){f=NULL;f=fopen(full_path3_rsd,"r");}
if(modersd==4){f=NULL;f=fopen(full_path4_rsd,"r");}
}

for(i=0;i<Nall-Nlines;i++){
do{ch = fgetc(f);}while(ch != '\n');
}//skip header line


if(iteration==1){kminB0=kminB0bao;kmaxB0=kmaxB0bao;}
if(iteration==2){kminB0=kminB0rsd;kmaxB0=kmaxB0rsd;}

for(i=0;i<Nlines;i++)
              {
                    get_lineB(f,params,1);
                    k1=params[1];
                    k2=params[3];
                    k3=params[5];
                    b0=params[6];
                    q0=params[8];

if( strcmp(bispectrumBQ,"B") == 0){bispectrum=b0;}
if( strcmp(bispectrumBQ,"Q") == 0){
//printf("Sorry no reduced bispectrum. Exiting now...\n");exit(0);//comment line if Q available
bispectrum=q0;
}

                  if(k1>kminB0 && k1<kmaxB0 && k2>kminB0 && k2<kmaxB0 && k3>kminB0 && k3<kmaxB0)
                  {
                      if(iteration==1){
                      if(Ncounterbis-1>=Nrealizations || Ncounterbis-1<0 || i0bis>=NeffB0bao || i0bis<0){printf("Error within get_cov_from_mocks (B0bao). Exiting now...\n");exit(0);}
                                     B0bao[Ncounterbis-1][i0bis]=bispectrum;
                      }//bao
                      if(iteration==2){
                      if(Ncounterbis-1>=Nrealizations || Ncounterbis-1<0 || i0bis>=NeffB0rsd || i0bis<0){printf("Error within get_cov_from_mocks (B0rsd). Exiting now...\n");exit(0);}
                      B0rsd[Ncounterbis-1][i0bis]=bispectrum;
                      }//rsd



                      if(iteration==1){
                      if(Ncounterbis-1>=Nrealizations || Ncounterbis-1<0 || i0bis+Pshift>=Ncov || i0bis+Pshift<0){printf("Error within get_cov_from_mocks (Pvector (bao)). Exiting now...\n");exit(0);}
                      Pvector[Ncounterbis-1][i0bis+Pshift]=bispectrum;
                      i0bis++;                   
                      }
                      if(iteration==2){
                      if(Ncounterbis-1>=Nrealizations || Ncounterbis-1<0 || i0bis+Bbaoshift+Pshift>=Ncov || i0bis+Bbaoshift+Pshift<0){printf("Error within get_cov_from_mocks (Pvector (rsd)). Exiting now...\n");exit(0);}   
                      Pvector[Ncounterbis-1][i0bis+Bbaoshift+Pshift]=bispectrum;
                      i0bis++;                      
                      }

                   }

}//end of i-for
fclose(f);

if(iteration==1){i0bisbao=i0bis;}
if(iteration==2){i0bisrsd=i0bis;}

           }//lines !=0
           else
           {
               printf("Warning, realization %d with no data at %s or %s\n",j,path_to_mocks_bao,path_to_mocks_rsd);
           }

}//iteration=1,2 loop

        }//if fail=0


   }//end of j-for


for(j=0;j<Ncounterbis;j++)
{

if(bao==1){
    for(i=0;i<i0bisbao;i++)
    {
        B0baoav[i]=B0baoav[i]+B0bao[j][i]/Ncounterbis*1.;
    }
}

if(rsd==1){
    for(i=0;i<i0bisrsd;i++)
    {
        B0rsdav[i]=B0rsdav[i]+B0rsd[j][i]/Ncounterbis*1.;
    }
}

    for(i=l0bao+l2bao+l4bao+l0rsd+l2rsd+l4rsd;i<l0bao+l2bao+l4bao+l0rsd+l2rsd+l4rsd+i0bisbao+i0bisrsd;i++)
    {
        Pvectorav[i]=Pvectorav[i]+Pvector[j][i]/Ncounterbis*1.;
    }

}


if(Ncounterbis!=Nrealizations){printf("Warning, %d/%d realizations used for Bispectrum.\n",Ncounterbis,Nrealizations);}

}//end of bispectrum if

if( strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 && Ncounter!=Ncounterbis){printf("Error, number of used realizations is different for the power spectrum and bispectrum: %d != %d. Exiting now...\n",Ncounter,Ncounterbis);exit(0);} 


//Build errors

if( strcmp(do_power_spectrum, "yes") == 0)
{

if( bao==1 && strcmp(type_BAORSD, "FSalphasrecon") != 0){

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
        errP0bao[i]=sqrt(errP0bao[i]/(Ncounter-1.))/sqrt(scaling_factor);
    }
    for(i=0;i<i2bao;i++)
    {
        errP2bao[i]=sqrt(errP2bao[i]/(Ncounter-1.))/sqrt(scaling_factor);

    }
    for(i=0;i<i4bao;i++)
    {
        errP4bao[i]=sqrt(errP4bao[i]/(Ncounter-1.))/sqrt(scaling_factor);
    }
}

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
        errP0rsd[i]=sqrt(errP0rsd[i]/(Ncounter-1.))/sqrt(scaling_factor);
    }
    for(i=0;i<i2rsd;i++)
    {
        errP2rsd[i]=sqrt(errP2rsd[i]/(Ncounter-1.))/sqrt(scaling_factor);

    }
    for(i=0;i<i4rsd;i++)
    {
        errP4rsd[i]=sqrt(errP4rsd[i]/(Ncounter-1.))/sqrt(scaling_factor);
    }


}

}
//exit(0);

if( strcmp(do_bispectrum, "yes") == 0)
{

if( bao==1 && strcmp(type_BAORSD, "FSalphasrecon") != 0){
  for(i_real=0;i_real<Ncounterbis;i_real++)
  {
    for(i=0;i<i0bisbao;i++)
    {
        errB0bao[i]=errB0bao[i]+pow( B0baoav[i]-B0bao[i_real][i],2);
    }

  }

    for(i=0;i<i0bisbao;i++)
    {
        errB0bao[i]=sqrt(errB0bao[i]/(Ncounterbis-1.))/sqrt(scaling_factor);
    }
}

if( rsd==1){
  for(i_real=0;i_real<Ncounterbis;i_real++)
  {
    for(i=0;i<i0bisrsd;i++)
    {
        errB0rsd[i]=errB0rsd[i]+pow( B0rsdav[i]-B0rsd[i_real][i],2);
    }

  }

    for(i=0;i<i0bisrsd;i++)
    {
        errB0rsd[i]=sqrt(errB0rsd[i]/(Ncounterbis-1.))/sqrt(scaling_factor);
    }
}

}

//Build covariance
if(strcmp(do_bispectrum,"yes") == 0){Ncounter=Ncounterbis;}


for(i_real=0;i_real<Ncounter;i_real++)
{
for(j=0;j<Ncov;j++)
{
jelement=Pvectorav[j]-Pvector[i_real][j];
//printf("real=%d, P[%d]=%lf\n",i_real,j,Pvector[i_real][j]);
//if(i_real==3){exit(0);}
for(i=0;i<Ncov;i++)
{
ielement=Pvectorav[i]-Pvector[i_real][i];
Variance[j+i*Ncov]=Variance[j+i*Ncov]+ielement*jelement/(Ncounter*1.-1.)*(1./scaling_factor);
}
}
}

//correct variance here by FSa
if( strcmp(covarianceFSa_option,"varying") == 0 )
{
diag_apara=Variance[0+0*Ncov];
diag_aperp=Variance[1+1*Ncov];
for(j=0;j<Ncov;j++)
{
for(i=0;i<Ncov;i++)
{

if(i==0)
{
   if(j==0){Variance[j+i*Ncov]=paramsBAO[0];}
   if(j==1){Variance[j+i*Ncov]=sqrt(paramsBAO[0]*paramsBAO[1])*paramsBAO[2];}
   if(j>1){Variance[j+i*Ncov]=Variance[j+i*Ncov]*sqrt(paramsBAO[0]/diag_apara);}

}

if(i==1)
{
   if(j==0){Variance[j+i*Ncov]=sqrt(paramsBAO[0]*paramsBAO[1])*paramsBAO[2];}
   if(j==1){Variance[j+i*Ncov]=paramsBAO[1];}
   if(j>1){Variance[j+i*Ncov]=Variance[j+i*Ncov]*sqrt(paramsBAO[1]/diag_aperp);}

}

if(i>1)
{
   if(j==0){Variance[j+i*Ncov]=Variance[j+i*Ncov]*sqrt(paramsBAO[0]/diag_apara);}
   if(j==1){Variance[j+i*Ncov]=Variance[j+i*Ncov]*sqrt(paramsBAO[1]/diag_aperp);}
//   if(j>1){Variance[j+i*Ncov]=;}//do nothing, auto k-part not modified

}

}
}

}


sprintf(covfile_in,"%s.cov",covfile);
write_cov_file(covfile_in,Variance,Ncov,Nrealizations,nP0bao,nP0rsd,nP2bao,nP2rsd,nP4bao,nP4rsd,nB0bao,nB0rsd);

//Apply Compression here
if( matrix_compression !=NULL)
{

Variance_comp= (double*) calloc( output*output, sizeof(double));
if(Ncov != input){printf("Error, input from compressed covariance (%d) does not correspond to the size of the uncompressed covariance (%d). Exiting  now...\n",Ncov,input);exit(0);}

  for(i=0;i<output;i++)//filas compression matrix
  {
     for(j=0;j<output;j++)//columnas compression matrix
     {

         Variance_comp[j+i*output]=0;      
         for(ii=0;ii<Ncov;ii++)//filas original matrix
         {
               for(jj=0;jj<Ncov;jj++)//columnas original matirx & columnas compression matrix
               {
                   Variance_comp[j+i*output]=Variance_comp[j+i*output]+(matrix_compression[i][jj]*Variance[jj+ii*Ncov]*matrix_compression[j][ii]);//check this
               }   
        }

     }

  }

//invert compressed covariance

        m = gsl_matrix_alloc (output, output);
        inverse = gsl_matrix_alloc (output, output);
        identity = gsl_matrix_alloc (output, output);

        perm = gsl_permutation_alloc (output);
         for(i=0;i<output;i++)
         {
         for(j=0;j<output;j++)
         {
                if(Variance_comp[i+i*output]==0 || Variance_comp[j+j*output]==0){printf("Error 0-element(s) in the diagonal of compressed covarince (V[%d]=%lf,V[%d]=%lf): Exiting now...\n",i,Variance_comp[i+i*output],j,Variance_comp[j+j*output]);exit(0);}
                gsl_matrix_set (m, i, j, Variance_comp[j+i*output]);
         }
         }

           free(Variance_comp);

           gsl_linalg_LU_decomp (m, perm, &s);
           gsl_linalg_LU_invert (m, perm, inverse);

         for(ii=0;ii<Ncov;ii++)//filas
         {
         for(jj=0;jj<Ncov;jj++)//col
         {
                 //here cov gets again dimensions of the uncompressed data-vector
          //     cov[j+Ncov*i]=1./((1.-(Ncov+1.)/(Ncounter*1.-1.))*gsl_matrix_get (inverse, i, j));
                cov[jj+Ncov*ii]=0;
                for(i=0;i<output;i++)//filas
                {
                for(j=0;j<output;j++)//col
                {
                 
                   cov[jj+ii*Ncov]=cov[jj+ii*Ncov]+(matrix_compression[i][jj]*gsl_matrix_get (inverse, i, j)*matrix_compression[j][ii]);

                }
                }

         }
         }


if( strcmp(covariance_correction,"none") == 0 || strcmp(covariance_correction,"Sellentin-Heavens") == 0)
{
hartlap_factor = 1.;
}
else{//Hartlap
hartlap_factor = 1./(1.-(output+1.)/(Ncounter*1.-1.));
}
out=Ncounter;

         for(ii=0;ii<Ncov;ii++)//filas
         {
         for(jj=0;jj<Ncov;jj++)//col
         {

//              cov[jj+Ncov*ii]=1./((1.-(output+1.)/(Ncounter*1.-1.))*cov[jj+Ncov*ii]);//invert it  and apply Hartlap factor
                cov[jj+Ncov*ii]=1./(cov[jj+Ncov*ii]/hartlap_factor);//invert it  and apply Hartlap factor
         }
         }


if( strcmp(covariance_correction,"Hartlap") == 0)
{
         printf("Hartlap correction of compressed covariance is %lf, %d x %d elements ( %d for P0-BAO, %d for P2-BAO, %d for P4-BAO, %d for B0-BAO, %d for P0-RSD, %d for P2-RSD, %d for P4-RSD, %d for B0-RSD; %d realizations)\n",1./((1.-(output+1.)/(Ncounter*1.-1.))),output,output,l0bao,l2bao,l4bao,i0bisbao,l0rsd,l2rsd,l4rsd,i0bisrsd,Ncounter);
}

         if(output>=Ncounter){printf("Error, not enough realizations to perform a reliable estimate of the covariance. Please reduce the dimension of your covariance or add more mock realizations. Exiting now...\n");exit(0);}

//inverse
//sprintf(covfile_in,"%s.invcov",covfile);
write_cov_file(covfile_in,cov,Ncov,Nrealizations,nP0bao,nP0rsd,nP2bao,nP2rsd,nP4bao,nP4rsd,nB0bao,nB0rsd);

}
else{

       //invert covariance

        m = gsl_matrix_alloc (Ncov, Ncov);
        inverse = gsl_matrix_alloc (Ncov, Ncov);
        identity = gsl_matrix_alloc (Ncov, Ncov);

        perm = gsl_permutation_alloc (Ncov);
         for(i=0;i<Ncov;i++)
         {
         for(j=0;j<Ncov;j++)
         {
                if(Variance[i+i*Ncov]==0 || Variance[j+j*Ncov]==0){printf("Error 0-element(s) in the diagonal of covarince (V[%d]=%lf,V[%d]=%lf): Exiting now...\n",i,Variance[i+i*Ncov],j,Variance[j+j*Ncov]);exit(0);}
                gsl_matrix_set (m, i, j, Variance[j+i*Ncov]);
         }
         }
           gsl_linalg_LU_decomp (m, perm, &s);
           gsl_linalg_LU_invert (m, perm, inverse);

if( strcmp(covariance_correction,"none") == 0 || strcmp(covariance_correction,"Sellentin-Heavens") == 0)
{
hartlap_factor = 1.;
}
else{//Hartlap
hartlap_factor = 1./(1.-(Ncov+1.)/(Ncounter*1.-1.));
}
out=Ncounter;

         for(i=0;i<Ncov;i++)
         {
         for(j=0;j<Ncov;j++)
         {
//           cov[j+Ncov*i]=1./((1.-(Ncov+1.)/(Ncounter*1.-1.))*gsl_matrix_get (inverse, i, j));
             cov[j+Ncov*i]=1./(gsl_matrix_get (inverse, i, j)/hartlap_factor);
         }
         }

if( strcmp(covariance_correction,"Hartlap") == 0)
{
         printf("Hartlap correction of covariance is %lf, %d x %d elements ( %d for P0-BAO, %d for P2-BAO, %d for P4-BAO, %d for B0-BAO, %d for P0-RSD, %d for P2-RSD, %d for P4-RSD, %d for B0-RSD; %d realizations)\n",1./((1.-(Ncov+1.)/(Ncounter*1.-1.))),Ncov,Ncov,l0bao,l2bao,l4bao,i0bisbao,l0rsd,l2rsd,l4rsd,i0bisrsd,Ncounter);
}

         if(Ncov>=Ncounter){printf("Error, not enough realizations to perform a reliable estimate of the covariance. Please reduce the dimension of your covariance or add more mock realizations. Exiting now...\n");exit(0);}

//inverse
sprintf(covfile_in,"%s.invcov",covfile);
write_cov_file(covfile_in,cov,Ncov,Nrealizations,nP0bao,nP0rsd,nP2bao,nP2rsd,nP4bao,nP4rsd,nB0bao,nB0rsd);

}//else NULL

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
if(NeffP0bao>0){free(P0baoav);}
if(NeffP2bao>0){free(P2baoav);}
if(NeffP4bao>0){free(P4baoav);}

if(NeffP0bao>0){freeTokens(P0bao,Nrealizations);}
if(NeffP2bao>0){freeTokens(P2bao,Nrealizations);}
if(NeffP4bao>0){freeTokens(P4bao,Nrealizations);}

}

if(rsd==1){
if(NeffP0rsd>0){free(P0rsdav);}
if(NeffP2rsd>0){free(P2rsdav);}
if(NeffP4rsd>0){free(P4rsdav);}
if(NeffP0rsd>0){freeTokens(P0rsd,Nrealizations);}
if(NeffP2rsd>0){freeTokens(P2rsd,Nrealizations);}
if(NeffP4rsd>0){freeTokens(P4rsd,Nrealizations);}

}


}


if( strcmp(do_bispectrum, "yes") == 0)
{
if(bao==1){
free(B0baoav);
freeTokens(B0bao,Nrealizations);
}
if(rsd==1){
free(B0rsdav);
freeTokens(B0rsd,Nrealizations);
}


}


return out;
}//end of get_cov

