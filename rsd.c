#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bispectrum.h"
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
#include "bao.h"
#include "cubature.h"

#define Pi (4.*atan(1.))

void do_window_matrix_for_P(double **matrix_window, char *type_of_analysis, int modeP0,int modeP2,int modeP4, double *k_theory,int N_Plin, int NeffP0, int NeffP2, int NeffP4,int factor_for_sampling, double *pos, double *W0,double *W2,double *W4,double *W6,double *W8,int Nmask,char *spacing_mask, double *k0, double *k2, double *k4,double k0min, double k0max, double k2min, double k2max, double k4min, double k4max, char *spacing_data,int Ncov,char *path_output,char *identifier, int ngcsgc, int baorsd)
{
fftw_plan plan1,plan2;
fftw_complex *a_pointer;a_pointer=NULL;
fftw_complex *b_pointer;b_pointer=NULL;
double kpad_up,kpad_down;
double params[5];
int params_int[2];
int Nplan;
double *k_theo;
int factor;
double *sum_matrix_in,*max_matrix_in;
double *P0,*P2,*P4;
double *P_theo0,*P_theo2,*P_theo4;
double pin,m;
int i_out,i_matrix,mode1,vector_input,N_in;
int modeP0_trial,modeP2_trial,modeP4_trial;
int i,j,c;
int j_k;
double *P_theo0win, *P_theo2win, *P_theo4win;
double ptheo,ptheo2,ptheo4;
double XMIN1[1]={-1};
double XMAX1[1]={+1};
int Nelements;
double kmax,kmin;
int NeffP;
int NeffP_max;
double div;
char bao_rsd[200];
if( baorsd == 0){sprintf(bao_rsd,"FS");}
if( baorsd == 1){sprintf(bao_rsd,"BAO");}

char out_Pin[2000];sprintf(out_Pin,"%s/Pin_Window_%s_%s_%d.txt",path_output,bao_rsd,identifier,ngcsgc);
char out_Pout[2000];sprintf(out_Pout,"%s/Pout_Window_%s_%s_%d.txt",path_output,bao_rsd,identifier,ngcsgc);
char out_M[2000];sprintf(out_M,"%s/Matrix_in_Window_%s_%s_%d.txt",path_output,bao_rsd,identifier,ngcsgc);
char out_M2[2000];sprintf(out_M2,"%s/Matrix_offdiag_in_Window_%s_%s_%d.txt",path_output,bao_rsd,identifier,ngcsgc);
FILE *fpin,*fpout,*fmatrix,*fmatrix2;
fpin=fopen(out_Pin,"w");
fpout=fopen(out_Pout,"w");
fmatrix=fopen(out_M,"w");
fmatrix2=fopen(out_M2,"w");

if(strcmp(spacing_data,"irregular") == 0){factor=1;}
else{factor=factor_for_sampling;}

get_kmin_kmax(params,modeP0,modeP2,modeP4,k0min,k0max,k2min,k2max,k4min,k4max);
kmax=params[1];
kmin=params[0];

get_NeffP(params_int,NeffP0,NeffP2,NeffP4,spacing_data,kmin,kmax,k0min,k0max,k2min,k2max,k4min,k4max, N_Plin );
NeffP_max=params_int[1];
NeffP=params_int[0]*factor;

NeffP0=NeffP0*modeP0;
NeffP2=NeffP2*modeP2;
NeffP4=NeffP4*modeP4;

Nelements=NeffP;
N_in=Nelements*(modeP0+modeP2+modeP4);

if(strcmp(spacing_data,"irregular") == 0)
{
k_theo = (double*) calloc( N_Plin, sizeof(double));
}
else
{
k_theo = (double*) calloc( NeffP, sizeof(double));
}

    for(j_k=0;j_k<Nelements;j_k++){k_theo[j_k]=get_ktheo(spacing_data,j_k,kmin,kmax,NeffP_max,factor_for_sampling,k_theory,NULL);}


set_mask_params(params,0.0,1.0,0.0,1.0,0.0);
Nplan=(int)(params[2]);
kpad_up=params[0];
kpad_down=params[1];

 plan1 = fftw_plan_dft_1d(Nplan,  a_pointer,  b_pointer,  -1, FFTW_ESTIMATE);//forward plan
 plan2 = fftw_plan_dft_1d(Nplan,  b_pointer,  b_pointer, +1, FFTW_ESTIMATE);//reverse plan

//matrix_window[N_in][N_out];
sum_matrix_in = (double*) calloc( Ncov, sizeof(double));
max_matrix_in = (double*) calloc( Ncov, sizeof(double));

//matrix_window = (double **) calloc(N_in, sizeof(double*));
//                 for(i=0;i<N_in;i++){matrix_window[i]= (double *) calloc(Ncov, sizeof(double));}

if(modeP0==1){P_theo0 = (double*) calloc( Nelements, sizeof(double));}
if(modeP2==1){P_theo2 = (double*) calloc( Nelements, sizeof(double));}
if(modeP4==1){P_theo4 = (double*) calloc( Nelements, sizeof(double));}

i_matrix=0;
for(mode1=0;mode1<=4;mode1=mode1+2){

    for(j_k=0;j_k<Nelements;j_k++)
    {

//printf("%d %d/%d\n",mode1,j_k,Nelements);
//k_theo[j_k]=get_ktheo(spacing_data,j_k,kmin,kmax,NeffP_max,factor_for_sampling,k_theory,NULL);

//if(strcmp(spacing_data,"linear") == 0){k_theo[j_k]=(j_k+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);}

//if(strcmp(spacing_data,"log") == 0){k_theo[j_k]=exp(  (j_k+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );}

//if(strcmp(spacing_data,"log10") == 0){k_theo[j_k]=pow(10,  (j_k+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );}

//if(strcmp(spacing_data,"irregular") == 0){k_theo[j_k]=k_theory[j_k];}

vector_input=0;

if(modeP0==1 && mode1==0)
{
         ptheo=10000;
           P_theo0[j_k]=ptheo;pin=ptheo;vector_input++;
             fprintf(fpin,"%d 0 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
             
}

if(modeP2==1 && mode1==2)
{
         ptheo2=10000;
         P_theo2[j_k]=ptheo2;pin=ptheo2;vector_input++;
         fprintf(fpin,"%d 2 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
        
}

if(modeP4==1 && mode1==4)
{
         ptheo4=10000;
         P_theo4[j_k]=ptheo4;pin=ptheo4;vector_input++;
         fprintf(fpin,"%d 4 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
}

if(vector_input>0){//some-vector-filled
if(modeP0==1){P_theo0win = (double*) calloc( Nelements, sizeof(double));}
if(modeP2==1){P_theo2win = (double*) calloc( Nelements, sizeof(double));}
if(modeP4==1){P_theo4win = (double*) calloc( Nelements, sizeof(double));}

if(modeP0==1 && mode1==0){
P_theo0win[j_k]=P_theo0[j_k];P_theo0[j_k]=0;
}

if(modeP2==1 && mode1==2){
P_theo2win[j_k]=P_theo2[j_k];P_theo2[j_k]=0;
}
if(modeP4==1 && mode1==4){
P_theo4win[j_k]=P_theo4[j_k];P_theo4[j_k]=0;
}

//printf("%lf %lf %d\n",k_theo[0], k_theo[NeffP-1],NeffP);
apply_mask(type_of_analysis, modeP0, modeP2, modeP4, k_theo,P_theo0win, P_theo2win, P_theo4win,NeffP, pos, W0, W2, W4, W6,W8,Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,0,baorsd);


if(modeP0==1){P0 = (double*) calloc( NeffP0, sizeof(double));}
if(modeP2==1){P2 = (double*) calloc( NeffP2, sizeof(double));}
if(modeP4==1){P4 = (double*) calloc( NeffP4, sizeof(double));}

j=0;
modeP0_trial=modeP0;
modeP2_trial=modeP2;
modeP4_trial=modeP4;

for(i_out=0;i_out<Ncov;i_out++)
{

       if(modeP4_trial==1 && modeP0_trial==0 && modeP2_trial==0){

            if(i_matrix==1){fprintf(fpout,"%d 4 %lf\n",i_out,k4[j]);}

       P4[j]= P_interpol(k4[j], k_theo, P_theo4win, Nelements);m=P4[j]/pin;
       
        j++;
        if(j==NeffP4){j=0;modeP4_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

        if(pin==0 || m==0){
if(k_theo[j_k]>kpad_down && k_theo[j_k]<kpad_up){

printf("Error generating the Window Matrix (P4, kout=%lf, kin[%d,%d]=%lf): matrix[%d][%d]= %lf / %lf\n",k4[j],j_k,mode1,k_theo[j_k],i_matrix-1,i_out,P4[j],pin);

//printf("Apply-mask parameters: %s (%d %d %d) %d %d %s (%lf,%lf) (%lf %lf) (%lf %lf) (%lf %lf) %s 0 %d \n",type_of_analysis, modeP0, modeP2, modeP4, NeffP, Nmask,spacing_mask, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,baorsd);

for(c=0;c<Nelements;c++){printf("outputed window power: k=%lf ",k_theo[c]);
if(modeP0==1){printf(" P0=%lf ",P_theo0win[c]);}
if(modeP2==1){printf(" P2=%lf ",P_theo2win[c]);}
if(modeP4==1){printf(" P4=%lf ",P_theo4win[c]);}
printf("\n");
}

exit(0);
}
else{
printf("Warning element kout[%d,4]=%lf, kin[%d,%d]=%lf in padding kpad_low=%lf, kpad_high=%lf\n",j,k4[j],j_k,mode1,k_theo[j_k],kpad_down,kpad_up);
}
}
        }

        if(modeP2_trial==1 && modeP0_trial==0){

            if(i_matrix==1){fprintf(fpout,"%d 2 %lf\n",i_out,k2[j]);}

        P2[j]= P_interpol(k2[j], k_theo, P_theo2win, Nelements);m=P2[j]/pin;
      
        j++;
        if(j==NeffP2){j=0;modeP2_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

if(pin==0 || m==0){

if(k_theo[j_k]>kpad_down && k_theo[j_k]<kpad_up){

printf("Error generating the Window Matrix (P2, kout=%lf, kin[%d,%d]=%lf): matrix[%d][%d]= %lf / %lf\n",k2[j],j_k,mode1,k_theo[j_k],i_matrix-1,i_out,P2[j],pin);
//printf("Apply-mask parameters: %s (%d %d %d) %d %d %s (%lf,%lf) (%lf %lf) (%lf %lf) (%lf %lf) %s 0 %d \n",type_of_analysis, modeP0, modeP2, modeP4, NeffP, Nmask,spacing_mask, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,baorsd);
for(c=0;c<Nelements;c++){printf("outputed window power: k=%lf ",k_theo[c]);
if(modeP0==1){printf(" P0=%lf ",P_theo0win[c]);}
if(modeP2==1){printf(" P2=%lf ",P_theo2win[c]);}
if(modeP4==1){printf(" P4=%lf ",P_theo4win[c]);}
printf("\n");
}
exit(0);
}
else{
printf("Warning element kout[%d,2]=%lf, kin[%d,%d]=%lf in padding kpad_low=%lf, kpad_high=%lf\n",j,k2[j],j_k,mode1,k_theo[j_k],kpad_down,kpad_up);
}
}


        }

        if(modeP0_trial==1){
        if(i_matrix==1){fprintf(fpout,"%d 0 %lf\n",i_out,k0[j]);}

        P0[j]= P_interpol(k0[j], k_theo, P_theo0win, Nelements);m=P0[j]/pin;
        
        j++;
        if(j==NeffP0){j=0;modeP0_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

        if(pin==0 || m==0){

if(k_theo[j_k]>kpad_down && k_theo[j_k]<kpad_up){

printf("Error generating the Window Matrix (P0, kout=%lf, kin[%d,%d]=%lf): matrix[%d][%d]= %lf / %lf\n",k0[j],j_k,mode1,k_theo[j_k],i_matrix-1,i_out,P0[j],pin);

//printf("Apply-mask parameters: %s (%d %d %d) %d %d %s (%lf,%lf) (%lf %lf) (%lf %lf) (%lf %lf) %s 0 %d \n",type_of_analysis, modeP0, modeP2, modeP4, NeffP, Nmask,spacing_mask, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,baorsd);

for(c=0;c<Nelements;c++){printf("outputed window power: k=%lf ",k_theo[c]);
if(modeP0==1){printf(" P0=%lf ",P_theo0win[c]);}
if(modeP2==1){printf(" P2=%lf ",P_theo2win[c]);}
if(modeP4==1){printf(" P4=%lf ",P_theo4win[c]);}
printf("\n");
}
exit(0);
}
else{
printf("Warning element kout[%d,0]=%lf, kin[%d,%d]=%lf in padding kpad_low=%lf, kpad_high=%lf\n",j,k0[j],j_k,mode1,k_theo[j_k],kpad_down,kpad_up);
}
}

        }

}//i_out loop

if(modeP0==1){free(P0);free(P_theo0win);}
if(modeP2==1){free(P2);free(P_theo2win);}
if(modeP4==1){free(P4);free(P_theo4win);}

}//if-some-vector-filled

}//j_k-Nelements
}//mode1

//printf("%d %d\n",i_matrix,Ncov);
for(i=0;i<i_matrix;i++)
{
    for(j=0;j<Ncov;j++)
    {
       fprintf(fmatrix,"%e\t",matrix_window[i][j]);
//    div=sum_matrix_in[j];//normalize by the sum of all elements of the line
       div=max_matrix_in[j];//normalize by the maximum value of the element of the line
        fprintf(fmatrix2,"%e\t",matrix_window[i][j]/div );
       if(j==Ncov-1){fprintf(fmatrix,"\n");}
       if(j==Ncov-1){fprintf(fmatrix2,"\n");}
    }

}

//freeTokens(matrix_window,N_in);
fftw_destroy_plan(plan1);
fftw_destroy_plan(plan2);
free(k_theo);
free(sum_matrix_in);
free(max_matrix_in);
fclose(fpin);
fclose(fpout);
fclose(fmatrix);
fclose(fmatrix2);
if(modeP0==1){free(P_theo0);}
if(modeP2==1){free(P_theo2);}
if(modeP4==1){free(P_theo4);}
}

/*
void do_Ptheo_RSD_window_wo_model(char *type_of_analysis, int modeP0,int modeP2,int modeP4,double **Theory,int N_Plin,  double k_theo[],int NeffP0, int NeffP2, int NeffP4,int factor_for_sampling, double *pos, double *W0,double *W2,double *W4,double *W6,double *W8,int Nmask,char *spacing_mask, double *k0, double *k2, double *k4, char *path_to_mask1,fftw_plan plan1,fftw_plan plan2,double k0min, double k0max, double k2min, double k2max, double k4min, double k4max, char *spacing_data,int Ncov,char *path_output,char *identifier, int ngcsgc)//delete
{
double **matrix_window;//Nelements*(modeP0+modeP2+modep4) x Ncov
double *sum_matrix_in,*max_matrix_in;
double *P0,*P2,*P4;
double *P_theo0,*P_theo2,*P_theo4;
double pin,m;
int i_out,i_matrix,mode1,vector_input,N_in;
int modeP0_trial,modeP2_trial,modeP4_trial;
int i,j;
int j_k;
double *P_theo0win, *P_theo2win, *P_theo4win;
double ptheo,ptheo2,ptheo4;
double XMIN1[1]={-1};
double XMAX1[1]={+1};
int Nelements;
double kmax,kmin;
int NeffP;
int NeffP_max;
double div;

char out_Pin[2000];sprintf(out_Pin,"%s/Pin_Window_%s_%d.txt",path_output,identifier,ngcsgc);
char out_Pout[2000];sprintf(out_Pout,"%s/Pout_Window_%s_%d.txt",path_output,identifier,ngcsgc);
char out_M[2000];sprintf(out_M,"%s/Matrix_in_Window_%s_%d.txt",path_output,identifier,ngcsgc);
char out_M2[2000];sprintf(out_M2,"%s/Matrix_offdiag_in_Window_%s_%d.txt",path_output,identifier,ngcsgc);
FILE *fpin,*fpout,*fmatrix,*fmatrix2;
fpin=fopen(out_Pin,"w");
fpout=fopen(out_Pout,"w");
fmatrix=fopen(out_M,"w");
fmatrix2=fopen(out_M2,"w");


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

Nelements=NeffP;
N_in=Nelements*(modeP0+modeP2+modeP4);
//matrix_window[N_in][N_out];
sum_matrix_in = (double*) calloc( Ncov, sizeof(double));
max_matrix_in = (double*) calloc( Ncov, sizeof(double));

matrix_window = (double **) calloc(N_in, sizeof(double*));
                 for(i=0;i<N_in;i++){matrix_window[i]= (double *) calloc(Ncov, sizeof(double));}

if(modeP0==1){P_theo0 = (double*) calloc( Nelements, sizeof(double));}
if(modeP2==1){P_theo2 = (double*) calloc( Nelements, sizeof(double));}
if(modeP4==1){P_theo4 = (double*) calloc( Nelements, sizeof(double));}

i_matrix=0;
for(mode1=0;mode1<=4;mode1=mode1+2){

    for(j_k=0;j_k<Nelements;j_k++)
    {

if(strcmp(spacing_data,"linear") == 0){k_theo[j_k]=(j_k+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);}

if(strcmp(spacing_data,"log") == 0){k_theo[j_k]=exp(  (j_k+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );}

if(strcmp(spacing_data,"log10") == 0){k_theo[j_k]=pow(10,  (j_k+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );}

if(strcmp(spacing_data,"irregular") == 0){k_theo[j_k]=Theory[j_k][0];}

vector_input=0;

if(modeP0==1 && mode1==0)
{
         ptheo=10000;
           P_theo0[j_k]=ptheo;pin=ptheo;vector_input++;
             fprintf(fpin,"%d 0 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
             
}

if(modeP2==1 && mode1==2)
{
         ptheo2=10000;
         P_theo2[j_k]=ptheo2;pin=ptheo2;vector_input++;
         fprintf(fpin,"%d 2 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
        
}

if(modeP4==1 && mode1==4)
{
         ptheo4=10000;
         P_theo4[j_k]=ptheo4;pin=ptheo4;vector_input++;
         fprintf(fpin,"%d 4 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
}

if(vector_input>0){//some-vector-filled
if(modeP0==1){P_theo0win = (double*) calloc( Nelements, sizeof(double));}
if(modeP2==1){P_theo2win = (double*) calloc( Nelements, sizeof(double));}
if(modeP4==1){P_theo4win = (double*) calloc( Nelements, sizeof(double));}

if(modeP0==1 && mode1==0){
P_theo0win[j_k]=P_theo0[j_k];P_theo0[j_k]=0;
}

if(modeP2==1 && mode1==2){
P_theo2win[j_k]=P_theo2[j_k];P_theo2[j_k]=0;
}
if(modeP4==1 && mode1==4){
P_theo4win[j_k]=P_theo4[j_k];P_theo4[j_k]=0;
}


apply_mask(type_of_analysis, modeP0, modeP2, modeP4, k_theo,P_theo0win, P_theo2win, P_theo4win,NeffP, pos, W0, W2, W4, W6,W8,Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,0,0);


if(modeP0==1){P0 = (double*) calloc( NeffP0, sizeof(double));}
if(modeP2==1){P2 = (double*) calloc( NeffP2, sizeof(double));}
if(modeP4==1){P4 = (double*) calloc( NeffP4, sizeof(double));}

j=0;
modeP0_trial=modeP0;
modeP2_trial=modeP2;
modeP4_trial=modeP4;

for(i_out=0;i_out<Ncov;i_out++)
{

       if(modeP4_trial==1 && modeP0_trial==0 && modeP2_trial==0){

            if(i_matrix==1){fprintf(fpout,"%d 4 %lf\n",i_out,k4[j]);}

       P4[j]= P_interpol(k4[j], k_theo, P_theo4win, Nelements);m=P4[j]/pin;
       
        j++;
        if(j==NeffP4){j=0;modeP4_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

        if(pin==0 || m==0){printf("Error generating the Window Matrix (P4, k=%lf): matrix[%d][%d]= %lf / %lf\n",k4[j],i_matrix,i_out,P4[j],pin);}
        }

        if(modeP2_trial==1 && modeP0_trial==0){

            if(i_matrix==1){fprintf(fpout,"%d 2 %lf\n",i_out,k2[j]);}

        P2[j]= P_interpol(k2[j], k_theo, P_theo2win, Nelements);m=P2[j]/pin;
      
        j++;
        if(j==NeffP2){j=0;modeP2_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

if(pin==0 || m==0){printf("Error generating the Window Matrix (P2, k=%lf): matrix[%d][%d]= %lf / %lf\n",k2[j],i_matrix,i_out,P2[j],pin);}

        }

        if(modeP0_trial==1){
        if(i_matrix==1){fprintf(fpout,"%d 0 %lf\n",i_out,k0[j]);}

        P0[j]= P_interpol(k0[j], k_theo, P_theo0win, Nelements);m=P0[j]/pin;
        
        j++;
        if(j==NeffP0){j=0;modeP0_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

        if(pin==0 || m==0){printf("Error generating the Window Matrix (P0, k=%lf): matrix[%d][%d]= %lf / %lf\n",k0[j],i_matrix,i_out,P0[j],pin);}

        }

}//i_out loop

if(modeP0==1){free(P0);free(P_theo0win);}
if(modeP2==1){free(P2);free(P_theo2win);}
if(modeP4==1){free(P4);free(P_theo4win);}

}//if-some-vector-filled

}//j_k-Nelements
}//mode1


for(i=0;i<i_matrix;i++)
{
    for(j=0;j<Ncov;j++)
    {
       fprintf(fmatrix,"%e\t",matrix_window[i][j]);
//    div=sum_matrix_in[j];//normalize by the sum of all elements of the line
       div=max_matrix_in[j];//normalize by the maximum value of the element of the line
        fprintf(fmatrix2,"%e\t",matrix_window[i][j]/div );
       if(j==Ncov-1){fprintf(fmatrix,"\n");}
       if(j==Ncov-1){fprintf(fmatrix2,"\n");}
    }

}

freeTokens(matrix_window,N_in);
free(sum_matrix_in);
free(max_matrix_in);
fclose(fpin);
fclose(fpout);
fclose(fmatrix);
fclose(fmatrix2);
if(modeP0==1){free(P_theo0);}
if(modeP2==1){free(P_theo2);}
if(modeP4==1){free(P_theo4);}
}
*/

/*
void do_Ptheo_RSD_window(char *ptmodel,char *rsdmodel_ps,char *fogmodel_ps, char *type_of_analysis, char *RSD_fit, char *fit_BAO,int modeP0,int modeP2,int modeP4,double **Theory,int N_Plin,  double k_theo[], double k_theo0[], double k_theo2[], double k_theo4[], double P_theo0[], double P_theo2[], double P_theo4[], double Pnoise,int NeffP0, int NeffP2, int NeffP4,int factor_for_sampling, double *parameters1,double *pos, double *W0,double *W2,double *W4,double *W6,double *W8,int Nmask,char *spacing_mask, double *k0, double *k2, double *k4, char *path_to_mask1,fftw_plan plan1,fftw_plan plan2,double k0min, double k0max, double k2min, double k2max, double k4min, double k4max, char *spacing_data,char *spacing_theory,int Neffmax,int Ncov, int interpolation_order,char *path_output,char *identifier, int ngcsgc)//delete
{
double **matrix_window;//Nelements*(modeP0+modeP2+modep4) x Ncov
double *sum_matrix_in,*max_matrix_in;
double *P0,*P2,*P4;
double pin,m;
int Ninterpol;
int i_out,i_matrix,mode1,vector_input,N_in;
int modeP0_trial,modeP2_trial,modeP4_trial;
int i,j;
int j_k;
double *P_theo0win, *P_theo2win, *P_theo4win;
double b1,b2,Anoise,sigma8_scaling,f,mBGV,m2BGV;
double bs2,b3nl,sigmaP,apara,aperp;
double ptheo,ptheo2,ptheo4,junk;
double XMIN1[1]={-1};
double XMAX1[1]={+1};
double precision=1e-1;
int Nelements;
double kmax,kmin;
int NeffP;
int NeffP_max;
double w1,w2,w0;
double div;
int shiftN;
char out_Pin[2000];sprintf(out_Pin,"%s/Pin_Window_%s_%d.txt",path_output,identifier,ngcsgc);
char out_Pout[2000];sprintf(out_Pout,"%s/Pout_Window_%s_%d.txt",path_output,identifier,ngcsgc);
char out_M[2000];sprintf(out_M,"%s/Matrix_in_Window_%s_%d.txt",path_output,identifier,ngcsgc);
char out_M2[2000];sprintf(out_M2,"%s/Matrix_offdiag_in_Window_%s_%d.txt",path_output,identifier,ngcsgc);
FILE *fpin,*fpout,*fmatrix,*fmatrix2;
fpin=fopen(out_Pin,"w");
fpout=fopen(out_Pout,"w");
fmatrix=fopen(out_M,"w");
fmatrix2=fopen(out_M2,"w");

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


b1=parameters1[6];
b2=parameters1[7];
Anoise=parameters1[8];
sigmaP=parameters1[11];
sigma8_scaling=pow(parameters1[5]/Theory[0][41],2);
f=parameters1[4];
mBGV=parameters1[2];
m2BGV=parameters1[3];
apara=parameters1[0];
aperp=parameters1[1];
bs2=parameters1[9];//-4./7.*(b1-1)
b3nl=parameters1[10];//32./315.*(b1-1)

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
(*function_parameters).mBGV=mBGV;
(*function_parameters).m2BGV=m2BGV;
(*function_parameters).sigma8=sigma8_scaling;
(*function_parameters).a_parallel=apara;
(*function_parameters).a_perpendicular=aperp;
(*function_parameters).theory=Theory;
(*function_parameters).spacing=spacing_theory;

Nelements=NeffP;
N_in=Nelements*(modeP0+modeP2+modeP4);
//matrix_window[N_in][N_out];
sum_matrix_in = (double*) calloc( Ncov, sizeof(double));
max_matrix_in = (double*) calloc( Ncov, sizeof(double));

matrix_window = (double **) calloc(N_in, sizeof(double*));
                 for(i=0;i<N_in;i++){matrix_window[i]= (double *) calloc(Ncov, sizeof(double));}


i_matrix=0;
for(mode1=0;mode1<=4;mode1=mode1+2){

    for(j_k=0;j_k<Nelements;j_k++)
    {

if(strcmp(spacing_data,"linear") == 0){k_theo[j_k]=(j_k+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j_k];}

if(strcmp(spacing_data,"log") == 0){k_theo[j_k]=exp(  (j_k+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j_k];}

if(strcmp(spacing_data,"log10") == 0){k_theo[j_k]=pow(10,  (j_k+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j_k];}

if(strcmp(spacing_data,"irregular") == 0){(*function_parameters).kinput=Theory[j_k][0];k_theo[j_k]=Theory[j_k][0];}

vector_input=0;

if(modeP0==1 && mode1==0)
{
         (*function_parameters).mode=0;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);
         //ptheo=10000;
          if(ptheo!=0 ){
           P_theo0[j_k]=ptheo-Pnoise;pin=ptheo-Pnoise;vector_input++;
             fprintf(fpin,"%d 0 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
              }
}

if(modeP2==1 && mode1==2)
{
         (*function_parameters).mode=2;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo2,&junk);
         //ptheo2=10000;
         if(ptheo2!=0){
         P_theo2[j_k]=ptheo2;pin=ptheo2;vector_input++;
         fprintf(fpin,"%d 2 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
        }
}

if(modeP4==1 && mode1==4)
{
         (*function_parameters).mode=4;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo4,&junk);
         //ptheo4=10000;
         if(ptheo4!=0){
         P_theo4[j_k]=ptheo4;pin=ptheo4;vector_input++;
         fprintf(fpin,"%d 4 %lf\n",i_matrix,k_theo[j_k]);i_matrix++;
         }
}

if(vector_input>0){//some-vector-filled
if(modeP0==1){P_theo0win = (double*) calloc( Nelements, sizeof(double));}
if(modeP2==1){P_theo2win = (double*) calloc( Nelements, sizeof(double));}
if(modeP4==1){P_theo4win = (double*) calloc( Nelements, sizeof(double));}

if(modeP0==1 && mode1==0){
P_theo0win[j_k]=P_theo0[j_k];P_theo0[j_k]=0;
}

if(modeP2==1 && mode1==2){
P_theo2win[j_k]=P_theo2[j_k];P_theo2[j_k]=0;
}
if(modeP4==1 && mode1==4){
P_theo4win[j_k]=P_theo4[j_k];P_theo4[j_k]=0;
}


apply_mask(type_of_analysis, modeP0, modeP2, modeP4, k_theo,P_theo0win, P_theo2win, P_theo4win,NeffP, pos, W0, W2, W4, W6,W8,Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,0,0);


if(modeP0==1){P0 = (double*) calloc( NeffP0, sizeof(double));}
if(modeP2==1){P2 = (double*) calloc( NeffP2, sizeof(double));}
if(modeP4==1){P4 = (double*) calloc( NeffP4, sizeof(double));}

j=0;
modeP0_trial=modeP0;
modeP2_trial=modeP2;
modeP4_trial=modeP4;

for(i_out=0;i_out<Ncov;i_out++)
{

       if(modeP4_trial==1 && modeP0_trial==0 && modeP2_trial==0){

            if(i_matrix==1){fprintf(fpout,"%d 4 %lf\n",i_out,k4[j]);}

       P4[j]= P_interpol(k4[j], k_theo, P_theo4win, Nelements);m=P4[j]/pin;
       
        j++;
        if(j==NeffP4){j=0;modeP4_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

        if(pin==0 || m==0){printf("Error generating the Window Matrix (P4, k=%lf): matrix[%d][%d]= %lf / %lf\n",k4[j],i_matrix,i_out,P4[j],pin);}
        }

        if(modeP2_trial==1 && modeP0_trial==0){

            if(i_matrix==1){fprintf(fpout,"%d 2 %lf\n",i_out,k2[j]);}

        P2[j]= P_interpol(k2[j], k_theo, P_theo2win, Nelements);m=P2[j]/pin;
      
        j++;
        if(j==NeffP2){j=0;modeP2_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

if(pin==0 || m==0){printf("Error generating the Window Matrix (P2, k=%lf): matrix[%d][%d]= %lf / %lf\n",k2[j],i_matrix,i_out,P2[j],pin);}

        }

        if(modeP0_trial==1){
        if(i_matrix==1){fprintf(fpout,"%d 0 %lf\n",i_out,k0[j]);}

        P0[j]= P_interpol(k0[j], k_theo, P_theo0win, Nelements);m=P0[j]/pin;
        
        j++;
        if(j==NeffP0){j=0;modeP0_trial=0;}
        matrix_window[i_matrix-1][i_out]=m;
        sum_matrix_in[i_out]=sum_matrix_in[i_out]+m;
        if(m>max_matrix_in[i_out]){max_matrix_in[i_out]=m;}

        if(pin==0 || m==0){printf("Error generating the Window Matrix (P0, k=%lf): matrix[%d][%d]= %lf / %lf\n",k0[j],i_matrix,i_out,P0[j],pin);}

        }

}//i_out loop

if(modeP0==1){free(P0);free(P_theo0win);}
if(modeP2==1){free(P2);free(P_theo2win);}
if(modeP4==1){free(P4);free(P_theo4win);}

}//if-some-vector-filled

}//j_k-Nelements
}//mode1
free(function_parameters);


for(i=0;i<i_matrix;i++)
{
    for(j=0;j<Ncov;j++)
    {
       fprintf(fmatrix,"%e\t",matrix_window[i][j]);
//    div=sum_matrix_in[j];//normalize by the sum of all elements of the line
       div=max_matrix_in[j];//normalize by the maximum value of the element of the line
        fprintf(fmatrix2,"%e\t",matrix_window[i][j]/div );
       if(j==Ncov-1){fprintf(fmatrix,"\n");}
       if(j==Ncov-1){fprintf(fmatrix2,"\n");}
    }

}

freeTokens(matrix_window,N_in);
free(sum_matrix_in);
free(max_matrix_in);
fclose(fpin);
fclose(fpout);
fclose(fmatrix);
fclose(fmatrix2);
}
*/

void do_Ptheo_RSD(char *ptmodel,char *rsdmodel_ps,char *fogmodel_ps, char *type_of_analysis, char *RSD_fit, char *fit_BAO,int modeP0,int modeP2,int modeP4,double **Theory,int N_Plin,  double k_theo[], double k_theo0[], double k_theo2[], double k_theo4[], double P_theo0[], double P_theo2[], double P_theo4[], double Pnoise,int NeffP0, int NeffP2, int NeffP4,int factor_for_sampling, double *parameters1,double *pos, double *W0,double *W2,double *W4,double *W6,double *W8,int Nmask,char *spacing_mask, double *k0, double *k2, double *k4, char *path_to_mask1,fftw_plan plan1,fftw_plan plan2,double k0min, double k0max, double k2min, double k2max, double k4min, double k4max, char *spacing_data,char *spacing_theory, char *mask_matrix, double **MatrixFS_mask, int noise_option)
{
int factor;
double params[2];
int params_int[2];
int i,j;
double b1,b2,Anoise,sigma8_scaling,f,mBGV,m2BGV,avir;
double bs2,b3nl,sigmaP,apara,aperp;
double ptheo,ptheo2,ptheo4,junk;
double XMIN1[1]={-1};
double XMAX1[1]={+1};
double precision=1e-1;
int Nelements;
double kmax,kmin;
int NeffP;
int NeffP_max;
double *vector_matrix_in;
int i_vector_matrix,Nin,Nout;

//printf("%lf %lf %lf %lf\n",k0[0],k0[1],k0[2],k0[3]);

if(strcmp(spacing_data,"irregular") == 0){factor=1;}
else{factor=factor_for_sampling;}

get_kmin_kmax(params,modeP0,modeP2,modeP4,k0min,k0max,k2min,k2max,k4min,k4max);
kmax=params[1];
kmin=params[0];
get_NeffP(params_int,NeffP0,NeffP2,NeffP4,spacing_data,kmin,kmax,k0min,k0max,k2min,k2max,k4min,k4max, N_Plin );
NeffP_max=params_int[1];
NeffP=params_int[0]*factor;

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}//mask
else{Nelements=NeffP0*modeP0;}//no-mask

if(  strcmp(mask_matrix,"yes") == 0){
Nin=Nelements*(modeP0+modeP2+modeP4);
Nout=NeffP0*modeP0+NeffP2*modeP2+NeffP4*modeP4;
vector_matrix_in = (double*) calloc( Nin, sizeof(double));
}

b1=parameters1[6];
b2=parameters1[7];
Anoise=parameters1[8];
sigmaP=parameters1[11];
avir=parameters1[13];
sigma8_scaling=pow(parameters1[5]/Theory[0][41],2);
f=parameters1[4];
mBGV=parameters1[2];
m2BGV=parameters1[3];
apara=parameters1[0];
aperp=parameters1[1];
bs2=parameters1[9];//-4./7.*(b1-1)
b3nl=parameters1[10];//32./315.*(b1-1)

        f_params *function_parameters;

        function_parameters = (f_params *) malloc(sizeof(f_params));

(*function_parameters).type_fog=fogmodel_ps;
(*function_parameters).type_ptmodel=ptmodel;
(*function_parameters).type_rsdmodel=rsdmodel_ps;

(*function_parameters).N=N_Plin;
(*function_parameters).sigmaP=sigmaP;
(*function_parameters).avir=avir;
(*function_parameters).b1=b1;
(*function_parameters).b2=b2;
(*function_parameters).bs2=bs2;
(*function_parameters).b3nl=b3nl;
(*function_parameters).A=Anoise;
(*function_parameters).Pnoise=Pnoise;
(*function_parameters).f=f;
(*function_parameters).mBGV=mBGV;
(*function_parameters).m2BGV=m2BGV;
(*function_parameters).sigma8=sigma8_scaling;
(*function_parameters).a_parallel=apara;
(*function_parameters).a_perpendicular=aperp;
(*function_parameters).theory=Theory;
(*function_parameters).spacing=spacing_theory;
(*function_parameters).noise_option=noise_option;

    for(j=0;j<Nelements;j++)
    {
        if(strcmp(path_to_mask1, "none") != 0){

k_theo[j]=get_ktheo(spacing_data,j,kmin,kmax,NeffP_max,factor_for_sampling,NULL,Theory);(*function_parameters).kinput=k_theo[j];
        }

    }

//for(i=0;i<NeffP0;i++)printf("data: %d, %lf\n",i,k0[i]);
//exit(0);
i_vector_matrix=0;
if(modeP0==1)
{
    for(j=0;j<Nelements;j++)
    {
        if(strcmp(path_to_mask1, "none") != 0){

//k_theo[j]=get_ktheo(spacing_data,j,kmin,kmax,NeffP_max,factor_for_sampling,NULL,Theory);
(*function_parameters).kinput=k_theo[j];
//if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"irregular") == 0){(*function_parameters).kinput=Theory[j][0];k_theo[j]=Theory[j][0];}


}
        else{k_theo0[j]=k0[j];(*function_parameters).kinput=k_theo0[j];}
         
         (*function_parameters).mode=0;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo,&junk);
         if(  strcmp(mask_matrix,"no") == 0){P_theo0[j]=ptheo;}
         if(  strcmp(mask_matrix,"yes") == 0){vector_matrix_in[i_vector_matrix]=ptheo;i_vector_matrix++;}
    }
}


if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP2*modeP2;}

if(modeP2==1)
{
    for(j=0;j<Nelements;j++)
    {
        if(strcmp(path_to_mask1, "none") != 0){

//k_theo[j]=get_ktheo(spacing_data,j,kmin,kmax,NeffP_max,factor_for_sampling,NULL,Theory);
     (*function_parameters).kinput=k_theo[j];

//if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"irregular") == 0){(*function_parameters).kinput=Theory[j][0];k_theo[j]=Theory[j][0];}


}
        else{k_theo2[j]=k2[j];(*function_parameters).kinput=k_theo2[j];}

         (*function_parameters).mode=2;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo2,&junk);
         if(  strcmp(mask_matrix,"no") == 0){P_theo2[j]=ptheo2;}
         if(  strcmp(mask_matrix,"yes") == 0){vector_matrix_in[i_vector_matrix]=ptheo2;i_vector_matrix++;}
    }
}

if(strcmp(path_to_mask1, "none") != 0){Nelements=NeffP;}
else{Nelements=NeffP4*modeP4;}

if(modeP4==1)
{
    for(j=0;j<Nelements;j++)
    {
        if(strcmp(path_to_mask1, "none") != 0){

//k_theo[j]=get_ktheo(spacing_data,j,kmin,kmax,NeffP_max,factor_for_sampling,NULL,Theory);
(*function_parameters).kinput=k_theo[j];

//if(strcmp(spacing_data,"linear") == 0){k_theo[j]=(j+0.5)*(kmax-kmin)/((NeffP_max-1.)*factor_for_sampling);  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"log") == 0){k_theo[j]=exp(  (j+0.5)*(log(kmax)-log(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"log10") == 0){k_theo[j]=pow(10,  (j+0.5)*(log10(kmax)-log10(kmin))/((NeffP_max-1.)*factor_for_sampling)  );  (*function_parameters).kinput=k_theo[j];}

//if(strcmp(spacing_data,"irregular") == 0){(*function_parameters).kinput=Theory[j][0];k_theo[j]=Theory[j][0];}

}

        else{k_theo4[j]=k4[j];(*function_parameters).kinput=k_theo4[j];}

         (*function_parameters).mode=4;
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&ptheo4,&junk);
         if(  strcmp(mask_matrix,"no") == 0){P_theo4[j]=ptheo4;}
         if(  strcmp(mask_matrix,"yes") == 0){vector_matrix_in[i_vector_matrix]=ptheo4;i_vector_matrix++;}
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
if(  strcmp(mask_matrix,"no") == 0){apply_mask(type_of_analysis, modeP0, modeP2, modeP4, k_theo,P_theo0, P_theo2, P_theo4,NeffP, pos, W0, W2, W4, W6,W8,Nmask,spacing_mask, plan1, plan2, k_theo[0], k_theo[NeffP-1], k0min, k0max, k2min, k2max, k4min, k4max,spacing_data,0,0);}

if(  strcmp(mask_matrix,"yes") == 0){

//apply matrix

apply_mask_matrix( P_theo0,P_theo2,P_theo4,MatrixFS_mask,vector_matrix_in,Nin,Nout, modeP0, modeP2, modeP4,NeffP0,NeffP2,NeffP4,type_of_analysis);
free(vector_matrix_in);
}


}

//for(i=0;i<Nelements;i++){printf("%lf %lf %lf\n",k_theo[i],P_theo0[i],P_theo2[i]);}


//exit(0);
}

double chi2_bao_rsd(char *type_BAO_fit, char *type_of_analysis,char *fit_BAO,char *fit_RSD, double *parameters2_bao,double  *parameters2_rsd, double *k_Plin, double *Plin,int N_Plin,double *k_Olin,double *Olin,int N_Olin,double **Theory,int Ntheory, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0bao, double *k0rsd, double *P0bao, double *P0rsd, double Pnoise,int NeffP0bao,int NeffP0rsd, double *k2bao, double *k2rsd, double *P2bao, double *P2rsd,int NeffP2bao,int NeffP2rsd, double *k4bao, double *k4rsd, double *P4bao, double *P4rsd, int NeffP4bao, int NeffP4rsd, double *k0baoSGC, double *k0rsdSGC, double *P0baoSGC, double *P0rsdSGC,double PnoiseSGC,int NeffP0baoSGC,int NeffP0rsdSGC, double *k2baoSGC,double *k2rsdSGC, double *P2baoSGC,double *P2rsdSGC,int NeffP2baoSGC,int NeffP2rsdSGC, double *k4baoSGC, double *k4rsdSGC, double *P4baoSGC,double *P4rsdSGC, int NeffP4baoSGC, int NeffP4rsdSGC,double *cov, double *covSGC, char *Sigma_def_type, char *Sigma_independent,  double ffactor,double *Sigma_type,  double *Sigma_nl_mean,  double *Sigma_nl_stddev, int Npolynomial, char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free, char *fog_bs, int Nchunks, fftw_plan plan1bao, fftw_plan plan2bao, fftw_plan plan1rsd, fftw_plan plan2rsd ,char *do_power_spectrum, char *do_bispectrum, int Nalphas, int Nsigmas_tot,int Nsigmas_free, double Sigma_smooth,int factor_sampling_mask_in,char *spacing_dataNGC_bao,char *spacing_dataNGC_rsd,char *spacing_dataSGC_bao,char *spacing_dataSGC_rsd, char *spacing_theory_bao, char *spacing_theory_rsd,char *bispectrum_BQ, char *mask_matrix, double **MatrixBAO_mask_NGC, double **MatrixBAO_mask_SGC, double **MatrixFS_mask_NGC, double **MatrixFS_mask_SGC, double *FSprior_type, double *FSprior_mean, double *FSprior_stddev,int noise_option,char *covariance_correction, int NrealNGC, int NrealSGC)
{
double prior;
double Anoise,apara,aperp;
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
prior=0;
if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}

factor_sampling_mask=factor_sampling_mask_in;
if( strcmp(spacing_dataNGC_bao,"irregular") == 0  ){factor_sampling_mask=1;}

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

if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0)//no mask
{

if(strcmp(path_to_mask1, "none") == 0){
if(NeffP0bao*modeP0bao>0){k_theo0bao = (double*) calloc( NeffP0bao, sizeof(double));}
if(NeffP2bao*modeP2bao>0){k_theo2bao = (double*) calloc( NeffP2bao, sizeof(double));}
if(NeffP4bao*modeP4bao>0){k_theo4bao = (double*) calloc( NeffP4bao, sizeof(double));}
}
if(strcmp(mask_matrix, "yes") == 0)
{
Neffmax=get_Neffmax(spacing_dataNGC_bao, modeP0bao, modeP2bao, modeP4bao, NeffP0bao, NeffP2bao, NeffP4bao, k0bao[0], k0bao[NeffP0bao-1], k2bao[0], k2bao[NeffP2bao-1], k4bao[0], k4bao[NeffP4bao-1], N_Plin );
k_theobao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}

if(NeffP0bao*modeP0bao>0){P_theo0bao = (double*) calloc( NeffP0bao, sizeof(double));}
if(NeffP2bao*modeP2bao>0){P_theo2bao = (double*) calloc( NeffP2bao, sizeof(double));}
if(NeffP4bao*modeP4bao>0){P_theo4bao = (double*) calloc( NeffP4bao, sizeof(double));}
}
else
{
Neffmax=get_Neffmax(spacing_dataNGC_bao, modeP0bao, modeP2bao, modeP4bao, NeffP0bao, NeffP2bao, NeffP4bao, k0bao[0], k0bao[NeffP0bao-1], k2bao[0], k2bao[NeffP2bao-1], k4bao[0], k4bao[NeffP4bao-1], N_Plin );
/*
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
*/

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


do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0bao,NeffP2bao,NeffP4bao,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin,pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0bao,k2bao,k4bao, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0bao[0],k0bao[NeffP0bao-1],k2bao[0],k2bao[NeffP2bao-1], k4bao[0],k4bao[NeffP4bao-1], 1,spacing_dataNGC_bao,spacing_theory_bao,Sigma_smooth, mask_matrix, MatrixBAO_mask_NGC);
}

if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){

//for(i=0;i<(modeP0bao+modeP2bao+modeP4bao)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1;i++){printf("BAO-NGC: %d %lf\n",i,parameters1[i]);}

do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0bao,NeffP2bao,NeffP4bao,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, Sigma_smooth, pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, path_to_mask1,k0bao,k2bao,k4bao, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0bao[0],k0bao[NeffP0bao-1],k2bao[0],k2bao[NeffP2bao-1], k4bao[0],k4bao[NeffP4bao-1], 1,spacing_dataNGC_bao,spacing_theory_bao, mask_matrix, MatrixBAO_mask_NGC);
//if(strcmp(path_to_mask1, "none") != 0 ){printf("%lf %lf %lf\n",k_theobao[300],P_theo0bao[300],P_theo2bao[300]);}
//if(strcmp(path_to_mask1, "none") == 0 ){printf("%lf %lf\n",P_theo0bao[10],P_theo2bao[10]);}

}

free(parameters1);

if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0)
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
if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0){offsetrsd++;}//avir free

//RSD shape2 this has to be 14
dimension=14;//This is a fixed quantity.


if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0)
{

if(strcmp(path_to_mask1, "none") == 0){
if(NeffP0rsd*modeP0rsd>0){k_theo0rsd = (double*) calloc( NeffP0rsd, sizeof(double));}
if(NeffP2rsd*modeP2rsd>0){k_theo2rsd = (double*) calloc( NeffP2rsd, sizeof(double));}
if(NeffP4rsd*modeP4rsd>0){k_theo4rsd = (double*) calloc( NeffP4rsd, sizeof(double));}
}
if(strcmp(mask_matrix, "yes") == 0)
{
Neffmax=get_Neffmax(spacing_dataNGC_rsd, modeP0rsd, modeP2rsd, modeP4rsd, NeffP0rsd, NeffP2rsd, NeffP4rsd, k0rsd[0], k0rsd[NeffP0rsd-1], k2rsd[0], k2rsd[NeffP2rsd-1], k4rsd[0], k4rsd[NeffP4rsd-1], Ntheory );
k_theorsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}

if(NeffP0rsd*modeP0rsd>0){P_theo0rsd = (double*) calloc( NeffP0rsd, sizeof(double));}
if(NeffP2rsd*modeP2rsd>0){P_theo2rsd = (double*) calloc( NeffP2rsd, sizeof(double));}
if(NeffP4rsd*modeP4rsd>0){P_theo4rsd = (double*) calloc( NeffP4rsd, sizeof(double));}
}
else{

Neffmax=get_Neffmax(spacing_dataNGC_rsd, modeP0rsd, modeP2rsd, modeP4rsd, NeffP0rsd, NeffP2rsd, NeffP4rsd, k0rsd[0], k0rsd[NeffP0rsd-1], k2rsd[0], k2rsd[NeffP2rsd-1], k4rsd[0], k4rsd[NeffP4rsd-1], Ntheory );

/*
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
*/

k_theorsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0rsd==1){P_theo0rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2rsd==1){P_theo2rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4rsd==1){P_theo4rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}


}

//difference = (double*) calloc( Ncov, sizeof(double));

parameters1 = (double*) calloc( dimension, sizeof(double));

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<5){
if(i==2){parameters1[i]=0;}//m1BGV
if(i==3){parameters1[i]=0;}//m2BGV
if(i!=2 && i!=3){parameters1[i]=parameters2_rsd[i1];
i1++;}
}

if(strcmp(RSD_fit, "no") == 0 && i<5){
if(i==0){parameters1[i]=1.;}//apara
if(i==1){parameters1[i]=1.;}//aperp
if(i==2){parameters1[i]=0;}//m1BGV
if(i==3){parameters1[i]=0;}//m2BGV
if(i==4){parameters1[i]=0;}//f
}

if(strcmp(RSD_fit, "shape") == 0 && i<5){
if(i==3){parameters1[i]=0;}//m2BGV
else{parameters1[i]=parameters2_rsd[i1];
i1++;}
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){//i==5
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==5 && strcmp(sigma8_free, "no") == 0 ){//i==5
parameters1[i]=Theory[0][41];
}
if(i>=6 && i<=8){//6, 8
parameters1[i]=parameters2_rsd[i1];
if(i==8){Anoise=parameters1[i]/(parameters1[0]*parameters1[1]*parameters1[1]);}

if(i==8 && FSprior_type[0] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[0],FSprior_stddev[0]);}//Anoise
if(i==7 && FSprior_type[1] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[1],FSprior_stddev[1]);}//b2

i1++;
}

if(i==9 && strcmp(local_b2s2, "no") == 0){//9
parameters1[i]=parameters2_rsd[i1];
if(FSprior_type[2] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[2],FSprior_stddev[2]);}//b2s2
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") != 0){//9
if(strcmp(local_b2s2, "yes") == 0){parameters1[i]=-4./7.*(parameters1[i-3]-1);}
if(strcmp(local_b2s2, "off") == 0){parameters1[i]=0;}

}

if(i==10 && strcmp(local_b3nl, "no") == 0){//10
parameters1[i]=parameters2_rsd[i1];
if(FSprior_type[3] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[3],FSprior_stddev[3]);}//b3nl
i1++;
}

if(i==10 && strcmp(local_b3nl, "no") != 0){//10
if(strcmp(local_b3nl, "yes") == 0){parameters1[i]=32./315.*(parameters1[i-4]-1);}
if(strcmp(local_b3nl, "off") == 0){parameters1[i]=0;}

}

if(i==11 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){//11
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(i==11 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){//11
parameters1[i]=0;
}

if(i==11 && strcmp(do_power_spectrum, "no") == 0){//11
parameters1[i]=0;
}

if(i==12 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){//12
parameters1[i]=parameters2_rsd[i1];
i1++;
}
if(i==12 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){//12
parameters1[i]=parameters1[i-1];
}

if(i==12 && strcmp(do_bispectrum, "no") == 0){//12
parameters1[i]=0;
}

if(i==12 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){//12
parameters1[i]=0;
}

if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") != 0){
parameters1[i]=0;
}
if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") == 0){

if( strcmp(fog_free, "yes") == 0 )
{
parameters1[i]=parameters2_rsd[i1];
i1++;
}
else{parameters1[i]=0;}

}


}

//for(i=0;i<dimension;i++){printf("RSD-NGC: %d %lf\n",i,parameters1[i]);}

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_RSD,modeP0rsd,modeP2rsd,modeP4rsd,Theory,Ntheory, k_theorsd,k_theo0rsd,k_theo2rsd,k_theo4rsd, P_theo0rsd, P_theo2rsd, P_theo4rsd,Pnoise, NeffP0rsd,NeffP2rsd,NeffP4rsd,factor_sampling_mask, parameters1,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC, k0rsd,k2rsd,k4rsd,path_to_mask1,plan1rsd, plan2rsd, k0rsd[0], k0rsd[NeffP0rsd-1], k2rsd[0], k2rsd[NeffP2rsd-1], k4rsd[0], k4rsd[NeffP4rsd-1],spacing_dataNGC_rsd,spacing_theory_rsd, mask_matrix, MatrixFS_mask_NGC,noise_option);
//if(strcmp(path_to_mask1, "none") != 0 ){printf("%lf %lf %lf %lf\n",k_theorsd[300],P_theo0rsd[300],P_theo2rsd[300],P_theo4rsd[300]);}
//if(strcmp(path_to_mask1, "none") == 0 ){printf("%lf %lf %lf\n",P_theo0rsd[10],P_theo2rsd[10],P_theo4rsd[10]);}

free(parameters1);

if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0)
{
i=-1;

   if(modeP0rsd==1){
        for(j=0;j<NeffP0rsd;j++){
        ptheo=P_theo0rsd[j]-Pnoise+Pnoise*Anoise*noise_option;
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

if(strcmp(covariance_correction,"Sellentin-Heavens") == 0)
{
ch2=NrealNGC*log(1+ch2/(NrealNGC*1.-1.));
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
if(strcmp(path_to_mask1, "none") == 0  ){free(k_theo0bao);}
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
if(strcmp(path_to_mask1, "none") == 0  ){free(k_theo4bao);}
}

if(strcmp(fit_RSD, "P4") == 0 || strcmp(fit_RSD, "P24") == 0 || strcmp(fit_RSD, "P04") == 0 || strcmp(fit_RSD, "P024") == 0)
{
free(P_theo4rsd);
if(strcmp(path_to_mask1, "none") == 0  ){free(k_theo4rsd);}
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

if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0 )
{
if(strcmp(path_to_mask2, "none") == 0){
if(NeffP0baoSGC*modeP0bao>0){k_theo0bao = (double*) calloc( NeffP0baoSGC, sizeof(double));}
if(NeffP2baoSGC*modeP2bao>0){k_theo2bao = (double*) calloc( NeffP2baoSGC, sizeof(double));}
if(NeffP4baoSGC*modeP4bao>0){k_theo4bao = (double*) calloc( NeffP4baoSGC, sizeof(double));}
}

if(strcmp(mask_matrix, "yes") == 0 )
{
Neffmax=get_Neffmax(spacing_dataSGC_bao, modeP0bao, modeP2bao, modeP4bao, NeffP0baoSGC, NeffP2baoSGC, NeffP4baoSGC, k0baoSGC[0], k0baoSGC[NeffP0baoSGC-1], k2baoSGC[0], k2baoSGC[NeffP2baoSGC-1], k4baoSGC[0], k4baoSGC[NeffP4baoSGC-1], N_Plin );
k_theobao = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}


if(NeffP0baoSGC*modeP0bao>0){P_theo0bao = (double*) calloc( NeffP0baoSGC, sizeof(double));}
if(NeffP2baoSGC*modeP2bao>0){P_theo2bao = (double*) calloc( NeffP2baoSGC, sizeof(double));}
if(NeffP4baoSGC*modeP4bao>0){P_theo4bao = (double*) calloc( NeffP4baoSGC, sizeof(double));}
}
else{
Neffmax=get_Neffmax(spacing_dataSGC_bao, modeP0bao, modeP2bao, modeP4bao, NeffP0baoSGC, NeffP2baoSGC, NeffP4baoSGC, k0baoSGC[0], k0baoSGC[NeffP0baoSGC-1], k2baoSGC[0], k2baoSGC[NeffP2baoSGC-1], k4baoSGC[0], k4baoSGC[NeffP4baoSGC-1], N_Plin );
/*
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
*/

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
do_Ptheo_multiple_iso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0baoSGC,NeffP2baoSGC,NeffP4baoSGC,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, spacing_maskSGC,path_to_mask2,k0baoSGC,k2baoSGC,k4baoSGC, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0baoSGC[0],k0baoSGC[NeffP0baoSGC-1],k2baoSGC[0],k2baoSGC[NeffP2baoSGC-1], k4baoSGC[0],k4baoSGC[NeffP4baoSGC-1], 1,spacing_dataSGC_bao,spacing_theory_bao,Sigma_smooth,mask_matrix, MatrixBAO_mask_SGC);
}

if( strcmp(type_of_analysis, "FSBAOANISO") == 0 ){

//for(i=0;i<(modeP0bao+modeP2bao+modeP4bao)*(Npolynomial)+1+Nalphas_for_param1+Nsigmas_for_param1+1;i++){printf("BAO-SGC: %d %lf\n",i,parameters1[i]);}


do_Ptheo_multiple_aniso(type_BAO_fit,type_of_analysis,fit_BAO,modeP0bao,modeP2bao,modeP4bao, k_theobao,k_theo0bao,k_theo2bao,k_theo4bao, P_theo0bao, P_theo2bao, P_theo4bao,NeffP0baoSGC,NeffP2baoSGC,NeffP4baoSGC,factor_sampling_mask, parameters1,k_Plin,Plin,N_Plin, k_Olin, Olin, N_Olin, Sigma_smooth,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, spacing_maskSGC, path_to_mask2,k0baoSGC,k2baoSGC,k4baoSGC, Npolynomial, plan1bao, plan2bao, k_Plin[0], k_Plin[N_Plin-1],k0baoSGC[0],k0baoSGC[NeffP0baoSGC-1],k2baoSGC[0],k2baoSGC[NeffP2baoSGC-1], k4baoSGC[0],k4baoSGC[NeffP4baoSGC-1], 1,spacing_dataSGC_bao,spacing_theory_bao, mask_matrix, MatrixBAO_mask_SGC);
}

free(parameters1);

if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0)
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

if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0)
{

if(strcmp(path_to_mask2, "none") == 0){
if(NeffP0rsdSGC*modeP0rsd>0){k_theo0rsd = (double*) calloc( NeffP0rsdSGC, sizeof(double));}
if(NeffP2rsdSGC*modeP2rsd>0){k_theo2rsd = (double*) calloc( NeffP2rsdSGC, sizeof(double));}
if(NeffP4rsdSGC*modeP4rsd>0){k_theo4rsd = (double*) calloc( NeffP4rsdSGC, sizeof(double));}
}
if(strcmp(mask_matrix, "yes") == 0)
{
Neffmax=get_Neffmax(spacing_dataSGC_rsd, modeP0rsd, modeP2rsd, modeP4rsd, NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, k0rsdSGC[0], k0rsdSGC[NeffP0rsdSGC-1], k2rsdSGC[0], k2rsdSGC[NeffP2rsdSGC-1], k4rsdSGC[0], k4rsdSGC[NeffP4rsdSGC-1], Ntheory );
k_theorsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}

if(NeffP0rsdSGC*modeP0rsd>0){P_theo0rsd = (double*) calloc( NeffP0rsdSGC, sizeof(double));}
if(NeffP2rsdSGC*modeP2rsd>0){P_theo2rsd = (double*) calloc( NeffP2rsdSGC, sizeof(double));}
if(NeffP4rsdSGC*modeP4rsd>0){P_theo4rsd = (double*) calloc( NeffP4rsdSGC, sizeof(double));}
}
else{
Neffmax=get_Neffmax(spacing_dataSGC_rsd, modeP0rsd, modeP2rsd, modeP4rsd, NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, k0rsdSGC[0], k0rsdSGC[NeffP0rsdSGC-1], k2rsdSGC[0], k2rsdSGC[NeffP2rsdSGC-1], k4rsdSGC[0], k4rsdSGC[NeffP4rsdSGC-1], Ntheory );
/*
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
*/

k_theorsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0rsd==1){P_theo0rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2rsd==1){P_theo2rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4rsd==1){P_theo4rsd = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
}


parameters1 =  (double*) calloc( dimension, sizeof(double));

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<5){
if(i==2){parameters1[i]=0;}//m1_BGV
if(i==3){parameters1[i]=0;}//m2_BGV
if(i!=2 && i!=3){
parameters1[i]=parameters2_rsd[i1];
i1++;}
}

if(strcmp(RSD_fit, "shape") == 0 && i<5){
if(i==3){parameters1[i]=0;}//m2_BGV
else{
parameters1[i]=parameters2_rsd[i1];
i1++;}
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
parameters1[i]=parameters2_rsd[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<5){
if(i==0){parameters1[i]=1.;}//apara
if(i==1){parameters1[i]=1.;}//aperp
if(i==2){parameters1[i]=0;}//m1
if(i==3){parameters1[i]=0;}//m2
if(i==4){parameters1[i]=0;}//f
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2_rsd[i1];
i1++;
}
if(i==5 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=6 && i<=8){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
if(i==8){Anoise=parameters1[i]/(parameters1[0]*parameters1[1]*parameters1[1]);}
if(i==8 && FSprior_type[0] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[0],FSprior_stddev[0]);}//Anoise
if(i==7 && FSprior_type[1] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[1],FSprior_stddev[1]);}//b2
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
if(FSprior_type[2] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[2],FSprior_stddev[2]);}//bs2
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") != 0){
if(strcmp(local_b2s2, "yes") == 0){parameters1[i]=-4./7.*(parameters1[i-3]-1);}
if(strcmp(local_b2s2, "off") == 0){parameters1[i]=0;}
}

if(i==10 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
if(FSprior_type[3] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[3],FSprior_stddev[3]);}//b3nl
i1++;
}
if(i==10 && strcmp(local_b3nl, "no") != 0){
if(strcmp(local_b3nl, "yes") == 0){parameters1[i]=32./315.*(parameters1[i-4]-1);}
if(strcmp(local_b3nl, "off") == 0){parameters1[i]=0;}
}

if(i==11 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}

if(i==11 && strcmp(fog_free, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=0;
}

if(i==11 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0  ){
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}

if(i==12 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0  ){
parameters1[i]=parameters1[i-1];
}

if(i==12 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(fog_bs, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "no") == 0  ){
parameters1[i]=0;
}

if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") != 0){
parameters1[i]=0;
}
if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") == 0){

if( strcmp(fog_free, "yes") == 0 )
{
parameters1[i]=parameters2_rsd[i1+offsetrsd];
i1++;
}
else{parameters1[i]=0;}


}

}

//for(i=0;i<dimension;i++){printf("RSD-SGC: %d %lf\n",i,parameters1[i]);}

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps,type_of_analysis, RSD_fit, fit_RSD,modeP0rsd,modeP2rsd,modeP4rsd,Theory,Ntheory, k_theorsd,k_theo0rsd,k_theo2rsd,k_theo4rsd, P_theo0rsd, P_theo2rsd, P_theo4rsd,PnoiseSGC,NeffP0rsdSGC,NeffP2rsdSGC,NeffP4rsdSGC,factor_sampling_mask, parameters1,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC,k0rsdSGC,k2rsdSGC,k4rsdSGC, path_to_mask2,plan1rsd, plan2rsd, k0rsdSGC[0], k0rsdSGC[NeffP0rsdSGC-1], k2rsdSGC[0], k2rsdSGC[NeffP2rsdSGC-1], k4rsdSGC[0], k4rsdSGC[NeffP4rsdSGC-1],spacing_dataSGC_rsd,spacing_theory_rsd, mask_matrix, MatrixFS_mask_SGC,noise_option);

free(parameters1);

if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0)
{
i=-1;

        if(modeP0rsd==1){
        for(j=0;j<NeffP0rsdSGC;j++){
        ptheo=P_theo0rsd[j]-PnoiseSGC+PnoiseSGC*Anoise*noise_option;
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

if(strcmp(covariance_correction,"Sellentin-Heavens") == 0)
{
ch2SGC=NrealSGC*log(1+ch2SGC/(NrealSGC*1.-1.));
}

if(strcmp(path_to_mask2, "none") != 0){free(k_theobao);}
if(strcmp(path_to_mask2, "none") != 0){free(k_theorsd);}

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

ch2=ch2+prior;//RSD priors

//printf("%lf %lf %lf (%d)\n",ch2-ch2SGC-prior_chi2,ch2SGC,prior_chi2,Ncov);
//printf("%lf %lf, %lf %lf, %lf %lf\n",ch2,ch2SGC,chi2_bao,chi2_baoSGC,chi2_rsd,chi2_rsdSGC);
//printf("%lf %lf %lf\n",ch2,ch2-ch2SGC,ch2SGC);
//exit(0);
return ch2;
}

double chi2_rsd_mcmc(char *type_of_analysis,double *parameters2, double **Theory,int N_Plin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC,double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0, double *P0, double Pnoise, double *k2, double *P2, double *k4, double *P4, double *k0SGC, double *P0SGC, double PnoiseSGC, double *k2SGC, double *P2SGC, double *k4SGC, double *P4SGC, int NeffP0, int NeffP2, int NeffP4,int NeffP0SGC, int NeffP2SGC, int NeffP4SGC, double *cov, double *covSGC,  char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free, char *fog_bs, int Nchunks, fftw_plan plan1, fftw_plan plan2, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum, double redshift_in, int factor_sampling_mask_in,char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory, double knl, double *n_func_final, double *sigma8_x, double *knl_y,int Nknl,double *k11,double *k22,double *k33, double *B0,double *Bnoise,int NeffB0,double *k11SGC,double *k22SGC,double *k33SGC,double *B0SGC,double *BnoiseSGC,int NeffB0SGC,char *bispectrum_BQ, char *mask_matrix, double **MatrixFS_mask_NGC, double **MatrixFS_mask_SGC,double *FSprior_type, double *FSprior_mean, double *FSprior_stddev,int noise_option,char *covariance_correction, int NrealNGC, int NrealSGC, double *alphasNGC, double *alphasSGC)
{
//double delete=0;
double Anoise,apara,aperp,prior;
double ch2,ch2SGC,prior_chi2;
double *difference,*parameters1,*k_theo,*P_theo0,*P_theo2,*P_theo4;
double *k_theo0,*k_theo2,*k_theo4;
int i,j,Ndifference;
int i1;
double ptheo,pobs;
int dimension;
int Ncov;
int modeP0,modeP2,modeP4,modeB0,modealphas;
int Ntheo=N_Plin;
int offset;
int Neffmax;
int factor_sampling_mask;
int interpolation_order,shiftN;
double w1,w2,w0;
int Ninterpol;
double *b_theo;
double knl_eff;
interpolation_order=1;
prior=0;
if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}
//for(i=0;i<4;i++){printf("FSprior_type[%d]=%lf\n",i,FSprior_type[i]);}


factor_sampling_mask=factor_sampling_mask_in;
if( strcmp(spacing_dataNGC,"irregular") == 0  ){factor_sampling_mask=1;}

offset=3;
if(strcmp(local_b2s2, "no") == 0){offset++;}//b2s2
if(strcmp(local_b3nl, "no") == 0){offset++;}//b3nl
if(strcmp(fog_free, "yes") == 0){offset++;}//fog_ps
if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){offset++;}//fog_bs
if(strcmp(fogmodel_ps, "Exponential_avir") == 0 && strcmp(fog_free, "yes") == 0){offset++;}//avir free

dimension=14;//This is a fixed quantity. b1,b2,b2s,b3nl,apara,aperp,f,s8,m1BGV,m2BGV,sigmaP,Anoise,sigmaB,avir. 13 parameters maximum

Ncov=0;
modeP0=0;
modeP2=0;
modeP4=0;
modeB0=0;
modealphas=0;
if( strcmp(do_power_spectrum,"yes") == 0 ){
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}
if(alphasNGC!=NULL){modealphas=1;}

Ncov=Ncov+NeffP0*modeP0+NeffP2*modeP2+NeffP4*modeP4+2*modealphas;
}
if( strcmp(do_bispectrum,"yes") == 0 ){
Ncov=Ncov+NeffB0;
b_theo = (double*) calloc( NeffB0, sizeof(double));
modeB0=1;
}

if( strcmp(do_power_spectrum,"yes") == 0 ){

if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix,"yes") == 0)
{

if(strcmp(path_to_mask1, "none") == 0){
if(NeffP0*modeP0>0){k_theo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2*modeP2>0){k_theo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4*modeP4>0){k_theo4 = (double*) calloc( NeffP4, sizeof(double));}
}
if(strcmp(mask_matrix,"yes") == 0){
Neffmax=get_Neffmax(spacing_dataNGC, modeP0, modeP2, modeP4, NeffP0, NeffP2, NeffP4, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1], Ntheo );
k_theo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}

if(NeffP0*modeP0>0){P_theo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2*modeP2>0){P_theo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4*modeP4>0){P_theo4 = (double*) calloc( NeffP4, sizeof(double));}
}
else{
//k_theo = (double*) calloc( Ntheo, sizeof(double));
//if(modeP0==1){P_theo0 = (double*) calloc( Ntheo, sizeof(double));}
//if(modeP2==1){P_theo2 = (double*) calloc( Ntheo, sizeof(double));}
//if(modeP4==1){P_theo4 = (double*) calloc( Ntheo, sizeof(double));}

Neffmax=get_Neffmax(spacing_dataNGC, modeP0, modeP2, modeP4, NeffP0, NeffP2, NeffP4, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1], Ntheo );

/*
if( strcmp(spacing_dataNGC,"linear") == 0  ){

//printf("%lf %lf %d\n",k0[0],k0[NeffP0-1],NeffP0);
//printf("%lf %lf %d\n",k2[0],k2[NeffP2-1],NeffP2);
//printf("%lf %lf %d\n",k4[0],k4[NeffP4-1],NeffP4);

if(NeffP0*modeP0>=NeffP2*modeP2){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}
if(NeffP0*modeP0>=NeffP4*modeP4){Neffmax=NeffP0-2+25+(int)(k0[0]/(k0[NeffP0-1]-k0[0])*(NeffP0-1));}

if(NeffP2*modeP2>=NeffP0*modeP0){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}
if(NeffP2*modeP2>=NeffP4*modeP4){Neffmax=NeffP2-2+25+(int)(k2[0]/(k2[NeffP2-1]-k2[0])*(NeffP2-1));}

if(NeffP4*modeP4>=NeffP0*modeP0){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}
if(NeffP4*modeP4>=NeffP2*modeP2){Neffmax=NeffP4-2+25+(int)(k4[0]/(k4[NeffP4-1]-k4[0])*(NeffP4-1));}

//printf("Neffmax=%d,factor=%d\n",Neffmax,factor_sampling_mask);


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
*/

k_theo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0==1){P_theo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2==1){P_theo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4==1){P_theo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}


}//else

}//do_power_spectrum

difference = (double*) calloc( Ncov, sizeof(double));

parameters1 = (double*) calloc( dimension, sizeof(double));

//parameters1 are ordered as follows: alphapara,alphaperp,f,s8,b1,b2,A,b2s,b3nl,sigmaP,sigmaB

i1=0;
for(i=0;i<dimension;i++){

//apara,aperp,f
if(strcmp(RSD_fit, "yes") == 0 && i<5){
if(i==2){parameters1[i]=0;}//m1BGV
if(i==3){parameters1[i]=0;}//m2BGV
if(i!=2 && i!=3){parameters1[i]=parameters2[i1];
i1++;}
}

if(strcmp(RSD_fit, "shape") == 0 && i<5){
if(i==3){parameters1[i]=0;}//mBGV
else{parameters1[i]=parameters2[i1];
i1++;}
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
parameters1[i]=parameters2[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<5){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}//m1BGV
if(i==3){parameters1[i]=0;}//m2BGV
if(i==4){parameters1[i]=0;}//f
}

//sigma8
if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==5 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=6 && i<=8){//b1,b2,A (these are always present in parameters2)
parameters1[i]=parameters2[i1];
if(i==8 && FSprior_type[0] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[0],FSprior_stddev[0]);/*printf("A: %lf, %lf, %lf (%lf)\n",parameters1[i],FSprior_mean[0],FSprior_stddev[0],prior);*/}//Anoise
if(i==7 && FSprior_type[1] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[1],FSprior_stddev[1]);/*printf("b2: %lf, %lf, %lf\n", parameters1[i],FSprior_mean[1],FSprior_stddev[1]);*/}//b2
if(i==8){Anoise=parameters1[i]/(parameters1[0]*parameters1[1]*parameters1[1]);}

i1++;
}

//b2s2
if(i==9 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2[i1];
if(FSprior_type[2] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[2],FSprior_stddev[2]);/*printf("b2s2: %lf, %lf, %lf\n", parameters1[i],FSprior_mean[2],FSprior_stddev[2]);*/}//b2s2
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") != 0){
if(strcmp(local_b2s2, "yes") == 0){parameters1[i]=-4./7.*(parameters1[i-3]-1);}
if(strcmp(local_b2s2, "off") == 0){parameters1[i]=0;}
}

//b3nl
if(i==10 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2[i1];
if(FSprior_type[3] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[3],FSprior_stddev[3]);/*printf("b3nl: %lf, %lf, %lf\n", parameters1[i],FSprior_mean[3],FSprior_stddev[3]);*/}//b3nl
i1++;
}

if(i==10 && strcmp(local_b3nl, "no") != 0){
if(strcmp(local_b3nl, "yes") == 0){parameters1[i]=32./315.*(parameters1[i-4]-1);}
if(strcmp(local_b3nl, "off") == 0){parameters1[i]=0;}
}

//sigmaP
if(i==11 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==11 && strcmp(fog_free, "no") == 0){
parameters1[i]=0;
}

if(i==11 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

//sigmaB
if(i==12 && strcmp(fog_free, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "yes") == 0){//if fog_bs = yes implies power_spectrum = yes
parameters1[i]=parameters1[i-1];
}

if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") != 0){
parameters1[i]=0;
}
if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") == 0){

if( strcmp(fog_free, "yes") == 0 )
{
parameters1[i]=parameters2[i1];
i1++;
}
else{parameters1[i]=0;}

}


}//for i


//printf("Noise %lf\n",Pnoise);
//for(i=0;i<dimension;i++){printf("NGC %d %lf\n",i,parameters1[i]);}
//exit(0);
//call here the Pmodel 

if( strcmp(do_power_spectrum,"yes") == 0){

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,Pnoise, NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC, k0,k2,k4,path_to_mask1,plan1, plan2, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1],spacing_dataNGC,spacing_theory, mask_matrix, MatrixFS_mask_NGC,noise_option);
//printf("%lf %lf %lf\n",P_theo0[0],P_theo2[0],P_theo4[0]);
//printf("%lf %lf %lf\n",P_theo0[NeffP0-1],P_theo2[NeffP2-1],P_theo4[NeffP4-1]);
}
//printf("\n ok \n");
//exit(0);
if( strcmp(do_bispectrum,"yes") == 0){

//get knl for that s8
if( strcmp(ptmodel_bs,"GilMarin14") ==  0){

if( strcmp(sigma8_free, "no") == 0 ){knl_eff=knl;}
else{knl_eff=P_interpol_extrapol(parameters1[3], sigma8_x, knl_y, Nknl);}//need to reverse the order of sigma8?

}else{knl_eff=0;/*printf("Warning knl_eff=0\n");exit(0);*/}


//printf("knl=%lf\n",knl_eff);
do_Btheo_RSD(ptmodel_ps,fogmodel_ps,ptmodel_bs,b_theo,NeffB0,Bnoise,parameters1,Theory,N_Plin,pos, W0, W2, W4, W6, W8,Nmask,spacing_maskNGC,path_to_mask1,k11,k22,k33,plan1, plan2,k11[0],k11[NeffB0-1],spacing_theory,knl_eff,n_func_final,bispectrum_BQ,mask_matrix, MatrixFS_mask_NGC,noise_option,redshift_in);

}

//if(strcmp(path_to_mask1, "none") == 0 ){for(i=0;i<NeffP0;i++){printf("%lf %lf\n",k_theo0[i],P_theo0[i]-Pnoise);}}
//if(strcmp(path_to_mask1, "none") != 0 ){for(i=0;i< Neffmax*factor_sampling_mask;i++){printf("%lf %lf\n",k_theo[i],P_theo0[i]-Pnoise);}}
//for(i=0;i<NeffB0;i++){printf("%d %lf %lf %lf %e %e %e\n",i,k11[i],k22[i],k33[i],b_theo[i]-Bnoise[i],Bnoise[i],B0[i]);}
//exit(0);

//for(i=0;i<dimension;i++){printf("NGC: %d %lf\n",i,parameters1[i]);}

//printf("%lf %lf %lf %lf\n",k_theo[300],P_theo0[300],P_theo2[300],P_theo4[300]);
//printf("%lf %lf %lf\n",P_theo0[10],P_theo2[10],P_theo4[10]);
free(parameters1);

if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix,"yes") == 0)
{
i=-1;

   if(modealphas==1)
   {
     for(j=0;j<2;j++){
     ptheo=parameters2[j];
     pobs=alphasNGC[j];
     i++;
     difference[i]=ptheo-pobs;//printf("alphas %d %lf %lf %lf\n",i,difference[i],ptheo,pobs);
     }
   }

   if(modeP0==1){
        for(j=0;j<NeffP0;j++){
        ptheo=P_theo0[j]-Pnoise+Pnoise*Anoise*noise_option;
        pobs=P0[j];
        i++;
        difference[i]=ptheo-pobs;//printf("0 %d %lf %lf %lf\n",i,difference[i],ptheo,pobs);//printf("%lf %lf %lf\n",k_theo0[j],P0[j],P_theo0[j]);
//printf("(%d,%d,%d,%d), %d %lf,%lf, %lf (%d/%d)\n",modeP0,modeP2,modeP4,modeB0,i,ptheo,pobs,ptheo-pobs,j,NeffP0);
        }
        }

        if(modeP2==1){
        for(j=0;j<NeffP2;j++){
        ptheo=P_theo2[j];
        pobs=P2[j];
        i++;
        difference[i]=ptheo-pobs;//printf("2 %d %lf %lf %lf\n",i,difference[i],ptheo,pobs);
//printf("(%d,%d,%d,%d), %d %lf,%lf, %lf (%d/%d)\n",modeP0,modeP2,modeP4,modeB0,i,ptheo,pobs,ptheo-pobs,j,NeffP2);
        }
        }

        if(modeP4==1){
        for(j=0;j<NeffP4;j++){
        ptheo=P_theo4[j];
        pobs=P4[j];
        i++;
        difference[i]=ptheo-pobs;//printf("4 %d %lf %lf %lf\n",i,difference[i],ptheo,pobs);
//printf("(%d,%d,%d,%d), %d %lf,%lf, %lf (%d/%d)\n",modeP0,modeP2,modeP4,modeB0,i,ptheo,pobs,ptheo-pobs,j,NeffP4);
        }
        }
   
        //add bispectrum here
       if( modeB0 == 1){

           for(j=0;j<NeffB0;j++){
           ptheo=b_theo[j]-Bnoise[j];
           pobs=B0[j];
           i++;
           difference[i]=ptheo-pobs;
//printf("(%d,%d,%d,%d), %d %lf,%lf, %lf (%d/%d)\n",modeP0,modeP2,modeP4,modeB0,i,ptheo,pobs,ptheo-pobs,j,NeffB0);
           }
       }


}
else{//mask

j=0;
    for(i=0;i<Ncov;i++)
    {

       

        if(modeB0==1 && modeP0==0 && modeP2==0 && modeP4== 0 && modealphas == 0){

           ptheo=b_theo[j]-Bnoise[j];
           pobs=B0[j];
           j++;
           if(j==NeffB0){j=0;modeB0=0;}

           }

        if(modeP4==1 && modeP0==0 && modeP2==0 && modealphas==0){

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

        if(modeP2==1 && modeP0==0 && modealphas==0){
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

        if(modeP0==1 && modealphas==0){
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
ptheo=P_interpol_fast(k0[j],P_theo0,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2)-Pnoise+Pnoise*Anoise*noise_option;
}


        pobs=P0[j];
        j++;
        if(j==NeffP0){j=0;modeP0=0;}
        }

        if(modealphas==1){
        ptheo=parameters2[j];
        pobs=alphasNGC[j];
        j++;
        if(j==2){j=0;modealphas=0;}
        }

                     
//printf("(%d,%d,%d,%d), %d %lf,%lf, %lf (%d) (%d,%d,%d,%d)\n",modeP0,modeP2,modeP4,modeB0,i,ptheo,pobs,ptheo-pobs,j-1,NeffP0,NeffP2,NeffP4,NeffB0);
        difference[i]=ptheo-pobs;//printf("%d, (%d,%d,%d,%d) %lf = %lf - %lf\n",i,modeP0,modeP2,modeP4,modealphas,difference[i],ptheo,pobs);
   }//for i

}//else  mask

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
//if(i<2 && j<2){delete=delete+difference[i]*1./cov[i+Ncov*j]*difference[j];printf("prior (%d,%d) %lf %lf\n",i,j,delete,difference[i]*1./cov[i+Ncov*j]*difference[j]);}
//if(i==j){printf("%d %lf %lf\n",i,difference[i],sqrt(cov[i+Ncov*j]));}

      }

}

if(strcmp(covariance_correction,"Sellentin-Heavens") == 0)
{
ch2=NrealNGC*log(1+ch2/(NrealNGC*1.-1.));
}

//printf("chi2=%lf\n",ch2);
//exit(0);

free(difference);

if( strcmp(do_power_spectrum,"yes") == 0){

if(strcmp(path_to_mask1, "none") != 0 ){free(k_theo);}

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo0);
if(strcmp(path_to_mask1, "none") == 0){free(k_theo0);}
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

}

if( strcmp(do_bispectrum,"yes") == 0){
free(b_theo);
}

if(Nchunks==2)//repetate for Nchunk=2
{
factor_sampling_mask=factor_sampling_mask_in;
modeP0=0;
modeP2=0;
modeP4=0;
modeB0=0;
modealphas=0;
Ncov=0;
if (strcmp(do_power_spectrum,"yes") ==  0){

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}
if(alphasSGC!=NULL ){modealphas=1;}

Ncov=Ncov+NeffP0SGC*modeP0+NeffP2SGC*modeP2+NeffP4SGC*modeP4+2*modealphas;

if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix,"yes") == 0 )
{
if(strcmp(path_to_mask2, "none") == 0)
{
if(NeffP0SGC*modeP0>0){k_theo0 = (double*) calloc( NeffP0SGC, sizeof(double));}
if(NeffP2SGC*modeP2>0){k_theo2 = (double*) calloc( NeffP2SGC, sizeof(double));}
if(NeffP4SGC*modeP4>0){k_theo4 = (double*) calloc( NeffP4SGC, sizeof(double));}
}
if(strcmp(mask_matrix,"yes") == 0 )
{
Neffmax=get_Neffmax(spacing_dataSGC, modeP0, modeP2, modeP4, NeffP0SGC, NeffP2SGC, NeffP4SGC, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1], Ntheo );
k_theo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}

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
Neffmax=get_Neffmax(spacing_dataSGC, modeP0, modeP2, modeP4, NeffP0SGC, NeffP2SGC, NeffP4SGC, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1], Ntheo );
/*
if(NeffP0SGC*modeP0>=NeffP2SGC*modeP2){Neffmax=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}
if(NeffP0SGC*modeP0>=NeffP4SGC*modeP4){Neffmax=NeffP0SGC-2+25+(int)(k0SGC[0]/(k0SGC[NeffP0SGC-1]-k0SGC[0])*(NeffP0SGC-1));}

if(NeffP2SGC*modeP2>=NeffP0SGC*modeP0){Neffmax=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}
if(NeffP2SGC*modeP2>=NeffP4SGC*modeP4){Neffmax=NeffP2SGC-2+25+(int)(k2SGC[0]/(k2SGC[NeffP2SGC-1]-k2SGC[0])*(NeffP2SGC-1));}

if(NeffP4SGC*modeP4>=NeffP0SGC*modeP0){Neffmax=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
if(NeffP4SGC*modeP4>=NeffP2SGC*modeP2){Neffmax=NeffP4SGC-2+25+(int)(k4SGC[0]/(k4SGC[NeffP4SGC-1]-k4SGC[0])*(NeffP4SGC-1));}
*/

/*
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
*/

k_theo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0==1){P_theo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2==1){P_theo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4==1){P_theo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}

}//else mask

}//do power


if(strcmp(do_bispectrum,"yes") == 0){
modeB0=1;
b_theo = (double*) calloc( NeffB0SGC, sizeof(double));
Ncov=Ncov+NeffB0SGC;
}


difference = (double*) calloc( Ncov, sizeof(double));
parameters1 =  (double*) calloc( dimension, sizeof(double));

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<5){
if(i==2){parameters1[i]=0;}
if(i==3){parameters1[i]=0;}
if(i!=2 && i!=3){
parameters1[i]=parameters2[i1];
i1++;}
}

if(strcmp(RSD_fit, "shape") == 0 && i<5){
if(i==3){parameters1[i]=0;}
else{
parameters1[i]=parameters2[i1];
i1++;}
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
parameters1[i]=parameters2[i1];
i1++;
}

if(strcmp(RSD_fit, "no") == 0 && i<5){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}
if(i==3){parameters1[i]=0;}
if(i==4){parameters1[i]=0;}
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==5 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=6 && i<=8){
parameters1[i]=parameters2[i1+offset];
if(i==8){Anoise=parameters1[i]/(parameters1[0]*parameters1[1]*parameters1[1]);}
if(i==8 && FSprior_type[0] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[0],FSprior_stddev[0]);}//Anoise
if(i==7 && FSprior_type[1] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[1],FSprior_stddev[1]);}//b2
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2[i1+offset];
if(FSprior_type[2] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[2],FSprior_stddev[2]);}//b2s2
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") != 0){
if(strcmp(local_b2s2, "yes") == 0){parameters1[i]=-4./7.*(parameters1[i-3]-1);}
if(strcmp(local_b2s2, "off") == 0){parameters1[i]=0;}
}

if(i==10 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2[i1+offset];
if(FSprior_type[3] == 1) {prior=prior+gauss(parameters1[i],FSprior_mean[3],FSprior_stddev[3]);}//3nl
i1++;
}
if(i==10 && strcmp(local_b3nl, "no") != 0){
if(strcmp(local_b3nl, "yes") == 0){parameters1[i]=32./315.*(parameters1[i-4]-1);}
if(strcmp(local_b3nl, "off") == 0){parameters1[i]=0;}
}

//sigmaP
if(i==11 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2[i1+offset];
i1++;
}

if(i==11 && strcmp(fog_free, "no") == 0){
parameters1[i]=0;
}

if(i==11 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

//sigmaB
if(i==12 && strcmp(fog_free, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "no") == 0){
parameters1[i]=parameters2[i1+offset];
i1++;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "yes") == 0){//if fog_bs = yes implies power_spectrum = yes
parameters1[i]=parameters1[i-1];
}

if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") != 0){
parameters1[i]=0;
}
if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") == 0){

if( strcmp(fog_free, "yes") == 0 )
{
parameters1[i]=parameters2[i1+offset];
i1++;
}
else{parameters1[i]=0;}

}


}

//for(i=0;i<dimension;i++){printf("SGC %d %lf\n",i,parameters1[i]);}

if( strcmp(do_power_spectrum,"yes") == 0){

do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps,type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, k_theo,k_theo0,k_theo2,k_theo4, P_theo0, P_theo2, P_theo4,PnoiseSGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC,k0SGC,k2SGC,k4SGC, path_to_mask2,plan1, plan2, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1],spacing_dataSGC,spacing_theory, mask_matrix, MatrixFS_mask_SGC,noise_option);

}
if( strcmp(do_bispectrum,"yes") == 0){

do_Btheo_RSD(ptmodel_ps,fogmodel_ps,ptmodel_bs,b_theo,NeffB0SGC,BnoiseSGC,parameters1,Theory,N_Plin,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC,NmaskSGC,spacing_maskSGC,path_to_mask2,k11SGC,k22SGC,k33SGC,plan1, plan2,k11SGC[0],k11SGC[NeffB0SGC-1],spacing_theory,knl_eff,n_func_final,bispectrum_BQ,mask_matrix, MatrixFS_mask_SGC,noise_option,redshift_in);

}

//for(i=0;i<dimension;i++){printf("SGC: %d %lf\n",i,parameters1[i]);}
//exit(0);

free(parameters1);

if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix,"yes") == 0)
{
i=-1;

   if(modealphas==1)
   {
     for(j=0;j<2;j++){
     ptheo=parameters2[j];
     pobs=alphasSGC[j];
     i++;
     difference[i]=ptheo-pobs;
     }
   }

        if(modeP0==1){
        for(j=0;j<NeffP0SGC;j++){
        ptheo=P_theo0[j]-PnoiseSGC+Anoise*PnoiseSGC*noise_option;
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
         if( modeB0 == 1){

           for(j=0;j<NeffB0SGC;j++){
           ptheo=b_theo[j]-BnoiseSGC[j];
           pobs=B0SGC[j];
           i++;
           difference[i]=ptheo-pobs;
           }
       }
}
else{


j=0;
    for(i=0;i<Ncov;i++)
    {

           if( modeB0 ==  1 && modeP0 == 0 && modeP2 == 0 && modeP4 == 0 && modealphas == 0){

           ptheo=b_theo[j]-BnoiseSGC[j];
           pobs=B0SGC[j];
           j++;
           if(j==NeffB0){modeB0=0;j=0;}

           }

        if(modeP4==1 && modeP0==0 && modeP2==0 && modealphas == 0){
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

        if(modeP2==1 && modeP0==0 && modealphas == 0){
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

        if(modeP0==1 && modealphas == 0){
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
ptheo=P_interpol_fast(k0SGC[j],P_theo0,Neffmax*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2)-PnoiseSGC+Anoise*PnoiseSGC*noise_option;
}


        pobs=P0SGC[j];//printf("%lf %lf %lf\n",k0SGC[j],P0SGC[j],ptheo);
        j++;
        if(j==NeffP0SGC){j=0;modeP0=0;}
        }

        if(modealphas == 1){
        ptheo=parameters2[j];
        pobs=alphasSGC[j];
        j++;
        if(j==2){j=0;modealphas=0;}
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

//if(i<2 && j<2){delete=delete+difference[i]*1./covSGC[i+Ncov*j]*difference[j];printf("prior (%d,%d) %lf %lf\n",i,j,delete,difference[i]*1./covSGC[i+Ncov*j]*difference[j]);}

//if(i==j){printf("%d %lf %lf\n",i,difference[i],sqrt(covSGC[i+Ncov*j]));}
      }

}

if(strcmp(covariance_correction,"Sellentin-Heavens") == 0)
{
ch2SGC=NrealSGC*log(1+ch2SGC/(NrealSGC*1.-1.));
}

free(difference);

if( strcmp(do_power_spectrum,"yes") == 0){

if(strcmp(path_to_mask2, "none") != 0){free(k_theo);}

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo0);
if(strcmp(path_to_mask2, "none") == 0){free(k_theo0);}
}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo2);
if(strcmp(path_to_mask2, "none") == 0){free(k_theo2);}
}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0)
{
free(P_theo4);
if(strcmp(path_to_mask2, "none") == 0){free(k_theo4);}
}

}
if( strcmp(do_bispectrum,"yes") == 0){

free(b_theo);

}

ch2=ch2+ch2SGC;
}//Nchunks==2
//printf("chi=%lf = %lf + %lf\n",ch2,ch2-ch2SGC,ch2SGC);
//exit(0);

ch2=ch2+prior;
//printf("ch2+p=%lf (%lf), AN=%lf, AS=%lf, Aprior=%lf pm %lf\n",ch2,prior, parameters2[6], parameters2[10],FSprior_mean[0],FSprior_stddev[0]);
//printf("ch2+p=%lf (%lf)\n",ch2,prior);
//printf("priors by alphas %lf\n",delete);
//exit(0);
return ch2;
}

void make_a_rsd_plot(char *type_of_analysis,double *parameters2,double chi2_min, double **Theory,int N_Plin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC,double *W6SGC,double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0, double *P0, double *errP0, double Pnoise, double *k2, double *P2, double *errP2, double *k4, double *P4, double *errP4, double *k0SGC, double *P0SGC, double *errP0SGC, double PnoiseSGC, double *k2SGC, double *P2SGC, double *errP2SGC, double *k4SGC, double *P4SGC, double *errP4SGC, int NeffP0, int NeffP2, int NeffP4,int NeffP0SGC, int NeffP2SGC, int NeffP4SGC, double *cov, double *covSGC,  char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free, char *fog_bs, int Nchunks,int Npoints, int Ndof, char *path_output, char *identifier, fftw_plan plan1, fftw_plan plan2, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum, double redshift_in, int factor_sampling_mask_in, char *spacing_dataNGC,char *spacing_dataSGC,char *spacing_theory,double knl, double *n_func_final, double *sigma8_x, double *knl_y, int Nknl, double *k11, double *k22, double *k33, double *B0, double *errB0, double *Bnoise,int NeffB0,double *k11SGC, double *k22SGC,double *k33SGC,double *B0SGC, double *errB0SGC, double *BnoiseSGC,int NeffB0SGC, char *bispectrum_BQ, char *mask_matrix, double **MatrixFS_mask_NGC, double **MatrixFS_mask_SGC,int noise_option)
{
double Anoise_NGC,Anoise_SGC;
FILE *f1,*f2,*f12;
int open,NeffP;
char name1[2000],name2[2000],name12[2000];
char nameP0_1[2000],nameP0_2[2000],nameP0_12[2000];
char nameP2_1[2000],nameP2_2[2000],nameP2_12[2000];
char nameP4_1[2000],nameP4_2[2000],nameP4_12[2000];
char nameB0_1[2000],nameB0_2[2000],nameB0_12[2000];

int Neffmax,NeffmaxSGC;
int NcovP;
double ptheo,ptheoSGC,ptheotot;
double P0tot,errP0tot;
double knl_eff;
double *parameters1,*ktheo,*Ptheo0,*Ptheo2,*Ptheo4;//*Ptheo0win,*Ptheo2win,*Ptheo4win;
double *ktheo0,*ktheo2,*ktheo4;
double *ktheoSGC,*Ptheo0SGC,*Ptheo2SGC,*Ptheo4SGC;//*Ptheo0SGCwin,*Ptheo2SGCwin,*Ptheo4SGCwin;
double *ktheo0SGC,*ktheo2SGC,*ktheo4SGC;
double *b_theo,*b_theoSGC;
int i,j,l,i1;
int dimension;
int modeP0,modeP2,modeP4,modeB0;
int Ntheo=N_Plin;
int offset;
double difference,difference_min;
int ik01_P0,ik01_P2,ik01_P4,ik01_B0;
int combine;
int factor_sampling_mask;
int interpolation_order,shiftN;
double w1,w2,w0;
int Ninterpol;

interpolation_order=1;

if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}

factor_sampling_mask=factor_sampling_mask_in;
if( strcmp(spacing_dataNGC,"irregular") == 0  ){factor_sampling_mask=1;}

combine=0;
dimension=14;
offset=3;//b1,b2,Anoise
if(strcmp(local_b2s2, "no") == 0){offset++;}//b2s2
if(strcmp(local_b3nl, "no") == 0){offset++;}//b3nl
if(strcmp(fog_free, "yes") == 0){offset++;}//fog_ps  or fog_bs
if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_power_spectrum, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 ){offset++;}//fog_bs
if(strcmp(fogmodel_ps, "Exponential_avir") == 0 && strcmp(fog_free, "yes") == 0){offset++;}//avir free

modeP0=0;
modeP2=0;
modeP4=0;
modeB0=0;

if( strcmp(do_power_spectrum,"yes") == 0){

if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}

NcovP=modeP0*NeffP0+modeP2*NeffP2+modeP4*NeffP4;

if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0 )//no mask
{

if(strcmp(path_to_mask1, "none") == 0){
if(NeffP0>0){ktheo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2>0){ktheo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4>0){ktheo4 = (double*) calloc( NeffP4, sizeof(double));}
}

if(strcmp(mask_matrix, "yes") == 0 )
{
Neffmax=get_Neffmax( spacing_dataNGC, modeP0, modeP2, modeP4, NeffP0, NeffP2, NeffP4, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1], Ntheo );
ktheo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
}



if(NeffP0>0){Ptheo0 = (double*) calloc( NeffP0, sizeof(double));}
if(NeffP2>0){Ptheo2 = (double*) calloc( NeffP2, sizeof(double));}
if(NeffP4>0){Ptheo4 = (double*) calloc( NeffP4, sizeof(double));}
}
else
{

Neffmax=get_Neffmax( spacing_dataNGC, modeP0, modeP2, modeP4, NeffP0, NeffP2, NeffP4, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1], Ntheo );
/*
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

*/


ktheo = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));
if(modeP0==1){Ptheo0 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP2==1){Ptheo2 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
if(modeP4==1){Ptheo4 = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}

//if(modeP0==1){Ptheo0win = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
//if(modeP2==1){Ptheo2win = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}
//if(modeP4==1){Ptheo4win = (double*) calloc( Neffmax*factor_sampling_mask, sizeof(double));}

}//else

}//do power

if( strcmp(do_bispectrum,"yes") == 0 )
{
modeB0=1;
b_theo = (double*) calloc( NeffB0, sizeof(double));

}

parameters1 = (double*) calloc( dimension, sizeof(double));
i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
parameters1[i]=parameters2[i1];
i1++;
}

if(strcmp(RSD_fit, "shape") == 0 && i<5){
if(i==3){parameters1[i]=0;}
else{parameters1[i]=parameters2[i1];
i1++;}
}

if(strcmp(RSD_fit, "yes") == 0 && i<5){
if(i==2){parameters1[i]=0;}
if(i==3){parameters1[i]=0;}
if(i!=2 && i!=3){parameters1[i]=parameters2[i1];
i1++;}
}


if(strcmp(RSD_fit, "no") == 0 && i<5){
if(i==0){parameters1[i]=1.;}
if(i==1){parameters1[i]=1.;}
if(i==2){parameters1[i]=0;}
if(i==3){parameters1[i]=0;}
if(i==4){parameters1[i]=0;}
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==5 && strcmp(sigma8_free, "no") == 0 ){
parameters1[i]=Theory[0][41];
}

if(i>=6 && i<=8){
parameters1[i]=parameters2[i1];if(i==8){Anoise_NGC=parameters1[i]/(parameters1[0]*parameters1[1]*parameters1[1]);}
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==9 && strcmp(local_b2s2, "no") != 0){
if(strcmp(local_b2s2, "yes") == 0){parameters1[i]=-4./7.*(parameters1[i-3]-1);}
if(strcmp(local_b2s2, "off") == 0){parameters1[i]=0;}
}

if(i==10 && strcmp(local_b3nl, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==10 && strcmp(local_b3nl, "no") != 0){
if(strcmp(local_b3nl, "yes") == 0){parameters1[i]=32./315.*(parameters1[i-4]-1);}
if(strcmp(local_b3nl, "off") == 0){parameters1[i]=0;}
}

if(i==11 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==11 && strcmp(fog_free, "no") == 0){
parameters1[i]=0;
}

if(i==11 && strcmp(do_power_spectrum, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(fog_free, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(do_bispectrum, "no") == 0){
parameters1[i]=0;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "no") == 0){
parameters1[i]=parameters2[i1];
i1++;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "yes") == 0){//if fog_bs = yes implies power_spectrum = yes
parameters1[i]=parameters1[i-1];
}

if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") != 0){
parameters1[i]=0;
}
if(i==13 && strcmp(fogmodel_ps, "Exponential_avir") == 0){

if( strcmp(fog_free, "yes") == 0 )
{
parameters1[i]=parameters2[i1];
i1++;
}
else{parameters1[i]=0;}

}


}

if( strcmp(do_power_spectrum,"yes") == 0){
do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, ktheo,ktheo0,ktheo2,ktheo4, Ptheo0, Ptheo2, Ptheo4,Pnoise,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC,k0,k2,k4, path_to_mask1,plan1, plan2, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1],spacing_dataNGC,spacing_theory, mask_matrix, MatrixFS_mask_NGC,noise_option);

//if(strcmp(path_to_mask1, "none") != 0 )//no mask
//{

//do_Ptheo_RSD_window(ptmodel_ps,rsdmodel_ps,fogmodel_ps,type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, ktheo,ktheo0,ktheo2,ktheo4, Ptheo0win, Ptheo2win, Ptheo4win,Pnoise,NeffP0,NeffP2,NeffP4,factor_sampling_mask, parameters1,pos, W0, W2, W4, W6, W8, Nmask,spacing_maskNGC,k0,k2,k4, path_to_mask1,plan1, plan2, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1],spacing_dataNGC,spacing_theory,Neffmax,NcovP,interpolation_order,path_output,identifier,1);

//do_Ptheo_RSD_window_wo_model(type_of_analysis, modeP0, modeP2, modeP4, Theory, N_Plin,  ktheo, NeffP0, NeffP2, NeffP4, factor_sampling_mask, pos, W0, W2, W4, W6, W8, Nmask, spacing_maskNGC, k0, k2, k4, path_to_mask1, plan1, plan2, k0[0], k0[NeffP0-1], k2[0], k2[NeffP2-1], k4[0], k4[NeffP4-1], spacing_dataNGC, NcovP, path_output, identifier, 1);
//}

}

if( strcmp(do_bispectrum,"yes") == 0){


if( strcmp(ptmodel_bs,"GilMarin14") ==  0){

if( strcmp(sigma8_free, "no") == 0 ){knl_eff=knl;}
else{knl_eff=P_interpol_extrapol(parameters1[3], sigma8_x, knl_y, Nknl);}//need to reverse the order of sigma8?

}else{knl_eff=0;/*printf("Warning knl_eff=0\n");exit(0);*/}

do_Btheo_RSD(ptmodel_ps,fogmodel_ps,ptmodel_bs,b_theo,NeffB0,Bnoise,parameters1,Theory,N_Plin,pos, W0, W2, W4, W6, W8,Nmask,spacing_maskNGC,path_to_mask1,k11,k22,k33,plan1, plan2,k11[0],k11[NeffB0-1],spacing_theory,knl_eff,n_func_final,bispectrum_BQ,mask_matrix, MatrixFS_mask_NGC,noise_option,redshift_in);

}

free(parameters1);


if(Nchunks==2)
{

if (strcmp(do_power_spectrum,"yes") ==  0){

factor_sampling_mask=factor_sampling_mask_in;

        parameters1 = (double*) calloc( dimension, sizeof(double));

        if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0 )//no mask
        {
                if(strcmp(path_to_mask2, "none") == 0){
                if(NeffP0SGC>0){ktheo0SGC = (double*) calloc( NeffP0SGC, sizeof(double));}
                if(NeffP2SGC>0){ktheo2SGC = (double*) calloc( NeffP2SGC, sizeof(double));}
                if(NeffP4SGC>0){ktheo4SGC = (double*) calloc( NeffP4SGC, sizeof(double));}
                }
                if(strcmp(mask_matrix, "yes") == 0 )
                {
                   NeffmaxSGC=get_Neffmax( spacing_dataSGC, modeP0, modeP2, modeP4, NeffP0SGC, NeffP2SGC, NeffP4SGC, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1], Ntheo );
                ktheoSGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));

                }

               if(NeffP0SGC>0){Ptheo0SGC = (double*) calloc( NeffP0SGC, sizeof(double));}
                if(NeffP2SGC>0){Ptheo2SGC = (double*) calloc( NeffP2SGC, sizeof(double));}
                if(NeffP4SGC>0){Ptheo4SGC = (double*) calloc( NeffP4SGC, sizeof(double));}
        }
       else
        {

NeffmaxSGC=get_Neffmax( spacing_dataSGC, modeP0, modeP2, modeP4, NeffP0SGC, NeffP2SGC, NeffP4SGC, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1], Ntheo );

/*
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
*/

                ktheoSGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));
                if(modeP0==1){Ptheo0SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}
                if(modeP2==1){Ptheo2SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}
                if(modeP4==1){Ptheo4SGC = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}

//                if(modeP0==1){Ptheo0SGCwin = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}
//                if(modeP2==1){Ptheo2SGCwin = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}
//                if(modeP4==1){Ptheo4SGCwin = (double*) calloc( NeffmaxSGC*factor_sampling_mask, sizeof(double));}


        }//else

}//do power

if (strcmp(do_bispectrum,"yes") ==  0){
b_theoSGC = (double*) calloc( NeffB0SGC, sizeof(double));

}


        i1=0;
        for(i=0;i<dimension;i++){

                if(strcmp(RSD_fit, "shape2") == 0 && i<5){
                        parameters1[i]=parameters2[i1];
                        i1++;
                }

                if(strcmp(RSD_fit, "shape") == 0 && i<5){

                        if(i==3){parameters1[i]=0;}
                        else{
                        parameters1[i]=parameters2[i1];
                        i1++;}
                }                

                if(strcmp(RSD_fit, "yes") == 0 && i<5){

                        if(i==2){parameters1[i]=0;}
                        if(i==3){parameters1[i]=0;}
                        if(i!=2 && i!=3){
                        parameters1[i]=parameters2[i1];
                        i1++;}
                }

                if(strcmp(RSD_fit, "no") == 0 && i<5){
                        if(i==0){parameters1[i]=1.;}
                        if(i==1){parameters1[i]=1.;}
                        if(i==2){parameters1[i]=0;}
                        if(i==3){parameters1[i]=0;}
                        if(i==4){parameters1[i]=0;}
                }

                if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
                        parameters1[i]=parameters2[i1];
                        i1++;
                }

                if(i==5 && strcmp(sigma8_free, "no") == 0 ){
                        parameters1[i]=Theory[0][41];
                }

                if(i>=6 && i<=8){
                        parameters1[i]=parameters2[i1+offset];Anoise_SGC=parameters1[i]/(parameters1[0]*parameters1[1]*parameters1[1]);
                        i1++;
                }

                if(i==9 && strcmp(local_b2s2, "no") == 0){
                        parameters1[i]=parameters2[i1+offset];
                        i1++;
                }

                if(i==9 && strcmp(local_b2s2, "no") != 0){
                        if(strcmp(local_b2s2, "yes") == 0){parameters1[i]=-4./7.*(parameters1[i-3]-1);}
                        if(strcmp(local_b2s2, "off") == 0){parameters1[i]=0;}
                }

                if(i==10 && strcmp(local_b3nl, "no") == 0){
              parameters1[i]=parameters2[i1+offset];
                        i1++;
                }

                if(i==10 && strcmp(local_b3nl, "no") != 0){
                        if(strcmp(local_b3nl, "yes") == 0){parameters1[i]=32./315.*(parameters1[i-4]-1);}
                        if(strcmp(local_b3nl, "off") == 0){parameters1[i]=0;}
                }

               if(i==11 && strcmp(fog_free, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){
                        parameters1[i]=parameters2[i1+offset];
                        i1++;
               }

               if(i==11 && strcmp(fog_free, "no") == 0){
                        parameters1[i]=0;
               }

               if(i==11 && strcmp(do_power_spectrum, "no") == 0){
                        parameters1[i]=0;
               }

               if(i==12 && strcmp(fog_free, "no") == 0){
                        parameters1[i]=0;
               }

               if(i==12 && strcmp(do_bispectrum, "no") == 0){
                        parameters1[i]=0;
               }

               if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "no") == 0){
                        parameters1[i]=parameters2[i1+offset];
                        i1++;
               }

               if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs, "yes") == 0){//if fog_bs = yes implies power_spectrum = yes
                        parameters1[i]=parameters1[i-1];
               }


}
if( strcmp(do_power_spectrum,"yes") == 0){
do_Ptheo_RSD(ptmodel_ps,rsdmodel_ps,fogmodel_ps, type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, ktheoSGC,ktheo0SGC,ktheo2SGC,ktheo4SGC, Ptheo0SGC, Ptheo2SGC, Ptheo4SGC,PnoiseSGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC,k0SGC,k2SGC,k4SGC, path_to_mask2,plan1, plan2, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4[NeffP4SGC-1],spacing_dataSGC,spacing_theory, mask_matrix, MatrixFS_mask_SGC,noise_option);

//if(strcmp(path_to_mask2, "none") != 0 )//no mask
//{
//do_Ptheo_RSD_window(ptmodel_ps,rsdmodel_ps,fogmodel_ps,type_of_analysis, RSD_fit, fit_BAO,modeP0,modeP2,modeP4,Theory,N_Plin, ktheoSGC,ktheo0SGC,ktheo2SGC,ktheo4SGC, Ptheo0SGCwin, Ptheo2SGCwin, Ptheo4SGCwin,PnoiseSGC,NeffP0SGC,NeffP2SGC,NeffP4SGC,factor_sampling_mask, parameters1,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC,spacing_maskSGC,k0SGC,k2SGC,k4SGC, path_to_mask2,plan1, plan2, k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1],spacing_dataSGC,spacing_theory,Neffmax,NcovP,interpolation_order,path_output,identifier,2);

//do_Ptheo_RSD_window_wo_model(type_of_analysis, modeP0,modeP2,modeP4, Theory, N_Plin,  ktheoSGC, NeffP0SGC, NeffP2SGC, NeffP4SGC, factor_sampling_mask, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, spacing_maskSGC, k0SGC, k2SGC, k4SGC, path_to_mask2, plan1, plan2,  k0SGC[0], k0SGC[NeffP0SGC-1], k2SGC[0], k2SGC[NeffP2SGC-1], k4SGC[0], k4SGC[NeffP4SGC-1],spacing_dataSGC, NcovP, path_output, identifier, 2);
//}

}

if( strcmp(do_bispectrum,"yes") == 0){

do_Btheo_RSD(ptmodel_ps,fogmodel_ps,ptmodel_bs,b_theoSGC,NeffB0SGC,BnoiseSGC,parameters1,Theory,N_Plin,posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC,NmaskSGC,spacing_maskSGC,path_to_mask2,k11SGC,k22SGC,k33SGC,plan1, plan2,k11SGC[0],k11SGC[NeffB0SGC-1],spacing_theory,knl_eff,n_func_final,bispectrum_BQ,mask_matrix, MatrixFS_mask_SGC,noise_option,redshift_in);

}

free(parameters1);



combine=0;
if(NeffP0!=NeffP0SGC || NeffP2!=NeffP2SGC || NeffP4!=NeffP4SGC || NeffB0!=NeffB0SGC)
{
combine=1;//can't combine
}
else{

if( strcmp(do_power_spectrum,"yes") == 0){
  for(i=0;i<NeffP0;i++){if(k0[i]!=k0SGC[i]){combine++;}}
}

if( strcmp(do_bispectrum,"yes") == 0){
  for(i=0;i<NeffB0;i++){if(k11[i]!=k11SGC[i]){combine++;}}
  for(i=0;i<NeffB0;i++){if(k22[i]!=k22SGC[i]){combine++;}}
  for(i=0;i<NeffB0;i++){if(k33[i]!=k33SGC[i]){combine++;}}
}


}

}//Nchuks

if(combine==0)//decide if they are going to be combined according to Area~1/err^2 for k=0.1 bin
{
if( strcmp(do_power_spectrum,"yes") == 0){

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

if( strcmp(do_bispectrum,"yes") == 0){
   ik01_B0=-1;
   for(i=0;i<NeffB0;i++)
   {
     if(i==0){difference_min=fabs(k11[i]-0.1)+fabs(k22[i]-0.1)+fabs(k33[i]-0.1);ik01_B0=i;}
     else
     {
        difference=fabs(k11[i]-0.1)+fabs(k22[i]-0.1)+fabs(k33[i]-0.1);
        if(difference<difference_min){difference_min=difference;ik01_B0=i;}

     }
  }


}

}//combine 0

sprintf(nameP0_1,"%s/Monopole1_FS_%s.txt",path_output,identifier);
sprintf(nameP0_2,"%s/Monopole2_FS_%s.txt",path_output,identifier);
sprintf(nameP0_12,"%s/Monopole12_FS_%s.txt",path_output,identifier);

sprintf(nameP2_1,"%s/Quadrupole1_FS_%s.txt",path_output,identifier);
sprintf(nameP2_2,"%s/Quadrupole2_FS_%s.txt",path_output,identifier);
sprintf(nameP2_12,"%s/Quadrupole12_FS_%s.txt",path_output,identifier);

sprintf(nameP4_1,"%s/Hexadecapole1_FS_%s.txt",path_output,identifier);
sprintf(nameP4_2,"%s/Hexadecapole2_FS_%s.txt",path_output,identifier);
sprintf(nameP4_12,"%s/Hexadecapole12_FS_%s.txt",path_output,identifier);

sprintf(nameB0_1,"%s/Bispectrum1_FS_%s.txt",path_output,identifier);
sprintf(nameB0_2,"%s/Bispectrum2_FS_%s.txt",path_output,identifier);
sprintf(nameB0_12,"%s/Bispectrum12_FS_%s.txt",path_output,identifier);

for(l=0;l<=6;l=l+2)
{
open=0;
if(modeP0==1 && l==0){string_copy(nameP0_1,name1);string_copy(nameP0_2,name2);string_copy(nameP0_12,name12);open=1;}
if(modeP2==1 && l==2){string_copy(nameP2_1,name1);string_copy(nameP2_2,name2);string_copy(nameP2_12,name12);open=1;}
if(modeP4==1 && l==4){string_copy(nameP4_1,name1);string_copy(nameP4_2,name2);string_copy(nameP4_12,name12);open=1;}
if(modeB0==1 && l==6){string_copy(nameB0_1,name1);string_copy(nameB0_2,name2);string_copy(nameB0_12,name12);open=1;}

if(open==1){

f1=fopen(name1,"w");
if(Nchunks==2)
{
f2=fopen(name2,"w");
if(combine==0){f12=fopen(name12,"w");}
}

if(l<=4){fprintf(f1,"#k \t Pdata \t err \t Pmodel\n");}
if(l==6){fprintf(f1,"#k1 \t k2 \t k3 \t Bdata \t err \t Bmodel\n");}


if(Nchunks==2)
{
if(l<=4){fprintf(f2,"#k \t Pdata \t err \t Pmodel\n");}
if(l==6){fprintf(f2,"#k1 \t k2 \t k3 \t Bdata \t err \t Bmodel\n");}
if(combine==0){
      if(l<=4){fprintf(f12,"#k \t Pdata \t err \t Pmodel\n");}
      if(l==6){fprintf(f12,"#k1 \t k2 \t k3 \t Bdata \t err \t Bmodel\n");}
      }
}

if(l==0){NeffP=NeffP0;}
if(l==2){NeffP=NeffP2;}
if(l==4){NeffP=NeffP4;}
if(l==6){NeffP=NeffB0;}

for(i=0;i<NeffP;i++)
{

if(l==6)
{

ptheo=b_theo[i]-Bnoise[i];
fprintf(f1,"%e %e %e %e %e %e\n",k11[i],k22[i],k33[i],B0[i],errB0[i],ptheo);

}

if(l==0)
{
if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0){ptheo=Ptheo0[i]-Pnoise+Anoise_NGC*Pnoise*noise_option;}
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
ptheo=P_interpol_fast(k0[i],Ptheo0,Neffmax*factor_sampling_mask,spacing_dataNGC,interpolation_order,Ninterpol,w0,w1,w2)-Pnoise+Anoise_NGC*Pnoise*noise_option;
}

}
fprintf(f1,"%e %e %e %e\n",k0[i],P0[i],errP0[i],ptheo);
}

if(l==2)
{
if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0 ){ptheo=Ptheo2[i];}
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

}

fprintf(f1,"%e %e %e %e\n",k2[i],P2[i],errP2[i],ptheo);
}

if(l==4)
{
if(strcmp(path_to_mask1, "none") == 0 || strcmp(mask_matrix, "yes") == 0 ){ptheo=Ptheo4[i];}
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

}
fprintf(f1,"%e %e %e %e\n",k4[i],P4[i],errP4[i],ptheo);
}



if(Nchunks==2)
{

if(combine==0)
{

if(l==6)
{
ptheoSGC=b_theoSGC[i]-BnoiseSGC[i];
}

if(l==0)
{
if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0){ptheoSGC=Ptheo0SGC[i]-PnoiseSGC+Anoise_SGC*PnoiseSGC*noise_option;}
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
ptheoSGC=P_interpol_fast(k0SGC[i],Ptheo0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2)-PnoiseSGC+Anoise_SGC*PnoiseSGC*noise_option;
}


}
}

if(l==2)
{
if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0){ptheoSGC=Ptheo2SGC[i];}
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


}
}

if(l==4)
{
if(strcmp(path_to_mask2, "none") == 0 || strcmp(mask_matrix, "yes") == 0){ptheoSGC=Ptheo4SGC[i];}
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


}
}

if(l==6)
{

P0tot=(B0[i]*pow(errB0[ik01_B0],-2)+B0SGC[i]*pow(errB0SGC[ik01_B0],-2))/( pow(errB0[ik01_B0],-2)+pow(errB0SGC[ik01_B0],-2) );

ptheotot=(ptheo*pow(errB0[ik01_B0],-2)+ptheoSGC*pow(errB0SGC[ik01_B0],-2))/( pow(errB0[ik01_B0],-2)+pow(errB0SGC[ik01_B0],-2) );

errP0tot=errB0[i]*sqrt( pow(errB0[ik01_B0],-2)/(pow(errB0[ik01_B0],-2)+pow(errB0SGC[ik01_B0],-2)) );//errNS=errNGC*sqrt(ANGC/(ANGC+ASGC))

fprintf(f12,"%e %e %e %e %e %e\n",k11[i],k22[i],k33[i],P0tot,errP0tot,ptheotot);

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
if(l==6){NeffP=NeffB0SGC;}

for(i=0;i<NeffP;i++)
{

if(l==6)
{
ptheoSGC=b_theoSGC[i]-BnoiseSGC[i];
fprintf(f2,"%e %e %e %e %e %e\n",k11SGC[i],k22SGC[i],k33SGC[i],B0SGC[i],errB0SGC[i],ptheoSGC);
}

if(l==0)
{
if(strcmp(path_to_mask2, "none") != 0 && strcmp(mask_matrix, "no") == 0){

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
ptheoSGC=P_interpol_fast(k0SGC[i],Ptheo0SGC,NeffmaxSGC*factor_sampling_mask,spacing_dataSGC,interpolation_order,Ninterpol,w0,w1,w2)-PnoiseSGC+Anoise_SGC*PnoiseSGC*noise_option;
}


}
else{ptheoSGC=Ptheo0SGC[i]-PnoiseSGC+Anoise_SGC*PnoiseSGC*noise_option;}
fprintf(f2,"%e %e %e %e\n",k0SGC[i],P0SGC[i],errP0SGC[i],ptheoSGC);
}

if(l==2)
{
if(strcmp(path_to_mask2, "none") != 0 && strcmp(mask_matrix, "no") == 0 ){

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


}
else{ptheoSGC=Ptheo2SGC[i];}
fprintf(f2,"%e %e %e %e\n",k2SGC[i],P2SGC[i],errP2SGC[i],ptheoSGC);
}

if(l==4)
{
if(strcmp(path_to_mask2, "none") != 0 && strcmp(mask_matrix, "no") == 0 ){

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
}
else{ptheoSGC=Ptheo4SGC[i];}
fprintf(f2,"%e %e %e %e\n",k4SGC[i],P4SGC[i],errP4SGC[i],ptheoSGC);
}


}//for NeffP

}//if Nchunks=2

fprintf(f1,"#");
if(strcmp(RSD_fit, "yes") == 0){fprintf(f1,"apara \t aperp \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "shape") == 0){fprintf(f1,"apara \t aperp \t mBGV \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "shape2") == 0){fprintf(f1,"apara \t aperp \t m1BGV \t m2BGV \t f \t");}//f,apara,aperp as free parameter
if(strcmp(sigma8_free, "yes") == 0){fprintf(f1,"s8 \t");}
fprintf(f1,"b1_NGC \t b2_NGC \t Anoise_NGC\t");
if(strcmp(local_b2s2, "yes") != 0){fprintf(f1,"b2s2_NGC \t");}//b2s2
if(strcmp(local_b3nl, "yes") != 0){fprintf(f1,"b3nl_NGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f1,"sigmaP_NGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){fprintf(f1,"sigmaB_NGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp( fogmodel_ps , "Exponential_avir") == 0){ fprintf(f1,"avir_NGC\t");}
fprintf(f1,"chi2 \t Npoints-Ndof\n");

i1=0;
fprintf(f1,"#");
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(strcmp(RSD_fit, "shape") == 0 && i<4){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(i>=6 && i<=8){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(i==9 && strcmp(local_b2s2, "yes") != 0){

if(strcmp(local_b2s2, "no") == 0){fprintf(f1,"%lf\t",parameters2[i1]);i1++;}
if(strcmp(local_b2s2, "off") == 0){fprintf(f1,"0\t");}
}

if(i==10 && strcmp(local_b3nl, "yes") != 0){
if(strcmp(local_b3nl, "no") == 0){fprintf(f1,"%lf\t",parameters2[i1]);i1++;}
if(strcmp(local_b3nl, "off") == 0){fprintf(f1,"0\t");}
}

if(i==11 && strcmp(fog_free, "yes") == 0){//can be either P or B (if only one of them is selected: eg. do_power_spectrum = no, do_bispectrum = yes
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}

if( i==13 && strcmp(fogmodel_ps,"Exponential_avir") == 0 && strcmp(fog_free, "yes") == 0)
{
fprintf(f1,"%lf\t",parameters2[i1]);
i1++;
}


}


fprintf(f1,"%lf \t %d-%d\n",chi2_min,Npoints,Ndof);

if(Nchunks==2)
{

fprintf(f2,"#");
if(strcmp(RSD_fit, "yes") == 0){fprintf(f2,"apara \t aperp \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "shape") == 0){fprintf(f2,"apara \t aperp \t mBGV \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "shape2") == 0){fprintf(f2,"apara \t aperp \t m1BGV \t m2BGV \t f \t");}//f,apara,aperp as free parameter
if(strcmp(sigma8_free, "yes") == 0){fprintf(f2,"s8 \t");}
fprintf(f2,"b1_SGC \t b2_SGC \t Anoise_SGC\t");
if(strcmp(local_b2s2, "yes") != 0){fprintf(f2,"b2s2_SGC \t");}//b2s2
if(strcmp(local_b3nl, "yes") != 0){fprintf(f2,"b3nl_SGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f2,"sigmaP_SGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){fprintf(f2,"sigmaB_SGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp( fogmodel_ps , "Exponential_avir") == 0){ fprintf(f2,"avir_SGC\t");}
fprintf(f2,"chi2 \t Npoints-Ndof\n");



i1=0;
fprintf(f2,"#");
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
fprintf(f2,"%lf\t",parameters2[i1]);
i1++;
}

if(strcmp(RSD_fit, "shape") == 0 && i<4){
fprintf(f2,"%lf\t",parameters2[i1]);
i1++;
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
fprintf(f2,"%lf\t",parameters2[i1]);
i1++;
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
fprintf(f2,"%lf\t",parameters2[i1]);
i1++;
}
if(i>=6 && i<=8){
fprintf(f2,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==9 && strcmp(local_b2s2, "yes") != 0){
if(strcmp(local_b2s2, "no") == 0){fprintf(f2,"%lf\t",parameters2[i1+offset]);i1++;}
if(strcmp(local_b2s2, "off") == 0){fprintf(f2,"0\t");}
}

if(i==10 && strcmp(local_b3nl, "yes") != 0){
if(strcmp(local_b3nl, "no") == 0){fprintf(f2,"%lf\t",parameters2[i1+offset]);i1++;}
if(strcmp(local_b3nl, "off") == 0){fprintf(f2,"0\t");}
}

if(i==11 && strcmp(fog_free, "yes") == 0){//can be either P or B (if only one of them is selected: eg. do_power_spectrum = no, do_bispectrum = yes
fprintf(f2,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){
fprintf(f2,"%lf\t",parameters2[i1+offset]);
i1++;
}

if( i==13 && strcmp(fogmodel_ps,"Exponential_avir") == 0 && strcmp(fog_free, "yes") == 0)
{
fprintf(f2,"%lf\t",parameters2[i1+offset]);
i1++;
}

}
fprintf(f2,"%lf \t %d-%d\n",chi2_min,Npoints,Ndof);




if(combine==0){


fprintf(f12,"#");
if(strcmp(RSD_fit, "yes") == 0){fprintf(f12,"apara \t aperp \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "shape") == 0){fprintf(f12,"apara \t aperp \t mBGV \t f \t");}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "shape2") == 0){fprintf(f12,"apara \t aperp \t m1BGV \t m2BGV \t f \t");}//f,apara,aperp as free parameter
if(strcmp(sigma8_free, "yes") == 0){fprintf(f12,"s8 \t");}
fprintf(f12,"b1_NGC \t b2_NGC \t Anoise_NGC\t");
if(strcmp(local_b2s2, "yes") != 0){fprintf(f12,"b2s2_NGC \t");}//b2s2
if(strcmp(local_b3nl, "yes") != 0){fprintf(f12,"b3nl_NGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f12,"sigmaP_NGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){fprintf(f12,"sigmaB_NGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp( fogmodel_ps , "Exponential_avir") == 0){ fprintf(f12,"avir_NGC\t");}
fprintf(f12,"b1_SGC \t b2_SGC \t Anoise_SGC\t");
if(strcmp(local_b2s2, "no") == 0){fprintf(f12,"b2s2_SGC \t");}//b2s2
if(strcmp(local_b3nl, "no") == 0){fprintf(f12,"b3nl_SGC \t");}//b3nl
if(strcmp(fog_free, "yes") == 0  && strcmp(do_power_spectrum, "yes") == 0){fprintf(f12,"sigmaP_SGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp(do_bispectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){fprintf(f12,"sigmaB_SGC\t");}
if(strcmp(fog_free, "yes") == 0  && strcmp( fogmodel_ps , "Exponential_avir") == 0){ fprintf(f2,"avir_SGC\t");}
fprintf(f12,"chi2 \t Npoints-Ndof\n");

i1=0;
fprintf(f12,"#");
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if(strcmp(RSD_fit, "shape") == 0 && i<4){
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}


if(i>=6 && i<=8){
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if(i==9 && strcmp(local_b2s2, "yes") != 0){
if(strcmp(local_b2s2, "no") == 0){fprintf(f12,"%lf\t",parameters2[i1]);i1++;}
if(strcmp(local_b2s2, "off") == 0){fprintf(f12,"0\t");}
}

if(i==10 && strcmp(local_b3nl, "yes") != 0){
if(strcmp(local_b3nl, "no") == 0){fprintf(f12,"%lf\t",parameters2[i1]);i1++;}
if(strcmp(local_b3nl, "off") == 0){fprintf(f12,"0\t");}
}

if(i==11 && strcmp(fog_free, "yes") == 0){//can be either P or B (if only one of them is selected: eg. do_power_spectrum = no, do_bispectrum = yes
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

if( i==13 && strcmp(fogmodel_ps,"Exponential_avir") == 0 && strcmp(fog_free, "yes") == 0)
{
fprintf(f12,"%lf\t",parameters2[i1]);
i1++;
}

}

i1=0;
for(i=0;i<dimension;i++){

if(strcmp(RSD_fit, "yes") == 0 && i<3){
i1++;
}
if(strcmp(RSD_fit, "shape") == 0 && i<4){
i1++;
}

if(strcmp(RSD_fit, "shape2") == 0 && i<5){
i1++;
}

if(i==5 && strcmp(sigma8_free, "yes") == 0 ){
i1++;
}

if(i>=6 && i<=8){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==9 && strcmp(local_b2s2, "yes") != 0){
if(strcmp(local_b2s2, "no") == 0){fprintf(f12,"%lf\t",parameters2[i1+offset]);i1++;}
if(strcmp(local_b2s2, "off") == 0){fprintf(f12,"0\t");}
}

if(i==10 && strcmp(local_b3nl, "yes") != 0){
if(strcmp(local_b3nl, "no") == 0){fprintf(f12,"%lf\t",parameters2[i1+offset]);i1++;}
if(strcmp(local_b3nl, "off") == 0){fprintf(f12,"0\t");}
}


if(i==11 && strcmp(fog_free, "yes") == 0){//can be either P or B (if only one of them is selected: eg. do_power_spectrum = no, do_bispectrum = yes
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

if(i==12 && strcmp(fog_free, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 && strcmp(fog_bs,"no") == 0){
fprintf(f12,"%lf\t",parameters2[i1+offset]);
i1++;
}

if( i==13 && strcmp(fogmodel_ps,"Exponential_avir") == 0 && strcmp(fog_free, "yes") == 0)
{
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
//else{free(Ptheo0win);}
free(Ptheo0);}
if(modeP2==1){
if(strcmp(path_to_mask1, "none") == 0){free(ktheo2);}
//else{free(Ptheo2win);}
free(Ptheo2);}
if(modeP4==1){
if(strcmp(path_to_mask1, "none") == 0){free(ktheo4);}
//else{free(Ptheo4win);}
free(Ptheo4);}

if(modeB0==1){
free(b_theo);
}

if(modeP0 == 1 || modeP2 == 1 || modeP4 == 1){
if(strcmp(path_to_mask1, "none") != 0  ){
free(ktheo);
}}

if(Nchunks==2)
{
if(modeP0==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo0SGC);}
//else{free(Ptheo0SGCwin);}
free(Ptheo0SGC);}
if(modeP2==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo2SGC);}
//else{free(Ptheo2SGCwin);}
free(Ptheo2SGC);}
if(modeP4==1){
if(strcmp(path_to_mask2, "none") == 0 ){free(ktheo4SGC);}
//else{free(Ptheo4SGCwin);}
free(Ptheo4SGC);}

if(modeP0 == 1 || modeP2 == 1 || modeP4 == 1){
if(strcmp(path_to_mask2, "none") != 0 ){
free(ktheoSGC);
}}

if(modeB0==1){
free(b_theoSGC);
}

}


}

void do_rsd_bao_mcmc(int nthreads,char *type_BAO_fit,char *type_of_analysis,char *fit_BAO,char *fit_RSD,double **Theory,int Ntheory,double *k_Plin,double *Plin,int N_Plin, double *k_Olin, double *Olin, int N_Olin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC, char *path_to_mask2, char *spacing_maskSGC, double *k0bao, double *k0rsd, double *P0bao, double *P0rsd, double *errP0bao, double *errP0rsd,double Pnoise_rsd, int NeffP0bao, int NeffP0rsd, double *k2bao, double *k2rsd, double *P2bao, double *P2rsd, double *errP2bao, double *errP2rsd, int NeffP2bao, int NeffP2rsd, double *k4bao, double *k4rsd, double *P4bao, double *P4rsd, double *errP4bao, double *errP4rsd, int NeffP4bao, int NeffP4rsd, double *k11bao, double *k11rsd, double *k22bao, double *k22rsd, double *k33bao, double *k33rsd, double *B0bao, double *B0rsd, double *errB0bao, double *errB0rsd, double *Bnoise_bao, double *Bnoise_rsd, int NeffB0bao, int NeffB0rsd, double *k0baoSGC, double *k0rsdSGC, double *P0baoSGC, double *P0rsdSGC, double *errP0baoSGC, double *errP0rsdSGC,double Pnoise_rsdSGC, int NeffP0baoSGC,int NeffP0rsdSGC, double *k2baoSGC, double *k2rsdSGC, double *P2baoSGC, double *P2rsdSGC, double *errP2baoSGC, double *errP2rsdSGC, int NeffP2baoSGC, int NeffP2rsdSGC, double *k4baoSGC, double *k4rsdSGC, double *P4baoSGC, double *P4rsdSGC, double *errP4baoSGC, double *errP4rsdSGC, int NeffP4baoSGC, int NeffP4rsdSGC, double *k11baoSGC,double *k11rsdSGC, double *k22baoSGC,double *k22rsdSGC, double *k33baoSGC,double *k33rsdSGC,double *B0baoSGC,double *B0rsdSGC, double *errB0baoSGC, double *errB0rsdSGC, double *Bnoise_baoSGC, double *Bnoise_rsdSGC,int NeffB0baoSGC,int NeffB0rsdSGC, double *cov, double *covSGC, double alpha_min, double alpha_max, char *Sigma_def_type, char *Sigma_independent, double ffactor, double Sigma_type[], double Sigma_nl_mean[], double Sigma_nl_stddev[], int  Npolynomial, char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit,char *sigma8_free,char *fog_free,char *fog_bs, int  Nchunks, char *path_output, char *identifier, char *do_plot, char *use_prop_cov, char *path_to_cov, int Nsteps, char *do_power_spectrum, char *do_bispectrum, double Sigma_smooth, char *spacing_dataNGC_bao, char *spacing_dataNGC_rsd, char *spacing_dataSGC_bao, char *spacing_dataSGC_rsd, char *spacing_theory,char *spacing_theory_rsd,char *type_of_analysis_BAO, char *type_of_analysis_FS,char *bispectrum_BQ, int factor_sampling_mask, char *mask_matrix, double **MatrixBAO_mask_NGC, double **MatrixBAO_mask_SGC, double **MatrixFS_mask_NGC, double **MatrixFS_mask_SGC, double *FSprior_type, double *FSprior_mean, double *FSprior_stddev,int noise_option, double step_size,char *covariance_correction, int NrealNGC, int NrealSGC)
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
int i_thread;

int Nplanrsd,Nplanbao;
double chi2;
double params[5];
int Ndof, Npoints, NpointsSGC;
int Npointsbao,Npointsrsd;
int NpointsbaoSGC,NpointsrsdSGC;
fftw_complex *a_pointerbao;a_pointerbao=NULL;
fftw_complex *b_pointerbao;b_pointerbao=NULL;
fftw_complex *a_pointerrsd;a_pointerrsd=NULL;
fftw_complex *b_pointerrsd;b_pointerrsd=NULL;
fftw_plan plan1bao,plan2bao,plan1rsd,plan2rsd;

int i1,i2,i;
int modeB0rsd,modeB0bao;
int dimension,dimensionbao,dimensionrsd;
int dimensionP;
double *parameters2;
double *parameters2_bao,*parameters2_rsd;
int baoiso_shift;

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
if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0){N_Cov_prop_rsd++;}
if(strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 ){N_Cov_prop_rsd++;}
N_Cov_prop_rsd=N_Cov_prop_rsd*Nchunks;//x2 
if(strcmp(sigma8_free, "yes") == 0){N_Cov_prop_rsd++;}
if(strcmp(RSD_fit, "yes") == 0){N_Cov_prop_rsd=N_Cov_prop_rsd+3;}//f, alpha-para, alpha-perp
if(strcmp(RSD_fit, "shape") == 0){N_Cov_prop_rsd=N_Cov_prop_rsd+4;}//f, alpha-para, alpha-perp, m_BGV
if(strcmp(RSD_fit, "shape2") == 0){N_Cov_prop_rsd=N_Cov_prop_rsd+5;}//f, alpha-para, alpha-perp, m1_BGV, m2_BGV

N_Cov_prop=N_Cov_prop_rsd+N_Cov_prop_bao-2;
//printf("%d,%d,%d %s\n",N_Cov_prop,N_Cov_prop_bao,N_Cov_prop_rsd,use_prop_cov);
//exit(0);
//read prop. cov.

if(Nsteps<=0)//in this case we just make  a plot
{

printf("Warning, invalid number of steps selected for the mcmc (%d). Only ploting function at prior set of parameters.\n",Nsteps);

vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized
if(strcmp(use_prop_cov, "yes") == 0)//full mcmc run with proposal
{
trial_mcmc=0;//no trial
Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
read_prop_cov(NULL,0,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);
free(Cov_prop);
}
if(strcmp(use_prop_cov, "no") == 0)//from set_propsal_mean
{
set_proposal_mean(vector_mean, N_Cov_prop, type_of_analysis, fit_BAO, do_power_spectrum, do_bispectrum, Sigma_def_type, Sigma_independent, Sigma_type, local_b2s2, local_b3nl, fog_free, fogmodel_ps,fog_bs, sigma8_free, RSD_fit, Nchunks, Npolynomial, path_output, identifier);
}

set_mask_params(params,0.0,1.0,0.0,1.0,0.0);

Nplanrsd=(int)(params[2]);
Nplanbao=(int)(params[2]);

    plan1bao = fftw_plan_dft_1d(Nplanbao,  a_pointerbao,  b_pointerbao,  -1, FFTW_ESTIMATE);//forward plan
    plan2bao = fftw_plan_dft_1d(Nplanbao,  b_pointerbao,  b_pointerbao, +1, FFTW_ESTIMATE);//reverse plan

    plan1rsd = fftw_plan_dft_1d(Nplanrsd,  a_pointerrsd,  b_pointerrsd,  -1, FFTW_ESTIMATE);//forward plan
    plan2rsd = fftw_plan_dft_1d(Nplanrsd,  b_pointerrsd,  b_pointerrsd, +1, FFTW_ESTIMATE);//reverse plan

Ndof=N_Cov_prop;

//define modeP0, modeP2,modeP4,modeB0
modeP0bao=0;modeP2bao=0;modeP4bao=0;modeB0bao=0;
modeP0rsd=0;modeP2rsd=0;modeP4rsd=0;modeB0rsd=0;
Npointsbao=0;
Npointsrsd=0;
if( strcmp(do_power_spectrum,"yes") == 0 ){

if(strcmp(type_of_analysis_BAO,"yes")==0){
if(strcmp(fit_BAO, "P0") == 0){Npointsbao=NeffP0bao;modeP0bao=1;}
if(strcmp(fit_BAO, "P2") == 0){Npointsbao=NeffP2bao;modeP2bao=1;}
if(strcmp(fit_BAO, "P4") == 0){Npointsbao=NeffP4bao;modeP4bao=1;}
if(strcmp(fit_BAO, "P02") == 0){Npointsbao=NeffP0bao+NeffP2bao;modeP0bao=1;modeP2bao=1;}
if(strcmp(fit_BAO, "P04") == 0){Npointsbao=NeffP0bao+NeffP4bao;modeP0bao=1;modeP4bao=1;}
if(strcmp(fit_BAO, "P24") == 0){Npointsbao=NeffP2bao+NeffP4bao;modeP2bao=1;modeP4bao=1;}
if(strcmp(fit_BAO, "P024") == 0){Npointsbao=NeffP0bao+NeffP2bao+NeffP4bao;modeP0bao=1;modeP2bao=1;modeP4bao=1;}}

if(strcmp(type_of_analysis_FS,"yes")==0){
if(strcmp(fit_RSD, "P0") == 0){Npointsrsd=NeffP0rsd;modeP0rsd=1;}
if(strcmp(fit_RSD, "P2") == 0){Npointsrsd=NeffP2rsd;modeP2rsd=1;}
if(strcmp(fit_RSD, "P4") == 0){Npointsrsd=NeffP4rsd;modeP4rsd=1;}
if(strcmp(fit_RSD, "P02") == 0){Npointsrsd=NeffP0rsd+NeffP2rsd;modeP0rsd=1;modeP2rsd=1;}
if(strcmp(fit_RSD, "P04") == 0){Npointsrsd=NeffP0rsd+NeffP4rsd;modeP0rsd=1;modeP4rsd=1;}
if(strcmp(fit_RSD, "P24") == 0){Npointsrsd=NeffP2rsd+NeffP4rsd;modeP2rsd=1;modeP4rsd=1;}
if(strcmp(fit_RSD, "P024") == 0){Npointsrsd=NeffP0rsd+NeffP2rsd+NeffP4rsd;modeP0rsd=1;modeP2rsd=1;modeP4rsd=1;}}

}

if( strcmp(do_bispectrum,"yes") == 0 ){

    if(strcmp(type_of_analysis_BAO,"yes")==0){
    Npointsbao=Npointsbao+NeffB0bao;modeB0bao=1;//printf("Warning!! %s %s\n",do_bispectrum,type_of_analysis_BAO);
    }
    if(strcmp(type_of_analysis_FS,"yes")==0){
    Npointsrsd=Npointsrsd+NeffB0rsd;modeB0rsd=1;
    }

}

Npoints=Npointsrsd+Npointsbao;
NpointsSGC=0;
if(Nchunks==2)
{
NpointsbaoSGC=0;
if( strcmp(do_power_spectrum,"yes") == 0 ){
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
}
if( strcmp(do_bispectrum,"yes") == 0 ){

   if(strcmp(type_of_analysis_BAO,"yes")==0){
  NpointsbaoSGC=NpointsbaoSGC+NeffB0baoSGC;
   }
   if(strcmp(type_of_analysis_FS,"yes")==0){
  NpointsrsdSGC=NpointsrsdSGC+NeffB0rsdSGC;
   }

}
NpointsSGC=NpointsbaoSGC+NpointsrsdSGC;
}


if(strcmp(type_of_analysis, "FSBAOISO") == 0){
dimension=0;
dimensionbao=0;
dimensionP=0;
if( strcmp(do_power_spectrum,"yes") == 0){
dimension=Nchunks*(Npolynomial+1)*(modeP0bao+modeP2bao+modeP4bao)+Nalphas+Nsigmas_tot;//for parameters plot
dimensionbao=Nchunks*(Npolynomial+1)*(modeP0bao+modeP2bao+modeP4bao)+Nalphas+Nsigmas_tot;//for parameters plot
dimensionP=dimension;
}
if( strcmp(do_bispectrum,"yes") == 0){
dimension=dimension+( Nchunks*(Npolynomial+1)+6.*Nchunks+1.0);//+1 is for sigma_B
dimensionbao=dimensionbao+( Nchunks*(Npolynomial+1)+6.*Nchunks+1.0);//+1 is for sigmaB

if( strcmp(do_power_spectrum,"no") == 0){
dimension++;//add one more for alpha0
dimensionbao++;//add one more for alpha0
}

}

if( strcmp(type_of_analysis, "FSBAOISO") == 0 && Nalphas==1){dimension=dimension+1;}//we add space for an extra alpha when FS-BAO and P0,P2,P4

//TBD for Bispectrum, not available so far

}//BAOISO o FSBAOISO

if(strcmp(type_of_analysis, "FSBAOANISO") == 0){

dimension=0;
dimensionbao=0;
dimensionP=0;
if( strcmp(do_power_spectrum,"yes") == 0){
dimension=Nchunks*(1+Npolynomial*(modeP0bao+modeP2bao+modeP4bao))+Nalphas+Nsigmas_tot+1;//for parameters plot (+1 is for beta)
dimensionbao=Nchunks*(1+Npolynomial*(modeP0bao+modeP2bao+modeP4bao))+Nalphas+Nsigmas_tot+1;//for parameters plot (+1 is for beta)
dimensionP=dimension;
}

if( strcmp(do_bispectrum,"yes") == 0){
dimension=dimension+( Nchunks*(Npolynomial+1)+6.*Nchunks+1.0);//+1 is for sigma_B
dimensionbao=dimensionbao+( Nchunks*(Npolynomial+1)+6.*Nchunks+1.0);//+1 is for sigma_B
if( strcmp(do_power_spectrum,"no") == 0){
dimension++;//add one more for alpha0
dimensionbao++;//add one more for alpha0
}

}

}//BAOANISO o FSBAOANISO


if(strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){

dimension=3;//b1,b2,A
if(strcmp(local_b2s2, "no") == 0){dimension++;}//b2s2
if(strcmp(local_b3nl, "no") == 0){dimension++;}//b3nl
if(strcmp(fog_free, "yes") == 0){dimension++;}//fog_ps
if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0){dimension++;}//fog_bs
dimension=dimension*Nchunks;//x2
if(strcmp(sigma8_free, "yes") == 0){dimension++;}//s8 as free parameter
if(strcmp(RSD_fit, "yes") == 0){dimension=dimension+3;}//f,apara,aperp as free parameter
if(strcmp(RSD_fit, "shape") == 0){dimension=dimension+4;}//f,apara,aperp,m_BGV as free parameter
if(strcmp(RSD_fit, "shape2") == 0){dimension=dimension+5;}//f,apara,aperp,m1_BGV, m2_BGV as free parameter
dimensionrsd=dimension;
}//FS o FSBAOISO o FSBAOANISO

if(strcmp(type_of_analysis, "FSBAOISO") == 0 || strcmp(type_of_analysis, "FSBAOANISO") == 0){dimension=dimensionrsd+dimensionbao-2;}


parameters2 =  (double*) calloc( dimension, sizeof(double));//printf("dim=%d\n", dimension);
if(strcmp(type_of_analysis_BAO,"yes")==0 && strcmp(type_of_analysis_FS,"yes")==0)
{
parameters2_bao =  (double*) calloc( dimensionbao, sizeof(double));
parameters2_rsd =  (double*) calloc( dimensionrsd, sizeof(double));
}

if(strcmp(type_of_analysis, "FSBAOISO") == 0){

parameters2[0]=vector_mean[0];//apara
parameters2[1]=vector_mean[1];//aperp

if(strcmp(RSD_fit, "shape") == 0){
parameters2[2]=vector_mean[2];//m_BGV
parameters2[3]=vector_mean[3];//f
i2=4;
}

if(strcmp(RSD_fit, "shape2") == 0){
parameters2[2]=vector_mean[2];//m1_BGV
parameters2[3]=vector_mean[3];//m2_BGV
parameters2[4]=vector_mean[4];//f
i2=5;
}

if(strcmp(RSD_fit, "yes") == 0){
parameters2[2]=vector_mean[2];//f
i2=3;
}

if(modeP0bao+modeP2bao+modeP4bao>1){
parameters2_bao[0]=vector_mean[0];//apara
parameters2_bao[1]=vector_mean[1];//aperp
baoiso_shift=1;
}
else
{
baoiso_shift=2;
if(modeP0bao==1){parameters2_bao[0]=pow(vector_mean[0],1./3.)*pow(vector_mean[1],2./3.);  }
if(modeP2bao==1){parameters2_bao[0]=pow(vector_mean[0],3./5.)*pow(vector_mean[1],2./5.);  }
if(modeP4bao==1){parameters2_bao[0]=pow(vector_mean[0],5./7.)*pow(vector_mean[1],2./7.);  }
}

i1=2;//i2=3;
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

for(i=i1+1;i<dimensionbao+baoiso_shift;i++){parameters2[i]=vector_mean[i2];parameters2_bao[i-baoiso_shift]=vector_mean[i2];i2++;}//set the rest of parameters (check we need dimension+1)

if(strcmp(RSD_fit, "yes") == 0){
parameters2_rsd[0]=vector_mean[0];
parameters2_rsd[1]=vector_mean[1];
parameters2_rsd[2]=vector_mean[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=vector_mean[i2];parameters2[i-3+dimensionbao+1]=vector_mean[i2];i2++;}
}

if(strcmp(RSD_fit, "shape") == 0){
parameters2_rsd[0]=vector_mean[0];
parameters2_rsd[1]=vector_mean[1];
parameters2_rsd[2]=vector_mean[2];
parameters2_rsd[3]=vector_mean[3];
for(i=4;i<dimensionrsd;i++){parameters2_rsd[i]=vector_mean[i2];parameters2[i-4+dimensionbao+1]=vector_mean[i2];i2++;}
}

if(strcmp(RSD_fit, "shape2") == 0){
parameters2_rsd[0]=vector_mean[0];
parameters2_rsd[1]=vector_mean[1];
parameters2_rsd[2]=vector_mean[2];
parameters2_rsd[3]=vector_mean[3];
parameters2_rsd[4]=vector_mean[4];
for(i=5;i<dimensionrsd;i++){parameters2_rsd[i]=vector_mean[i2];parameters2[i-5+dimensionbao+1]=vector_mean[i2];i2++;}
}

}

if(strcmp(type_of_analysis, "FSBAOANISO") == 0){

parameters2[0]=vector_mean[0];//apara
parameters2[1]=vector_mean[1];//aperp

if(strcmp(RSD_fit, "yes") == 0){
parameters2[2]=vector_mean[2];//f
i2=3;
}
if(strcmp(RSD_fit, "shape") == 0){
parameters2[2]=vector_mean[2];//m_BGV
parameters2[3]=vector_mean[3];//f
i2=4;
}

if(strcmp(RSD_fit, "shape2") == 0){
parameters2[2]=vector_mean[2];//m_BGV
parameters2[3]=vector_mean[3];//m_BGV
parameters2[4]=vector_mean[4];//f
i2=5;
}

parameters2_bao[0]=vector_mean[0];//apara
parameters2_bao[1]=vector_mean[1];//aperp

i1=2;//i2=3;
if(Nsigmas_free==0)
{

  if(strcmp(Sigma_independent, "yes") == 0)//Nsigma_tot=2
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2_bao[3]=Sigma_nl_mean[1];//sigmaperp
  parameters2[3]=Sigma_nl_mean[0];//sigmapara
  parameters2[4]=Sigma_nl_mean[1];//sigmaperp

  i1=4;//i2=3;

  }
  if(strcmp(Sigma_independent, "no") == 0)//Nsigma_tot=1
  {
  parameters2_bao[2]=Sigma_nl_mean[0];//sigmapara
  parameters2[3]=Sigma_nl_mean[0];//sigmapara

  i1=3;//i2=3;
  }

}

for(i=i1+1;i<dimensionbao+1;i++){parameters2[i]=vector_mean[i2];parameters2_bao[i-1]=vector_mean[i2];i2++;}//set the rest of parameters (check we need dimension+1)

if(strcmp(RSD_fit, "yes") == 0){
parameters2_rsd[0]=vector_mean[0];
parameters2_rsd[1]=vector_mean[1];
parameters2_rsd[2]=vector_mean[2];
for(i=3;i<dimensionrsd;i++){parameters2_rsd[i]=vector_mean[i2];parameters2[i-3+dimensionbao+1]=vector_mean[i2];i2++;}
}
if(strcmp(RSD_fit, "shape") == 0){
parameters2_rsd[0]=vector_mean[0];
parameters2_rsd[1]=vector_mean[1];
parameters2_rsd[2]=vector_mean[2];
parameters2_rsd[3]=vector_mean[3];
for(i=4;i<dimensionrsd;i++){parameters2_rsd[i]=vector_mean[i2];parameters2[i-4+dimensionbao+1]=vector_mean[i2];i2++;}
}

if(strcmp(RSD_fit, "shape2") == 0){
parameters2_rsd[0]=vector_mean[0];
parameters2_rsd[1]=vector_mean[1];
parameters2_rsd[2]=vector_mean[2];
parameters2_rsd[3]=vector_mean[3];
parameters2_rsd[4]=vector_mean[4];
for(i=5;i<dimensionrsd;i++){parameters2_rsd[i]=vector_mean[i2];parameters2[i-5+dimensionbao+1]=vector_mean[i2];i2++;}
}

}

free(vector_mean);

chi2=chi2_bao_rsd(type_BAO_fit,type_of_analysis,fit_BAO,fit_RSD,parameters2_bao,parameters2_rsd, k_Plin, Plin,N_Plin,k_Olin,Olin,N_Olin,Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0bao,k0rsd, P0bao,P0rsd,Pnoise_rsd,NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd,NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd,NeffP4bao,NeffP4rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC,Pnoise_rsdSGC,NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC,NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC,cov,covSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev, Npolynomial, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1bao, plan2bao, plan1rsd, plan2rsd, do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot,Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,bispectrum_BQ, mask_matrix, MatrixBAO_mask_NGC, MatrixBAO_mask_SGC, MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev, noise_option,covariance_correction, NrealNGC, NrealSGC);

    make_a_bao_plot(type_BAO_fit,type_of_analysis,fit_BAO,parameters2_bao,chi2, k0bao,P0bao,errP0bao,NeffP0bao,k0baoSGC,P0baoSGC,errP0baoSGC,NeffP0baoSGC, k2bao,P2bao,errP2bao,NeffP2bao,k2baoSGC,P2baoSGC,errP2baoSGC,NeffP2baoSGC,k4bao,P4bao,errP4bao,NeffP4bao,k4baoSGC,P4baoSGC,errP4baoSGC,NeffP4baoSGC, k_Plin, Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0,W2,W4,W6,W8, Nmask,path_to_mask1,spacing_maskNGC, posSGC, W0SGC,W2SGC,W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, Sigma_def_type, Sigma_independent,  ffactor,Sigma_type,  Sigma_nl_mean,  Sigma_nl_stddev,  Npolynomial, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier,plan1bao,plan2bao,do_power_spectrum, do_bispectrum, Nalphas,Nsigmas_tot, Nsigmas_free,Sigma_smooth,factor_sampling_mask,spacing_dataNGC_bao,spacing_dataSGC_bao,spacing_theory, 0,0,NULL,k11bao,k22bao,k33bao, B0bao,errB0bao,Bnoise_bao,NeffB0bao,k11baoSGC, k22baoSGC, k33baoSGC, B0baoSGC,errB0baoSGC, Bnoise_baoSGC,NeffB0baoSGC ,bispectrum_BQ,mask_matrix,MatrixBAO_mask_NGC, MatrixBAO_mask_SGC);


make_a_rsd_plot(type_of_analysis,parameters2_rsd,chi2, Theory,Ntheory, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0rsd, P0rsd, errP0rsd, Pnoise_rsd, k2rsd, P2rsd,errP2rsd, k4rsd, P4rsd,errP4rsd, k0rsdSGC, P0rsdSGC,errP0rsdSGC, Pnoise_rsdSGC, k2rsdSGC, P2rsdSGC,errP2rsdSGC, k4rsdSGC, P4rsdSGC,errP4rsdSGC, NeffP0rsd, NeffP2rsd, NeffP4rsd,NeffP0rsdSGC, NeffP2rsdSGC, NeffP4rsdSGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier, plan1rsd, plan2rsd, fit_RSD, do_power_spectrum, do_bispectrum, 0,factor_sampling_mask,spacing_dataNGC_rsd,spacing_dataSGC_rsd,spacing_theory_rsd,0, NULL, NULL, NULL,0, k11rsd,k22rsd,k33rsd,B0rsd,errB0rsd,Bnoise_rsd,NeffB0rsd,k11rsdSGC,k22rsdSGC,k33rsdSGC,B0rsdSGC,errB0rsdSGC,Bnoise_rsdSGC,NeffB0rsdSGC,bispectrum_BQ, mask_matrix, MatrixFS_mask_NGC, MatrixFS_mask_SGC,noise_option);

printf("Chi2 = %lf / (%d - %d).\n",chi2,Npoints+NpointsSGC,Ndof);
fftw_destroy_plan(plan1bao);
fftw_destroy_plan(plan2bao);
fftw_destroy_plan(plan1rsd);
fftw_destroy_plan(plan2rsd);
free(parameters2);
free(parameters2_rsd);
free(parameters2_bao);

}else{
if(strcmp(use_prop_cov, "yes") == 0)//full mcmc run with proposal
{
trial_mcmc=0;//no trial
fraction=1;
Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

read_prop_cov(NULL,0,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);
transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);

//generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc,type_of_analysis, fit_BAO, do_power_spectrum, do_bispectrum, Sigma_def_type, Sigma_independent, Sigma_type, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, Npolynomial, path_output, identifier);//writes transform and transform_inverse

free(Cov_prop);

mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,NULL,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,fit_RSD,k_Plin,Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao,k0rsd, P0bao, P0rsd, errP0bao,errP0rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, errP2bao,errP2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, errP4bao,errP4rsd, NeffP4bao,NeffP4rsd, k11bao,k11rsd, k22bao,k2rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, errP0baoSGC,errP0rsdSGC, NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, errP2baoSGC,errP2rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, errP4baoSGC,errP4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC, B0baoSGC,B0rsdSGC,errB0baoSGC,errB0rsdSGC,Bnoise_baoSGC,Bnoise_rsdSGC, NeffB0baoSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent,ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum, 0, Nalphas,Nsigmas_tot, Nsigmas_free,Theory,Ntheory, Pnoise_rsd,Pnoise_rsdSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , Sigma_smooth,factor_sampling_mask, spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS,0,0,NULL,NULL,NULL,0,bispectrum_BQ, mask_matrix, MatrixBAO_mask_NGC, MatrixBAO_mask_SGC, MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev, noise_option,step_size, covariance_correction, NrealNGC, NrealSGC);
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
//set_proposal_mean(vector_mean,N_Cov_prop);
set_proposal_mean(vector_mean, N_Cov_prop, type_of_analysis, fit_BAO, do_power_spectrum, do_bispectrum, Sigma_def_type, Sigma_independent, Sigma_type, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, Npolynomial, path_output, identifier);

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);

//generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc,type_of_analysis, fit_BAO, do_power_spectrum, do_bispectrum, Sigma_def_type, Sigma_independent, Sigma_type, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, Npolynomial, path_output, identifier);//writes transform and transform_inverse

vector_buffer = (double**) calloc(N_Cov_prop+2,sizeof(double*));

for(j=0;j<N_Cov_prop+2;j++)
{
   vector_buffer[j] = (double*) calloc(N_max,sizeof(double));
}

mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,fit_RSD,k_Plin,Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao,k0rsd, P0bao, P0rsd, errP0bao,errP0rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, errP2bao,errP2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, errP4bao,errP4rsd, NeffP4bao,NeffP4rsd, k11bao,k11rsd, k22bao,k2rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, errP0baoSGC,errP0rsdSGC, NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, errP2baoSGC,errP2rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, errP4baoSGC,errP4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC, B0baoSGC,B0rsdSGC,errB0baoSGC,errB0rsdSGC,Bnoise_baoSGC,Bnoise_rsdSGC, NeffB0baoSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent,ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum,0, Nalphas,Nsigmas_tot, Nsigmas_free,Theory,Ntheory, Pnoise_rsd,Pnoise_rsdSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , Sigma_smooth,factor_sampling_mask, spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS,0,0,NULL,NULL,NULL,0,bispectrum_BQ, mask_matrix, MatrixBAO_mask_NGC, MatrixBAO_mask_SGC, MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev, noise_option,step_size,covariance_correction, NrealNGC, NrealSGC);


Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized
//get covariance from buffer
read_prop_cov(vector_buffer,N_max,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);

freeTokens(vector_buffer,N_Cov_prop+2);

trial_mcmc=0;//no trial
fraction=1;

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
//generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc,type_of_analysis, fit_BAO, do_power_spectrum, do_bispectrum, Sigma_def_type, Sigma_independent, Sigma_type, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, Npolynomial, path_output, identifier);//writes transform and transform_inverse

free(Cov_prop);
mcmc_kernel(nthreads,type_BAO_fit,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, fit_BAO,fit_RSD,k_Plin,Plin, N_Plin, k_Olin, Olin, N_Olin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  k0bao,k0rsd, P0bao, P0rsd, errP0bao,errP0rsd, NeffP0bao,NeffP0rsd, k2bao,k2rsd, P2bao,P2rsd, errP2bao,errP2rsd, NeffP2bao,NeffP2rsd, k4bao,k4rsd, P4bao,P4rsd, errP4bao,errP4rsd, NeffP4bao,NeffP4rsd, k11bao,k11rsd, k22bao,k2rsd, k33bao,k33rsd, B0bao,B0rsd, errB0bao,errB0rsd, Bnoise_bao,Bnoise_rsd, NeffB0bao,NeffB0rsd, k0baoSGC,k0rsdSGC, P0baoSGC,P0rsdSGC, errP0baoSGC,errP0rsdSGC, NeffP0baoSGC,NeffP0rsdSGC, k2baoSGC,k2rsdSGC, P2baoSGC,P2rsdSGC, errP2baoSGC,errP2rsdSGC, NeffP2baoSGC,NeffP2rsdSGC, k4baoSGC,k4rsdSGC, P4baoSGC,P4rsdSGC, errP4baoSGC,errP4rsdSGC, NeffP4baoSGC,NeffP4rsdSGC, k11baoSGC,k11rsdSGC, k22baoSGC,k22rsdSGC, k33baoSGC,k33rsdSGC, B0baoSGC,B0rsdSGC,errB0baoSGC,errB0rsdSGC,Bnoise_baoSGC,Bnoise_rsdSGC, NeffB0baoSGC,NeffB0rsdSGC, cov, covSGC, alpha_min, alpha_max, Sigma_def_type, Sigma_independent,ffactor, Sigma_type, Sigma_nl_mean, Sigma_nl_stddev, Npolynomial, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum,0, Nalphas,Nsigmas_tot, Nsigmas_free,Theory,Ntheory, Pnoise_rsd,Pnoise_rsdSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , Sigma_smooth,factor_sampling_mask, spacing_dataNGC_bao,spacing_dataNGC_rsd,spacing_dataSGC_bao,spacing_dataSGC_rsd,spacing_theory,spacing_theory_rsd,type_of_analysis_BAO,type_of_analysis_FS,0,0,NULL,NULL,NULL,0,bispectrum_BQ, mask_matrix, MatrixBAO_mask_NGC, MatrixBAO_mask_SGC, MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev,noise_option,step_size,covariance_correction, NrealNGC, NrealSGC);

}

}

}


void do_rsd_mcmc(int nthreads, char *type_of_analysis, char *fit_RSD,double **Theory,int N_Plin,double *k_Plin_nw,double *Plin_nw,int N_Plin_nw, double *k_Olin, double *Olin, int N_Olin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *path_to_mask1, char *spacing_maskNGC, double *posSGC, double *W0SGC, double *W2SGC, double *W4SGC, double *W6SGC, double *W8SGC, int NmaskSGC,  char *path_to_mask2, char *spacing_maskSGC, double *k0, double *P0, double *errP0, double Pnoise, int NeffP0, double *k2, double *P2, double *errP2, int NeffP2, double *k4, double *P4, double *errP4, int NeffP4, double *k11, double *k22, double *k33, double *B0, double *errB0, double *Bnoise, int NeffB0, double *k0SGC, double *P0SGC, double *errP0SGC,double PnoiseSGC,int NeffP0SGC, double *k2SGC, double *P2SGC, double *errP2SGC,int NeffP2SGC, double *k4SGC, double *P4SGC, double *errP4SGC,int NeffP4SGC, double *k11SGC, double *k22SGC, double *k33SGC, double *B0SGC,double *errB0SGC, double *BnoiseSGC,int NeffB0SGC, double *cov, double *covSGC, double alpha_min, double alpha_max, char *ptmodel_ps, char *rsdmodel_ps, char *fogmodel_ps, char *ptmodel_bs, char *local_b2s2, char *local_b3nl,char *RSD_fit, char *sigma8_free, char *fog_free,char *fog_bs ,int Nchunks, char *path_output, char *identifier, char *do_plot, char *use_prop_cov, char *path_to_cov, long int Nsteps, char *do_power_spectrum, char *do_bispectrum, double redshift_in, char *spacing_dataNGC,char *spacing_dataSGC, char *spacing_theory,char *type_of_analysis_BAO,char *type_of_analysis_FS,char *bispectrum_BQ, int factor_sampling_mask, char *mask_matrix,double **MatrixFS_mask_NGC, double **MatrixFS_mask_SGC,double *FSprior_type, double *FSprior_mean, double *FSprior_stddev,int noise_option, double step_size,char *covariance_correction, int NrealNGC, int NrealSGC, double *P0bao, double *P0baoSGC)
{
double fraction;
int trial_mcmc;
long int N_max,N_print,N_burnout,j;
double **vector_buffer;
FILE *f;
double knl;
double *knl_y,*sigma8_x;
int Nknl;
double *Plin,*klin;
double *n_func,*k_n_func,*n_func_final;
double *Cov_prop;
double *vector_mean;
int N_Cov_prop;
gsl_matrix *transform_inverse;
gsl_matrix *transform;
char name_file[2000];
long int *params_mcmc;

int i_thread;
int i;

//double *parameters2;
//int dimension,dimensionP;
int Nplan;
double chi2;
double params[5];
int Ndof, Npoints, NpointsSGC;
fftw_complex *a_pointer;a_pointer=NULL;
fftw_complex *b_pointer;b_pointer=NULL;
fftw_plan plan1,plan2;
int i1,i2;
int modeP0,modeP2,modeP4,modeB0;

if (nthreads<2){
  sprintf(name_file,"%s/mcmcFS_output_%s.txt",path_output,identifier);
  f=fopen(name_file,"w");
  if(f==NULL){printf("Error writing %s. Exiting now...\n",name_file);exit(0);}
  fclose(f);
}else{
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
if(strcmp(fog_free, "yes") == 0){N_Cov_prop++;}
if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0){N_Cov_prop++;}
if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 ){N_Cov_prop++;}
N_Cov_prop=N_Cov_prop*Nchunks;//x2 
if(strcmp(sigma8_free, "yes") == 0){N_Cov_prop++;}
if(strcmp(RSD_fit, "yes") == 0){N_Cov_prop=N_Cov_prop+3;}//f, alpha-para, alpha-perp
if(strcmp(RSD_fit, "shape") == 0){N_Cov_prop=N_Cov_prop+4;}//f, alpha-para, alpha-perp, m_BGV
if(strcmp(RSD_fit, "shape2") == 0){N_Cov_prop=N_Cov_prop+5;}//f, alpha-para, alpha-perp, m1_BGV, m2_BGV

if(strcmp(do_bispectrum, "yes") == 0)
{
//get knl fid
//get n(k)
//get knl(s8) if sigma8_free
knl=0;
n_func_final=NULL;
knl_y=NULL;
sigma8_x=NULL;
if(strcmp(ptmodel_bs,"GilMarin14") == 0){

Nknl=1000;

Plin = malloc(sizeof(double)*N_Plin);
klin = malloc(sizeof(double)*N_Plin);
for(i=0;i<N_Plin;i++)
{
Plin[i]=Theory[i][1];
klin[i]=Theory[i][0];
}
knl=get_fid_knl(klin,Plin,N_Plin);

if(knl<=0)
{
if(knl==-1){printf("knl could not be computed from theory file (knl>kmax). Computing it from Olin times Psmooth....\n");}
if(knl==-2){printf("knl could not be computed from theory file (knl<kmin). Computing it from Olin times Psmooth....\n");}
if(N_Olin<=0){printf("Olin file not provided. Exiting now...\n");exit(0);}
else{
free(Plin);
free(klin);
Plin = malloc(sizeof(double)*N_Olin);
klin = malloc(sizeof(double)*N_Olin);
for(i=0;i<N_Olin;i++)
{
Plin[i]=Plin_nw[i]*Olin[i];
klin[i]=k_Olin[i];
}
knl=get_fid_knl(klin,Plin,N_Olin);

if(knl==-2){printf("Error determining knl, seems that knl is smaller than the kmin provided. Exiting now...\n");exit(0);}
if(knl==-1){printf("Warning, when determining knl the kmax in Plin is not large enough. Extrapolating knl...\n");

knl = (2*Pi*Pi*(klin[N_Olin-1]-klin[N_Olin-2])-pow(klin[N_Olin-1],3)*Plin[N_Olin-1]*(klin[N_Olin-1]-klin[N_Olin-2])+klin[N_Olin-1]*(pow(klin[N_Olin-1],3)*Plin[N_Olin-1]-pow(klin[N_Olin-2],3)*Plin[N_Olin-2]))/(pow(klin[N_Olin-1],3)*Plin[N_Olin-1]-pow(klin[N_Olin-2],3)*Plin[N_Olin-2]);
}
printf("Fiducial knl computed for bispectrum: knl=%lf\n",knl);

knl_y = malloc(sizeof(double)*Nknl);
sigma8_x = malloc(sizeof(double)*Nknl);
generate_knl_array(klin, Plin, N_Olin, knl_y, sigma8_x, Nknl, Theory[0][41]);//check this

//this is to has the same k-spacing as if knl comes from Theory when computing n_func_final later
free(klin);
klin = malloc(sizeof(double)*N_Plin);
for(i=0;i<N_Plin;i++){klin[i]=Theory[i][0];}

}

}
else{
printf("Fiducial knl computed for bispectrum: knl=%lf\n",knl);

knl_y = malloc(sizeof(double)*Nknl);
sigma8_x = malloc(sizeof(double)*Nknl);
generate_knl_array(klin, Plin, N_Plin, knl_y, sigma8_x, Nknl, Theory[0][41]);//check this
}

//for(i=0;i<Nknl;i++){printf("%lf %lf\n",sigma8_x[i],knl_y[i]);}

n_func=malloc(sizeof(double)*(N_Plin_nw-1));
k_n_func=malloc(sizeof(double)*(N_Plin_nw-1));
generate_n_array(k_Plin_nw, Plin_nw, N_Plin_nw, n_func, k_n_func);
//for(i=0;i<N_Plin_nw-1;i++){printf("%lf %lf\n",k_n_func[i],n_func[i]);}
n_func_final=malloc(sizeof(double)*N_Plin);
smooth_n_array(klin, N_Plin, n_func_final, n_func, k_n_func, N_Plin_nw-1);//this should be fine but check it

//for(i=0;i<N_Plin;i++){printf("%lf %lf\n",klin[i],n_func_final[i]);}

free(n_func);
free(k_n_func);
free(Plin);
free(klin);

}//GilMarin14

}//do_bispectrum
//exit(0);//remove this after checking things are ok

if(Nsteps<=0)//in this case we just make  a plot
{
printf("Warning, invalid number of steps selected for the mcmc (%ld). Only ploting function at prior set of parameters.\n",Nsteps);
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

if(strcmp(use_prop_cov, "yes") == 0)//full mcmc run with proposal
{
trial_mcmc=0;//no trial
Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
read_prop_cov(NULL,0,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);
free(Cov_prop);
}
if(strcmp(use_prop_cov, "no") == 0)//from set_propsal_mean
{
set_proposal_mean(vector_mean, N_Cov_prop, type_of_analysis, fit_RSD, do_power_spectrum, do_bispectrum, NULL, NULL, NULL, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, 0, path_output,identifier);
}


set_mask_params(params,0.0,1.0,0.0,1.0,0.0);

Nplan=(int)(params[2]);

    /*fftw_plan*/ plan1 = fftw_plan_dft_1d(Nplan,  a_pointer,  b_pointer,  -1, FFTW_ESTIMATE);//forward plan
    /*fftw_plan*/ plan2 = fftw_plan_dft_1d(Nplan,  b_pointer,  b_pointer, +1, FFTW_ESTIMATE);//reverse plan

Ndof=N_Cov_prop;

//define modeP0, modeP2,modeP4,modeB0
modeP0=0;modeP2=0;modeP4=0;modeB0=0;
if( strcmp(do_power_spectrum,"yes") ==0){

if( strcmp(fit_RSD,"P0") == 0 || strcmp(fit_RSD,"P02") == 0 || strcmp(fit_RSD,"P04") == 0 || strcmp(fit_RSD,"P024") == 0)
{
modeP0=1;
}
if( strcmp(fit_RSD,"P2") == 0 || strcmp(fit_RSD,"P02") == 0 || strcmp(fit_RSD,"P24") == 0 || strcmp(fit_RSD,"P024") == 0)
{
modeP2=1;
}

if( strcmp(fit_RSD,"P4") == 0 || strcmp(fit_RSD,"P04") == 0 || strcmp(fit_RSD,"P24") == 0 || strcmp(fit_RSD,"P024") == 0)
{
modeP4=1;
}

}
if( strcmp(do_bispectrum,"yes") ==0){
modeB0=1;
}

Npoints=NeffP0*modeP0+NeffP2*modeP2+NeffP4*modeP4+NeffB0*modeB0;
if(P0bao!=NULL){Npoints=Npoints+2;}
NpointsSGC=0;
if(Nchunks==2){NpointsSGC=NeffP0SGC*modeP0+NeffP2SGC*modeP2+NeffP4SGC*modeP4+NeffB0SGC*modeB0;if(P0baoSGC!=NULL){NpointsSGC=NpointsSGC+2;}}

//chi2=chi2_rsd();
chi2=chi2_rsd_mcmc(type_of_analysis,vector_mean, Theory,N_Plin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0, P0, Pnoise, k2, P2, k4, P4, k0SGC, P0SGC, PnoiseSGC, k2SGC, P2SGC, k4SGC, P4SGC, NeffP0, NeffP2, NeffP4,NeffP0SGC,NeffP2SGC,NeffP4SGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks, plan1, plan2, fit_RSD, do_power_spectrum, do_bispectrum,redshift_in,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory,knl, n_func_final, sigma8_x, knl_y,Nknl, k11,k22,k33,B0,Bnoise,NeffB0,k11SGC,k22SGC,k33SGC,B0SGC,BnoiseSGC,NeffB0SGC,bispectrum_BQ, mask_matrix, MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev, noise_option,covariance_correction, NrealNGC, NrealSGC,P0bao,P0baoSGC);

//makd_a_rsd_plot();
make_a_rsd_plot(type_of_analysis,vector_mean,chi2, Theory,N_Plin, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC,W6SGC,W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC, k0, P0, errP0, Pnoise, k2, P2,errP2, k4, P4,errP4, k0SGC, P0SGC,errP0SGC, PnoiseSGC, k2SGC, P2SGC,errP2SGC, k4SGC, P4SGC,errP4SGC, NeffP0, NeffP2, NeffP4,NeffP0SGC, NeffP2SGC, NeffP4SGC, cov, covSGC,  ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs, local_b2s2, local_b3nl,RSD_fit, sigma8_free, fog_free, fog_bs, Nchunks,Npoints+NpointsSGC,Ndof, path_output, identifier, plan1, plan2, fit_RSD, do_power_spectrum, do_bispectrum,redshift_in,factor_sampling_mask,spacing_dataNGC,spacing_dataSGC,spacing_theory,knl, n_func_final, sigma8_x, knl_y,Nknl, k11,k22,k33,B0,errB0,Bnoise,NeffB0,k11SGC,k22SGC,k33SGC,B0SGC,errB0SGC,BnoiseSGC,NeffB0SGC,bispectrum_BQ, mask_matrix, MatrixFS_mask_NGC, MatrixFS_mask_SGC,noise_option);

free(vector_mean);
fftw_destroy_plan(plan1);
fftw_destroy_plan(plan2);

printf("Chi2 = %lf / (%d - %d).\n",chi2,Npoints+NpointsSGC,Ndof);

}else{

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

//generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc,type_of_analysis, fit_RSD, do_power_spectrum, do_bispectrum, NULL, NULL, NULL, local_b2s2, local_b3nl, fog_free,fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, 0, path_output,identifier);//writes transform and transform_inverse

free(Cov_prop);

mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,NULL,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, NULL,fit_RSD,NULL,NULL, 0, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  NULL,k0, P0bao ,P0, NULL,errP0, 0,NeffP0, NULL,k2, NULL,P2, NULL,errP2, 0,NeffP2, NULL,k4, NULL,P4, NULL,errP4, 0,NeffP4, NULL,k11, NULL,k22, NULL,k33, NULL,B0, NULL,errB0, NULL,Bnoise, 0,NeffB0, NULL,k0SGC, P0baoSGC,P0SGC, NULL,errP0SGC, 0,NeffP0SGC, NULL,k2SGC, NULL,P2SGC, NULL,errP2SGC, 0,NeffP2SGC, NULL,k4SGC, NULL,P4SGC, NULL,errP4SGC, 0,NeffP4SGC, NULL,k11SGC, NULL,k22SGC, NULL,k33SGC, NULL,B0SGC, NULL,errB0SGC,NULL,BnoiseSGC, 0,NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum,redshift_in, 0,0, 0,Theory,N_Plin,Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , 0,factor_sampling_mask, NULL,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS,knl,0,n_func_final,sigma8_x,knl_y,Nknl,bispectrum_BQ, mask_matrix, NULL, NULL, MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev, noise_option,step_size,covariance_correction, NrealNGC, NrealSGC);



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
//set_proposal_mean(vector_mean,N_Cov_prop);
set_proposal_mean(vector_mean, N_Cov_prop, type_of_analysis, fit_RSD, do_power_spectrum, do_bispectrum, NULL, NULL, NULL, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, 0, path_output,identifier);

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);

//generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc, type_of_analysis, fit_RSD, do_power_spectrum, do_bispectrum, NULL, NULL, NULL, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, 0, path_output,identifier);//writes transform and transform_inverse

vector_buffer = (double**) calloc(N_Cov_prop+2,sizeof(double*));

for(j=0;j<N_Cov_prop+2;j++)
{
   vector_buffer[j] = (double*) calloc(N_max,sizeof(double));
}

mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, NULL,fit_RSD,NULL,NULL, 0, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  NULL,k0, P0bao,P0, NULL,errP0, 0,NeffP0, NULL,k2, NULL,P2, NULL,errP2, 0,NeffP2, NULL,k4, NULL,P4, NULL,errP4, 0,NeffP4, NULL,k11, NULL,k22, NULL,k33, NULL,B0, NULL,errB0, NULL,Bnoise, 0,NeffB0, NULL,k0SGC, P0baoSGC,P0SGC, NULL,errP0SGC, 0,NeffP0SGC, NULL,k2SGC, NULL,P2SGC, NULL,errP2SGC, 0,NeffP2SGC, NULL,k4SGC, NULL,P4SGC, NULL,errP4SGC, 0,NeffP4SGC, NULL,k11SGC, NULL,k22SGC, NULL,k33SGC, NULL,B0SGC, NULL,errB0SGC,NULL,BnoiseSGC, 0,NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum,redshift_in, 0,0, 0,Theory,N_Plin,Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , 0,factor_sampling_mask, NULL,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS,knl,0,n_func_final,sigma8_x,knl_y,Nknl,bispectrum_BQ, mask_matrix, NULL, NULL,MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev,noise_option,step_size,covariance_correction, NrealNGC, NrealSGC);

Cov_prop=(double*) calloc(N_Cov_prop*N_Cov_prop, sizeof(double));
vector_mean= (double*) calloc( N_Cov_prop, sizeof(double));//zero-inizialized

//get covariance from buffer
read_prop_cov(vector_buffer,N_max,trial_mcmc,path_to_cov,Cov_prop,vector_mean,N_Cov_prop);

freeTokens(vector_buffer,N_Cov_prop+2);

trial_mcmc=0;//no trial
fraction=1;

transform_inverse = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
transform = gsl_matrix_alloc (N_Cov_prop, N_Cov_prop);
//generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc);//writes transform and transform_inverse
generate_rotation_matrix(N_Cov_prop,Cov_prop,vector_mean,transform_inverse,transform,trial_mcmc,type_of_analysis, fit_RSD, do_power_spectrum, do_bispectrum, NULL, NULL, NULL, local_b2s2, local_b3nl, fog_free, fogmodel_ps, fog_bs, sigma8_free, RSD_fit, Nchunks, 0, path_output,identifier);//writes transform and transform_inverse

free(Cov_prop);

mcmc_kernel(nthreads,NULL,type_of_analysis,trial_mcmc,vector_buffer,fraction, N_Cov_prop, transform, transform_inverse, vector_mean, name_file, NULL,fit_RSD,NULL,NULL, 0, NULL, NULL, 0, pos, W0, W2, W4, W6, W8, Nmask, path_to_mask1,spacing_maskNGC, posSGC, W0SGC, W2SGC, W4SGC, W6SGC, W8SGC, NmaskSGC, path_to_mask2,spacing_maskSGC,  NULL,k0, P0bao,P0, NULL,errP0, 0,NeffP0, NULL,k2, NULL,P2, NULL,errP2, 0,NeffP2, NULL,k4, NULL,P4, NULL,errP4, 0,NeffP4, NULL,k11, NULL,k22, NULL,k33, NULL,B0, NULL,errB0, NULL,Bnoise, 0,NeffB0, NULL,k0SGC, P0baoSGC,P0SGC, NULL,errP0SGC, 0,NeffP0SGC, NULL,k2SGC, NULL,P2SGC, NULL,errP2SGC, 0,NeffP2SGC, NULL,k4SGC, NULL,P4SGC, NULL,errP4SGC, 0,NeffP4SGC, NULL,k11SGC, NULL,k22SGC, NULL,k33SGC, NULL,B0SGC, NULL,errB0SGC,NULL,BnoiseSGC, 0,NeffB0SGC, cov, covSGC, alpha_min, alpha_max, NULL, NULL, 0, NULL, NULL, NULL, 0, Nchunks, path_output, identifier, do_plot, Nsteps, do_power_spectrum, do_bispectrum,redshift_in, 0,0, 0,Theory,N_Plin,Pnoise,PnoiseSGC, ptmodel_ps, rsdmodel_ps, fogmodel_ps, ptmodel_bs,  local_b2s2,  local_b3nl, RSD_fit, sigma8_free, fog_free, fog_bs , 0,factor_sampling_mask, NULL,spacing_dataNGC,NULL,spacing_dataSGC,NULL,spacing_theory,type_of_analysis_BAO,type_of_analysis_FS,knl,0,n_func_final,sigma8_x,knl_y,Nknl,bispectrum_BQ,mask_matrix,NULL,NULL,MatrixFS_mask_NGC, MatrixFS_mask_SGC, FSprior_type, FSprior_mean, FSprior_stddev, noise_option,step_size,covariance_correction, NrealNGC, NrealSGC);

free(vector_mean);
}

}//Nsteps>0

if( strcmp(do_bispectrum,"yes") == 0){
free(knl_y);
free(sigma8_x);
free(n_func_final);
}

}


