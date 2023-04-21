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

double F2(double k1, double k2, double k3)
{
double f;
double ctheta;
ctheta=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
f=5./7.+1./2.*ctheta*(k1/k2+k2/k1)+2./7.*ctheta*ctheta;
return f;
}
double G2(double k1, double k2, double k3)
{
double f;
double ctheta;
ctheta=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
f=3./7.+1./2.*ctheta*(k1/k2+k2/k1)+4./7.*ctheta*ctheta;
return f;
}

double Q3(double n)
{
double f;
f=(4.-pow(2,n))/(1.+pow(2.,n+1.) );
return f;
}

double aG(double n, double k, double knl,double s8)
{

double A;
double sigma8=s8;//this is the value at that redshift
double a6=7.431407;
double a2=-3.878979;
double a1=3.599090;

double q=k/knl;

A=(1+pow(sigma8,a6)*pow(0.7*Q3(n),0.5)*pow(q*a1,n+a2))/(1.+pow(q*a1 ,n+a2) );
return A;
}

double bG(double n, double k, double knl){
double B;
double a3=0.518141;
double a8=-3.103696;
double a7=5.022235;
double q=k/knl;

B=(1.+0.2*a3*(n+3.)*pow(q*a7  ,n+3+a8)    )/( 1.+pow(q*a7,n+3.5+a8)   );

return B;
}

double cG(double n, double k, double knl)
{
double C;
double a4=-3.587849;
double a9=-0.483688;
double a5=0.336059;
double q=k/knl;

C=(1.+4.5*a4/(1.5+pow(n+3,4)  )*pow(q*a5,n+3.+a9))/( 1.+pow( q*a5 ,n+3.5+a9)   );

return C;
}

double S_ker(double k1, double k2, double k3)
{
double f,ctheta;
ctheta=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
f=ctheta*ctheta-1./3.;
return f;
}


double aF(double n, double k, double knl,double s8)
{
double A;
double sigma8=s8;
double a6=-0.575;
double a2=3.740;
double a1=0.484;
double q=k/knl;


A=(1+pow(sigma8,a6)*pow(0.7*Q3(n),0.5)*pow(q*a1,n+a2))/(1.+pow(q*a1 ,n+a2) );
return A;
}

double bF(double n, double k, double knl)
{
double B;
double a3=-0.849;
double a8=-0.722;
double a7=0.128;
double q=k/knl;

B=(1.+0.2*a3*(n+3.)*pow(q*a7  ,n+3+a8)    )/( 1.+pow(q*a7,n+3.5+a8)   );

return B;
}

double cF(double n, double k, double knl)
{
double C;
double a4=0.392;
double a9=-0.926;
double a5=1.013;
double q=k/knl;

C=(1.+4.5*a4/(1.5+pow(n+3,4)  )*pow(q*a5,n+3.+a9))/( 1.+pow( q*a5 ,n+3.5+a9)   );

return C;
}

double F2_eff(double k1, double k2, double k3, double knl, double **Theory, double *k_n, double *n,int N,double s8)
{
double f;
double ctheta,n1,n2,n3;
ctheta=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
if(Theory == NULL){
n1=P_interpol(k1,k_n,n,N);
n2=P_interpol(k2,k_n,n,N);
n3=P_interpol(k3,k_n,n,N);
}
else
{
n1=P_interpol_doublearray2(k1,Theory,n,N);
n2=P_interpol_doublearray2(k2,Theory,n,N);
n3=P_interpol_doublearray2(k3,Theory,n,N);
}
f=5./7.*aF(n1,k1,knl,s8)*aF(n2,k2,knl,s8)+1./2.*ctheta*(k1/k2+k2/k1)*bF(n1,k1,knl)*bF(n2,k2,knl)+2./7.*ctheta*ctheta*cF(n1,k1,knl)*cF(n2,k2,knl);
//printf("%lf %lf %lf, %lf %lf %lf, %lf %lf %lf,%lf,%lf %lf %lf\n",k1,k2,k3,n1,n2,n3, aF(n1,k1,knl,s8)*aF(n2,k2,knl,s8), bF(n1,k1,knl)*bG(n2,k2,knl), cF(n1,k1,knl)*cF(n2,k2,knl),f, 5./7.*aF(n1,k1,knl,s8)*aF(n2,k2,knl,s8), 1./2.*ctheta*(k1/k2+k2/k1)*bF(n1,k1,knl)*bF(n2,k2,knl), 2./7.*ctheta*ctheta*cF(n1,k1,knl)*cF(n2,k2,knl));exit(0);
return f;
}

double G2_eff(double k1, double k2, double k3, double knl, double **Theory, double *k_n, double *n, int N, double s8)
{
double f;
double ctheta,n1,n2,n3;
ctheta=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);
if(Theory == NULL){
n1=P_interpol(k1,k_n,n,N);
n2=P_interpol(k2,k_n,n,N);
n3=P_interpol(k3,k_n,n,N);
}
else
{
n1=P_interpol_doublearray2(k1,Theory,n,N);
n2=P_interpol_doublearray2(k2,Theory,n,N);
n3=P_interpol_doublearray2(k3,Theory,n,N);
}
f=3./7.*aG(n1,k1,knl,s8)*aG(n2,k2,knl,s8)+1./2.*ctheta*(k1/k2+k2/k1)*bG(n1,k1,knl)*bG(n2,k2,knl)+4./7.*ctheta*ctheta*cG(n1,k1,knl)*cG(n2,k2,knl);
//printf("%lf %lf %lf, %lf %lf %lf, %lf %lf %lf,%lf,%lf %lf %lf\n",k1,k2,k3,n1,n2,n3, aG(n1,k1,knl,s8)*aG(n2,k2,knl,s8), bG(n1,k1,knl)*bG(n2,k2,knl), cG(n1,k1,knl)*cG(n2,k2,knl),f, 3./7.*aG(n1,k1,knl,s8)*aG(n2,k2,knl,s8), 1./2.*ctheta*(k1/k2+k2/k1)*bG(n1,k1,knl)*bG(n2,k2,knl), 4./7.*ctheta*ctheta*cG(n1,k1,knl)*cG(n2,k2,knl));exit(0);
return f;
}

double M1(double beta, double k1, double k2, double k3)
{
double m;
double x;
x=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);

m=2./15.*(15.+10*beta+beta*beta*(1.+2.*x*x));

return m;
}


double M2(double beta, double k1, double k2, double k3)
{
double m,x,k12,k12x,k2sq,k1sq,k3sq;
x=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);

k12=k1*k2;
k12x=k1*k2*x;
k2sq=k2*k2;
k1sq=k1*k1;
k3sq=k3*k3;

m=2*beta/(105.*k3sq)*(35.*k1sq+28.*beta*k1sq+3*beta*beta*k1sq+35.*k2sq+28.*beta*k2sq+3.*beta*beta*k2sq+70.*k12x+84.*beta*k12x+18.*beta*beta*k12x+14.*beta*k1sq*x*x+12.*beta*beta*k1sq*x*x+14.*beta+k2sq*x*x+12.*beta*beta*k2sq*x*x+12.*beta*beta*k12x*x*x);
return m;
}

double M3(double beta, double k1, double k2, double k3)
{
double m,x,k12,k12x,k2sq,k1sq,k3sq;
x=(k3*k3-k1*k1-k2*k2)/(2.*k1*k2);

k12=k1*k2;
k12x=k1*k2*x;
k2sq=k2*k2;
k1sq=k1*k1;
k3sq=k3*k3;

m=beta/(315.*k12)*(210.*k12+210.*beta*k12+54.*beta*beta*k12+6.*beta*beta*beta*k12+105.*k1sq*x+189.*beta*k1sq*x+99.*beta*beta*k1sq*x+15.*beta*beta*beta*k1sq*x+105.*k2sq*x+189.*beta*k2sq*x+99*beta*beta*k2sq*x+15.*beta*beta*beta*k2sq*x+168.*beta*k12x*x+216.*beta*beta*k12x*x+48.*beta*beta*beta*k12x*x+36.*beta*beta*k1sq*x*x*x+20.*beta*beta*beta*k1sq*x*x*x+36.*beta*beta*k2sq*x*x*x+20.*beta*beta*beta*k2sq*x*x*x+16.*beta*beta*beta*k12x*x*x*x);
return m;
}

double Zeff(double k1, double k2, double k3, double knl, double *k_n, double *n,int N,double s8,double beta_F,double beta_G, double beta_S, double beta_mu,double C1,double C2)
{
double Z;
          Z=M1(beta_F,k1,k2,k3)*F2_eff(k1,k2,k3,knl,NULL,k_n,n,N,s8)+M2(beta_G,k1,k2,k3)*G2_eff(k1,k2,k3,knl,NULL,k_n,n,N,s8)+0.5*M1(beta_S,k1,k2,k3)*(C1*S_ker(k1,k2,k3)+C2)+M3(beta_mu,k1,k2,k3);

//printf("%lf %lf,%lf\n",F2_eff(k1,k2,k3,knl,NULL,k_n,n,N,s8),G2_eff(k1,k2,k3,knl,NULL,k_n,n,N,s8),Z);exit(0);

return Z;
}

double Zspt(double k1, double k2, double k3, double knl, double *k_n, double *n,int N,double s8,double beta_F,double beta_G, double beta_S, double beta_mu,double C1,double C2)
{
double Z;
          Z=M1(beta_F,k1,k2,k3)*F2(k1,k2,k3)+M2(beta_G,k1,k2,k3)*G2(k1,k2,k3)+0.5*M1(beta_S,k1,k2,k3)*(C1*S_ker(k1,k2,k3)+C2)+M3(beta_mu,k1,k2,k3);

return Z;
}

void do_Btheo_RSD(char *ptmodel, char *type_fog, char *type_bispectrummodel,double B_theo0[], int NeffB0,double *NoiseB, double *parameters1,double **Theory, int Nlin, double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *spacing_mask, char *path_to_mask1, double k_1[], double k_2[], double k_3[], fftw_plan plan1, fftw_plan plan2, double kmin_data0 , double kmax_data0, char *spacing_theory,double knl,double *n, char *bispectrum_BQ, char *mask_matrix, double **MatrixFS_mask, int noise_option, double redshift_in)//no mask-matrix yet in the bispectrum
{

 f_params *function_parameters;
 f_params *function_parameters_real;

  function_parameters = (f_params *) malloc(sizeof(f_params));
  function_parameters_real = (f_params *) malloc(sizeof(f_params));
char ptmodel_linear[200];sprintf(ptmodel_linear,"linear");
double B12,B13,B23,junk;
double XMIN2[2]={-1,0};
double XMAX2[2]={+1,2.*Pi};
double XMIN1[1]={-1};
double XMAX1[1]={+1};
double precision=1e-1;
int i;
double *Preal_nomask,*kreal_nomask;
int N_nomask;
double P1,P2,P3;
char type_of_analysis[2000];

double b1,b2,Anoise,sigma8_scaling,sigma8,f,mBGV,m2BGV;
double b3nl,sigmaP,bs2,sigmaB,apara,aperp,avir;

int interpolation_order,shiftN;
double w1_k1,w2_k1,w0_k1;
double w1_k2,w2_k2,w0_k2;
double w1_k3,w2_k3,w0_k3;
int Ninterpol1,Ninterpol2,Ninterpol3;
char spacing_data[2000];


sprintf(spacing_data,"linear");//always true for bispectrum!
interpolation_order=1;
if(interpolation_order==1){shiftN=1;}
if(interpolation_order==2){shiftN=2;}

apara=parameters1[0];
aperp=parameters1[1];
f=parameters1[4];
mBGV=parameters1[2];
m2BGV=parameters1[3];
sigma8=parameters1[5];
sigma8_scaling=pow(parameters1[5]/Theory[0][41],2);
b1=parameters1[6];
b2=parameters1[7];
Anoise=parameters1[8];
bs2=parameters1[9];//-4./7.*(b1-1)
b3nl=parameters1[10];//needed for Q
sigmaP=parameters1[11];//needed for Q
sigmaB=parameters1[12];//for B this is placed always here
avir=parameters1[13];


//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", apara,aperp,f,sigma8,sigma8_scaling,b1,b2,Anoise,bs2,sigmaB);
(*function_parameters).redshift_in=redshift_in;
(*function_parameters).type_ptmodel=ptmodel;

if( strcmp(type_bispectrummodel,"tree-level") == 0){(*function_parameters_real).type_ptmodel=ptmodel_linear;}
else{(*function_parameters_real).type_ptmodel=ptmodel;} //if tree-level, input real

(*function_parameters).type_fog=type_fog;
(*function_parameters).type_bispectrummodel=type_bispectrummodel;
(*function_parameters).knl=knl;//either fiducial (if s8 fixed), or the interpolated one
(*function_parameters).n=n;

(*function_parameters).b1=b1;
(*function_parameters).b2=b2;
(*function_parameters).A=Anoise;
(*function_parameters).bs2=bs2;
(*function_parameters).b3nl=b3nl;
(*function_parameters).sigma8_value=sigma8;
(*function_parameters).sigma8=sigma8_scaling;
(*function_parameters_real).sigma8=sigma8_scaling;

(*function_parameters).sigmaB=sigmaB;
(*function_parameters).sigmaP=sigmaP;
(*function_parameters).avir=avir;
(*function_parameters).f=f;
(*function_parameters).m2BGV=m2BGV;
(*function_parameters_real).m2BGV=m2BGV;
(*function_parameters).mBGV=mBGV;
(*function_parameters_real).mBGV=mBGV;
(*function_parameters).a_parallel=apara;
(*function_parameters).a_perpendicular=aperp;
(*function_parameters).theory=Theory;
(*function_parameters_real).theory=Theory;
(*function_parameters).N=Nlin;
(*function_parameters_real).N=Nlin;
(*function_parameters).spacing=spacing_theory;
(*function_parameters_real).spacing=spacing_theory;
(*function_parameters).noise_option=noise_option;


if( strcmp(path_to_mask1,"none") ==0){N_nomask=1;}
else{N_nomask=Nlin;}

                kreal_nomask = (double*) calloc(N_nomask, sizeof(double));
                Preal_nomask = (double*) calloc(N_nomask, sizeof(double));

   if( strcmp(path_to_mask1,"none") !=0)//mask app
   {
      for(i=0;i<N_nomask;i++)
      {
         kreal_nomask[i]=Theory[i][0];
      }

//apply first window, then distort AP (approximation, but much faster)
(*function_parameters_real).a_parallel=1;
(*function_parameters_real).a_perpendicular=1;

      Preal(kreal_nomask,Preal_nomask,N_nomask, function_parameters_real);

     
      //apply mask here
     sprintf(type_of_analysis,"BAOISO");//we apply the mask as if it were an isotropic bao, i.e. only through W0, ignoring W2,...
     apply_mask(type_of_analysis,1, 0, 0, kreal_nomask,Preal_nomask, NULL, NULL,N_nomask, pos, W0, W2, W4, W6, W8, Nmask,spacing_mask, plan1, plan2, kreal_nomask[0], kreal_nomask[N_nomask-1], kmin_data0, kmax_data0, kmin_data0, kmax_data0, kmin_data0, kmax_data0,spacing_theory,0,1);

//pass the real P(k) without apara, aperp distortions (but with m(apara=1, aperp=1) distortion)
(*function_parameters).Nreal=N_nomask;
(*function_parameters).kreal=kreal_nomask;
(*function_parameters).Preal=Preal_nomask;

}//if mask=none do nothing, just tell integral_B that mask is none

(*function_parameters).mask=path_to_mask1;

   for(i=0;i<NeffB0;i++)
   {
          (*function_parameters).k1input=k_1[i];
          (*function_parameters).k2input=k_2[i];
          (*function_parameters).k3input=k_3[i];
          adapt_integrate(1,integralB,function_parameters,2,XMIN2,XMAX2,0,precision,precision,&B12,&junk);

          (*function_parameters).k1input=k_1[i];
          (*function_parameters).k2input=k_3[i];
          (*function_parameters).k3input=k_2[i];
          adapt_integrate(1,integralB,function_parameters,2,XMIN2,XMAX2,0,precision,precision,&B13,&junk);

          (*function_parameters).k1input=k_2[i];
          (*function_parameters).k2input=k_3[i];
          (*function_parameters).k3input=k_1[i];
          adapt_integrate(1,integralB,function_parameters,2,XMIN2,XMAX2,0,precision,precision,&B23,&junk);

//if( strcmp(bispectrum_BQ,"B") == 0){
B_theo0[i]=(2.*(B12+B13+B23)+Anoise*NoiseB[i])*pow(apara,-2)*pow(aperp,-4);

if( strcmp(bispectrum_BQ,"Q") == 0 && N_nomask==1){//reduced bispectrum and without mask

      (*function_parameters).mode=0;//compute the monopole

         (*function_parameters).kinput=k_1[i];
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&P1,&junk);

         (*function_parameters).kinput=k_2[i];
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&P2,&junk);

         (*function_parameters).kinput=k_3[i];
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&P3,&junk);

         B_theo0[i]=B_theo0[i]/(P1*P2+P1*P3+P2*P3);

}



//}

/*
if( strcmp(bispectrum_BQ,"Q") == 0){
//B_theo0[i]=2.*(P1*P2*B12+P1*P3*B13+P2*P3*B23)/(P1*P2+P1*P3+P2*P3)+Anoise*NoiseB[i];//note that the volume element (apara*aperp^2)^2 cancells out in Q!

// k, mu integration of Preal 
//
// //mode == 0 (monopole)
P1=0;
            adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&P1,&junk);
P2=0;
            adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&P1,&junk);
P3=0;
            adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&P1,&junk);

//apply window to P's as iso-only for simplicity



B_theo0[i]=2.*(B12+B13+B23)/(P1*P2+P1*P3+P2*P3)+Anoise*NoiseB[i];//note that the volume element (apara*aperp^2)^2 cancells out in Q!

}*/

//B_theo0[i]=P2+Anoise*NoiseB[i];

///B_theo0[i]=(2.*(P1*P2*B12+P1*P3*B13+P2*P3*B23))*pow(apara,-2)*pow(aperp,-4);
//printf("%lf %lf %lf %lf (%lf,%lf,%lf) (%lf,%lf,%lf) (%lf,%lf,%lf)\n",k_1[i],k_2[i],k_3[i],B_theo0[i],P1,P2,B12,P1,P3,B13,P2,P3,B23);

}//loop over triangles



free(kreal_nomask);
free(Preal_nomask);

if( strcmp(bispectrum_BQ,"Q") == 0 && N_nomask>1){//only Q with mask

                N_nomask=Nlin;
                kreal_nomask = (double*) calloc(N_nomask, sizeof(double));
                Preal_nomask = (double*) calloc(N_nomask, sizeof(double));

      (*function_parameters).mode=0;//compute the monopole
      for(i=0;i<N_nomask;i++)
      {
         kreal_nomask[i]=Theory[i][0];
         (*function_parameters).kinput=kreal_nomask[i];
         adapt_integrate(1,integralP,function_parameters,1,XMIN1,XMAX1,0,precision,precision,&P1,&junk);
         Preal_nomask[i]=P1;
      }

     sprintf(type_of_analysis,"BAOISO");//we apply the mask as if it were an isotropic bao, i.e. only through W0, ignoring W2,...
     apply_mask(type_of_analysis,1, 0, 0, kreal_nomask,Preal_nomask, NULL, NULL,N_nomask, pos, W0, W2, W4, W6, W8, Nmask,spacing_mask, plan1, plan2, kreal_nomask[0], kreal_nomask[N_nomask-1], kmin_data0, kmax_data0, kmin_data0, kmax_data0, kmin_data0, kmax_data0,spacing_theory,0,1);
     
       for(i=0;i<NeffB0;i++)
       {
          if(i==0){
          P1=P_interpolLOG(k_1[i],kreal_nomask,Preal_nomask,N_nomask);
          P2=P_interpolLOG(k_2[i],kreal_nomask,Preal_nomask,N_nomask);
          P3=P_interpolLOG(k_3[i],kreal_nomask,Preal_nomask,N_nomask);
          }
          else
          {   //this is useful if they are ordered, but if they are not it's not a problem, it just does nothing
              if( k_1[i]!=k_1[i-1]){P1=P_interpolLOG(k_1[i],kreal_nomask,Preal_nomask,N_nomask);}
              if( k_2[i]!=k_2[i-1]){P2=P_interpolLOG(k_2[i],kreal_nomask,Preal_nomask,N_nomask);}
              if( k_3[i]!=k_3[i-1]){P3=P_interpolLOG(k_3[i],kreal_nomask,Preal_nomask,N_nomask);}
          }

         B_theo0[i]=B_theo0[i]/(P1*P2+P1*P3+P2*P3);

       }
     
free(kreal_nomask);
free(Preal_nomask);
}


}

void do_Btheo_iso(double B_theo0[], int NeffB0, double *parameters1,double *k_Plin,double *Plin,int Nlin, double *k_Olin, double *Olin, int NOlin,double *pos, double *W0, double *W2, double *W4, double *W6, double *W8, int Nmask, char *spacing_mask, char *path_to_mask1, double k_1[], double k_2[], double k_3[], int Npoly, fftw_plan plan1, fftw_plan plan2, double kmin_data0 , double kmax_data0, int wiggle,char *spacing_theory, double Sigma_smooth,double knl,double *k_n,double *n, double s8, char *bispectrum_BQ)
{
int i,j;
double func_1,func_2,func_3,olin_eff_1,olin_eff_2,olin_eff_3;
double P1,P2,P3,Z12,Z13,Z23;
double plin1,kused_1,plin2,kused_2,plin3,kused_3;
double kprime_1,kprime_2,kprime_3;
int interpolation_order,shiftN;
double alpha0, sigma_nl0;
int offset;
int Nlin1,Nlin2,Nlin3;
double w0lin_k1,w1lin_k1,w2lin_k1;
double w0lin_k2,w1lin_k2,w2lin_k2;
double w0lin_k3,w1lin_k3,w2lin_k3;
double beta_F,beta_G,beta_S,C1,C2,beta_mu;
double arg_12,arg_13,arg_23,Z_12,Z_13,Z_23;
double arg_1,arg_2,arg_3;
double *Pmask,*kmask;
char type_of_analysis[3000];
interpolation_order=1;

    if(interpolation_order==1){shiftN=1;}
    if(interpolation_order==2){shiftN=2;}


alpha0=parameters1[0];//always this way, even if multipoles are fitted, the input is aiso
sigma_nl0=parameters1[1];
beta_F=parameters1[2];
beta_G=parameters1[3];
beta_S=parameters1[4];
beta_mu=parameters1[5];
C1=parameters1[6];
C2=parameters1[7];
offset=8;



if(strcmp(path_to_mask1, "none") == 0){//no-mask

    for(j=0;j<NeffB0;j++)
    {

                  kused_1=k_1[j];
                  kused_2=k_2[j];
                  kused_3=k_3[j];
   
  Nlin1=determine_N_singlearray(k_Plin,kused_1,Nlin,spacing_theory);
  Nlin2=determine_N_singlearray(k_Plin,kused_2,Nlin,spacing_theory);
  Nlin3=determine_N_singlearray(k_Plin,kused_3,Nlin,spacing_theory);

if(interpolation_order==1){
w1lin_k1=determine_w1_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);
w1lin_k2=determine_w1_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);
w1lin_k3=determine_w1_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
}
if(interpolation_order==2){
w0lin_k1=determine_w0_2ndorder_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);
w1lin_k1=determine_w1_2ndorder_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);
w2lin_k1=determine_w2_2ndorder_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);

w0lin_k2=determine_w0_2ndorder_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);
w1lin_k2=determine_w1_2ndorder_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);
w2lin_k2=determine_w2_2ndorder_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);

w0lin_k3=determine_w0_2ndorder_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
w1lin_k3=determine_w1_2ndorder_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
w2lin_k3=determine_w2_2ndorder_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
}


plin1=P_interpol_fast(kused_1,Plin,Nlin,spacing_theory,interpolation_order,Nlin1,w0lin_k1,w1lin_k1,w2lin_k1);
plin2=P_interpol_fast(kused_2,Plin,Nlin,spacing_theory,interpolation_order,Nlin2,w0lin_k2,w1lin_k2,w2lin_k2);
plin3=P_interpol_fast(kused_3,Plin,Nlin,spacing_theory,interpolation_order,Nlin3,w0lin_k3,w1lin_k3,w2lin_k3);

       func_1=0;
       func_2=0;
       func_3=0;
       for(i=0;i<Npoly+1;i++)
       {
           if(i==0){func_1=plin1*parameters1[0+offset];func_2=plin2*parameters1[0+offset];func_3=plin3*parameters1[0+offset];}
           else{
                    func_1=func_1+parameters1[0+offset+i]*pow(kused_1,2-i);
                    func_2=func_2+parameters1[0+offset+i]*pow(kused_2,2-i);
                    func_3=func_3+parameters1[0+offset+i]*pow(kused_3,2-i);
           }
       }
          kprime_1=kused_1/alpha0;
          kprime_2=kused_2/alpha0;
          kprime_3=kused_3/alpha0;

  Nlin1=determine_N_singlearray(k_Olin,kprime_1,NOlin,spacing_theory);
  Nlin2=determine_N_singlearray(k_Olin,kprime_2,NOlin,spacing_theory);
  Nlin3=determine_N_singlearray(k_Olin,kprime_3,NOlin,spacing_theory);
if(interpolation_order==1){
w1lin_k1=determine_w1_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);
w1lin_k2=determine_w1_singlearray(k_Olin,kprime_2,Nlin2,spacing_theory);
w1lin_k3=determine_w1_singlearray(k_Olin,kprime_3,Nlin3,spacing_theory);

}
if(interpolation_order==2){
w0lin_k1=determine_w0_2ndorder_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);
w1lin_k1=determine_w1_2ndorder_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);
w2lin_k1=determine_w2_2ndorder_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);

w0lin_k2=determine_w0_2ndorder_singlearray(k_Olin,kprime_2,Nlin2,spacing_theory);
w1lin_k2=determine_w1_2ndorder_singlearray(k_Olin,kprime_2,Nlin2,spacing_theory);
w2lin_k2=determine_w2_2ndorder_singlearray(k_Olin,kprime_2,Nlin2,spacing_theory);

w0lin_k3=determine_w0_2ndorder_singlearray(k_Olin,kprime_3,Nlin3,spacing_theory);
w1lin_k3=determine_w1_2ndorder_singlearray(k_Olin,kprime_3,Nlin3,spacing_theory);
w2lin_k3=determine_w2_2ndorder_singlearray(k_Olin,kprime_3,Nlin3,spacing_theory);
}
olin_eff_1=P_interpol_fast(kprime_1,Olin,NOlin,spacing_theory,interpolation_order,Nlin1,w0lin_k1,w1lin_k1,w2lin_k1);
olin_eff_2=P_interpol_fast(kprime_2,Olin,NOlin,spacing_theory,interpolation_order,Nlin2,w0lin_k2,w1lin_k2,w2lin_k2);
olin_eff_3=P_interpol_fast(kprime_3,Olin,NOlin,spacing_theory,interpolation_order,Nlin3,w0lin_k3,w1lin_k3,w2lin_k3);

if(wiggle==0){olin_eff_1=1;olin_eff_2=1;olin_eff_3=1;}//wiggle=0 is only called at the very end for plotting reasons
//          arg_12=-0.5*sigma_nl0*sigma_nl0*(kused_1*kused_1+kused_2*kused_2);arg_12=exp(arg_12);
//          arg_13=-0.5*sigma_nl0*sigma_nl0*(kused_1*kused_1+kused_3*kused_3);arg_13=exp(arg_13);
//          arg_23=-0.5*sigma_nl0*sigma_nl0*(kused_2*kused_2+kused_3*kused_3);arg_23=exp(arg_23);

          arg_1=-0.5*sigma_nl0*sigma_nl0*(kused_1*kused_1);arg_1=exp(arg_1);
          arg_2=-0.5*sigma_nl0*sigma_nl0*(kused_2*kused_2);arg_2=exp(arg_2);
          arg_3=-0.5*sigma_nl0*sigma_nl0*(kused_3*kused_3);arg_3=exp(arg_3);

          P1=func_1*(1.+(olin_eff_1-1)*arg_1);
          P2=func_2*(1.+(olin_eff_2-1)*arg_2);
          P3=func_3*(1.+(olin_eff_3-1)*arg_3);
          
             Z12=Zeff(kused_1,kused_2,kused_3,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);
             Z13=Zeff(kused_1,kused_3,kused_2,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);
             Z23=Zeff(kused_2,kused_3,kused_1,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);

             //Z12=Zspt(kused_1,kused_2,kused_3,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);
             //Z13=Zspt(kused_1,kused_3,kused_2,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);
             //Z23=Zspt(kused_2,kused_3,kused_1,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);

if( strcmp(bispectrum_BQ,"B") == 0){
//            B_theo0[j]=P1*(P2*Z12*arg_12+P3*Z13*arg_13)+P2*P3*Z23*arg_23;
            B_theo0[j]=P1*(P2*Z12+P3*Z13)+P2*P3*Z23;
}
if( strcmp(bispectrum_BQ,"Q") == 0){
//            B_theo0[j]=(P1*(P2*Z12*arg_12+P3*Z13*arg_13)+P2*P3*Z23*arg_23)/(P1*P2*arg_12+P1*P3*arg_13+P2*P3*arg_23);
              B_theo0[j]=(P1*(P2*Z12+P3*Z13)+P2*P3*Z23)/(P1*P2+P1*P3+P2*P3);
}


    }
}//no-mask
else//mask
{

Pmask =   (double*) calloc( Nlin, sizeof(double));
kmask =   (double*) calloc( Nlin, sizeof(double));

     for(j=0;j<Nlin;j++)
     {

            kused_1=k_Plin[j];
            plin1=Plin[j];
       func_1=0;       
       for(i=0;i<Npoly+1;i++)
       {
           if(i==0){func_1=plin1*parameters1[0+offset];}
           else{
                    func_1=func_1+parameters1[0+offset+i]*pow(kused_1,2-i);
               }
       }
          kprime_1=kused_1/alpha0;

  Nlin1=determine_N_singlearray(k_Olin,kprime_1,NOlin,spacing_theory);
if(interpolation_order==1){
w1lin_k1=determine_w1_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);
}
if(interpolation_order==2){
w0lin_k1=determine_w0_2ndorder_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);
w1lin_k1=determine_w1_2ndorder_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);
w2lin_k1=determine_w2_2ndorder_singlearray(k_Olin,kprime_1,Nlin1,spacing_theory);
}

           olin_eff_1=P_interpol_fast(kprime_1,Olin,NOlin,spacing_theory,interpolation_order,Nlin1,w0lin_k1,w1lin_k1,w2lin_k1);
          if(wiggle==0){olin_eff_1=1;}//wiggle=0 is only called at the very end for plotting reason
          arg_12=-0.5*sigma_nl0*sigma_nl0*(kused_1*kused_1);arg_12=exp(arg_12);

          Pmask[j]=func_1*(1.+(olin_eff_1-1)*arg_12);//
//          Pmask[j]=func_1*(1.+(olin_eff_1-1))*arg_12;//old one
          kmask[j]=kused_1;

     }

//apply mask to Pmask
sprintf(type_of_analysis,"BAOISO");
apply_mask(type_of_analysis,1, 0, 0, kmask,Pmask, NULL, NULL,Nlin, pos, W0, W2, W4, W6, W8, Nmask,spacing_mask, plan1, plan2, kmask[0], kmask[Nlin-1], kmin_data0, kmax_data0, kmin_data0, kmax_data0, kmin_data0, kmax_data0,spacing_theory,0,1);

    for(j=0;j<NeffB0;j++)
    {

                  kused_1=k_1[j];
                  kused_2=k_2[j];
                  kused_3=k_3[j];

  Nlin1=determine_N_singlearray(k_Plin,kused_1,Nlin,spacing_theory);
  Nlin2=determine_N_singlearray(k_Plin,kused_2,Nlin,spacing_theory);
  Nlin3=determine_N_singlearray(k_Plin,kused_3,Nlin,spacing_theory);

if(interpolation_order==1){
w1lin_k1=determine_w1_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);
w1lin_k2=determine_w1_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);
w1lin_k3=determine_w1_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
}
if(interpolation_order==2){
w0lin_k1=determine_w0_2ndorder_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);
w1lin_k1=determine_w1_2ndorder_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);
w2lin_k1=determine_w2_2ndorder_singlearray(k_Plin,kused_1,Nlin1,spacing_theory);

w0lin_k2=determine_w0_2ndorder_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);
w1lin_k2=determine_w1_2ndorder_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);
w2lin_k2=determine_w2_2ndorder_singlearray(k_Plin,kused_2,Nlin2,spacing_theory);

w0lin_k3=determine_w0_2ndorder_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
w1lin_k3=determine_w1_2ndorder_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
w2lin_k3=determine_w2_2ndorder_singlearray(k_Plin,kused_3,Nlin3,spacing_theory);
}


P1=P_interpol_fast(kused_1,Pmask,Nlin,spacing_theory,interpolation_order,Nlin1,w0lin_k1,w1lin_k1,w2lin_k1);
P2=P_interpol_fast(kused_2,Pmask,Nlin,spacing_theory,interpolation_order,Nlin2,w0lin_k2,w1lin_k2,w2lin_k2);
P3=P_interpol_fast(kused_3,Pmask,Nlin,spacing_theory,interpolation_order,Nlin3,w0lin_k3,w1lin_k3,w2lin_k3);

             Z12=Zeff(kused_1,kused_2,kused_3,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);
             Z13=Zeff(kused_1,kused_3,kused_2,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);
             Z23=Zeff(kused_2,kused_3,kused_1,knl,k_n,n,Nlin-1,s8,beta_F,beta_G,beta_S,beta_mu,C1,C2);

//          arg_12=-0.5*sigma_nl0*sigma_nl0*(kused_1*kused_1+kused_2*kused_2);arg_12=exp(arg_12);
//          arg_13=-0.5*sigma_nl0*sigma_nl0*(kused_1*kused_1+kused_3*kused_3);arg_13=exp(arg_13);
//          arg_23=-0.5*sigma_nl0*sigma_nl0*(kused_2*kused_2+kused_3*kused_3);arg_23=exp(arg_23);

if( strcmp(bispectrum_BQ,"B") == 0){
            B_theo0[j]=P1*(P2*Z12+P3*Z13)+P2*P3*Z23;//exponentials already added in Ps
}
if( strcmp(bispectrum_BQ,"Q") == 0){
            B_theo0[j]=(P1*(P2*Z12+P3*Z13)+P2*P3*Z23)/(P1*P2+P1*P3+P2*P3);//exponentials already added in Ps
}

    }

free(Pmask);
free(kmask);
}



}

double Wth(double x)
{
double f;
f=3.*pow(x,-3)*(sin(x)-x*cos(x));
return f;
}

void integrals8(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;//cast from void to f_params
        double *k_L= params_function.k_Plin;
        double *Pk_L= params_function.Plin;
        int N= params_function.Nlin;

        double r=x[0];
        fval[0]=1./(2.*Pi*Pi)*r*r*P_interpolLOG(r,k_L,Pk_L,N)*Wth(r*8.)*Wth(r*8.);
}

double get_fid_sigma8(double *klin, double *Plin, int Nlin)
{
double s8,trash;
double precision=1e-6;
double XMIN[1]={klin[0]};
double XMAX[1]={klin[Nlin-1]};

        f_params *function_parameters;
        function_parameters = (f_params *) malloc(sizeof(f_params));        

       (*function_parameters).k_Plin=klin;
       (*function_parameters).Nlin=Nlin;
       (*function_parameters).Plin=Plin;


adapt_integrate(1, integrals8 , function_parameters, 1, XMIN, XMAX ,0, precision, precision, &s8, &trash);

free(function_parameters);

s8=sqrt(s8);
if(s8<=0){printf("Error computing fiducial sigma8, %lf. Exiting now...\n",s8);exit(0);}
return s8;
}

double get_fid_knl(double *klin, double *Plin, int Nlin)
{
int j;
double m,n;
double knl,trial_low,trial_high;
knl=-1;
        for(j=0;j<=Nlin;j++)
        {

           if(j==Nlin){break;}
           trial_high=pow(klin[j],3)*Plin[j]/(2.*Pi*Pi);
           if(trial_high>1){break;}//exit j-loop when crosses knl condition
       }
//printf("%lf %lf, %lf,%lf\n",klin[j],klin[j-1],pow(klin[j],3)*Plin[j]/(2.*Pi*Pi),pow(klin[j-1],3)*Plin[j-1]/(2.*Pi*Pi));
       if(j>0 && j<Nlin)//we should always within this
       {

          trial_low=pow(klin[j-1],3)*Plin[j-1]/(2.*Pi*Pi);
          m=(trial_low-trial_high)/(klin[j-1]-klin[j]);//printf("m=%lf\n",m);
          n=trial_high-m*klin[j];//printf("n=%lf\n",n);
          knl=(1.-n)/m;//printf("knl=%lf\n",knl);
       }

if(j==Nlin){knl=-1;}
if(j==0){knl=-2;}
return knl;
}

void generate_knl_array(double *klin, double *Plin, int Nlin, double *knl_y, double *s8_x, int Nknl, double s8ref)
{//this only runs one time so no need to be fast
int i,j;
double m,n;
double knl_min=klin[0];
double knl_max=klin[Nlin-1];//printf("%lf\n",knl_max);
double trial_low,trial_high;
double s8_min=0;//this should match the s8 priors
double s8_max=10;//this should match the s8 priors
double step=(s8_max-s8_min)/Nknl*1.;
double D2;


for(i=0;i<Nknl;i++)
{
  s8_x[i]=s8_min+(i)*step;//value of s8
  D2=pow(s8_x[i]/s8ref,2);//power spectrum re-scaling for the new s8 value

  if(s8_x[i]==0){knl_y[i]=knl_max;}//fixed value large enough
  else{
        for(j=0;j<=Nlin;j++)
        {
           if(j==Nlin){break;}
           //k3*P always grows.
           trial_high=pow(klin[j],3)*Plin[j]*D2/(2.*Pi*Pi);
           if(trial_high>1){break;}//exit j-loop when crosses knl condition
        }

       if(j>0 && j<Nlin)//we should always within this
       {
          
          trial_low=pow(klin[j-1],3)*Plin[j-1]*D2/(2.*Pi*Pi);
          m=(trial_low-trial_high)/(klin[j-1]-klin[j]);
          n=trial_high-m*klin[j];
          knl_y[i]=(1.-n)/m;
       }
       if(j==0)//s8[i] must be pretty large (knl has to be small)
       {
           knl_y[i]=knl_min;
       }
       if(j==Nlin)//s8[i] must be pretty small (knl has to  be prety large)
       {
           knl_y[i]=knl_max;
       }

}
}
//knl and s8 filled
}

void generate_n_array(double *klin, double *Plin, int Nlin, double *n_func, double *k_n_func)
{
int i;
for(i=0;i<Nlin-1;i++)
{
n_func[i]=( log10(Plin[i+1])-log10(Plin[i]))/( log10(klin[i+1])-log10(klin[i]) );
k_n_func[i]= pow(10, 0.5*log10(klin[i+1])+0.5*log10(klin[i]) );
//printf("%e %e\n",k_n_func[i],n_func[i]);
}

}

void smooth_n_array(double *k_in, int N_in, double *n_new, double *n_func, double *k_n_func, int Nfunc)
{
int i,n,j;
double *n_func_smooth,*k_n_func_smooth,*k_n_func_smooth_eff;
int *counts;
//n_new=malloc(sizeof(double)*N_in);
int n2;
n=50;//factor of smoothing (the larger less smoothing)
n_func_smooth=malloc(sizeof(double)*n);
k_n_func_smooth=malloc(sizeof(double)*n);
k_n_func_smooth_eff=malloc(sizeof(double)*n);
counts = malloc(sizeof(int)*n);
//smooth original n_func into k_n_func_smooth
for(i=0;i<n;i++){
k_n_func_smooth[i]=pow(10, (log10(k_in[N_in-1])-log10(k_in[0]))/n*1.*i+log10(k_in[0]) );
counts[i]=0;
n_func_smooth[i]=0;
}
n2=0;
for(i=1;i<n;i++){

  for(j=1;j<Nfunc;j++)
  {

      if(k_n_func[j]>=k_n_func_smooth[i-1] && k_n_func[j]<k_n_func_smooth[i])
      {
        k_n_func_smooth_eff[i]=k_n_func_smooth_eff[i]+log10(k_n_func[j]);
        n_func_smooth[i]=n_func_smooth[i]+n_func[j];
        counts[i]++;
      }
  }

}

j=1;
for(i=1;i<n;i++){

if(counts[i]>0){
k_n_func_smooth_eff[j]=pow(10,k_n_func_smooth_eff[i]/counts[i]);
n_func_smooth[j]=n_func_smooth[i]/counts[i];
//printf("%e %e\n",k_n_func_smooth_eff[j],n_func_smooth[j]);
j++;
}

}
free(counts);
free(k_n_func_smooth);
n2=j;

//final, cast new smooth array into k_in vector
for(i=0;i<N_in;i++)
{
n_new[i]=P_interpol_extrapol(k_in[i],k_n_func_smooth_eff,n_func_smooth,n2);//it extrapolates below and beyond
//printf("%e %e\n",k_in[i],n_new[i]);
}
free(k_n_func_smooth_eff);
free(n_func_smooth);
}

double interpol_f(double z, double **f,int i)
{
double g,a,a05,a1,a2;
a=1./(1+z);
a05=1./1.5;
a1=1./2.;
a2=1./3.;

//2nd order Lagrange interpolation formula
////http://www-classes.usc.edu/engr/ce/108/lagrange.pdf
g=f[i][0]+(a-a05)*( f[i][0]-f[i][1])/(a05-a1)+(a-a05)*(a-a1)/(a05-a2)*( (f[i][0]-f[i][1])/(a05-a1) - (f[i][1] - f[i][2])/(a1-a2) );

return g;
}

double geo_factor(double z, double A, double ctheta_min, double ctheta_med, double ctheta_max)
{
double geo;

int i;
double **f;//f[5][3];
double Anorm=0.001;// in (h/Mpc)^2
 f = (double **) calloc(5, sizeof(double*));
 for(i=0;i<5;i++){f[i]= (double *) calloc(3, sizeof(double));}

double feff[5];
//redshift 0.5
f[0][0]=1.0033;
f[1][0]=-0.0040;
f[2][0]=0.0240;
f[3][0]=-0.0568;
f[4][0]=0.01325;

//redshift 1.0
f[0][1]=1.0180;
f[1][1]=-0.0041;
f[2][1]=0.0149;
f[3][1]=-0.0547;
f[4][1]=0.011144;

//redshift 2.0
f[0][2]=1.037;
f[1][2]=0.0020;
f[2][2]=0.0056;
f[3][2]=-0.048;
f[4][2]=0.0081;

//interpol
for(i=0;i<5;i++)
{
feff[i]=interpol_f(z,f,i);
}

geo=feff[0]+feff[1]*ctheta_med/ctheta_max+feff[2]*ctheta_min/ctheta_max+feff[3]*A/Anorm+feff[4]*A*A/(Anorm*Anorm);

freeTokens(f, 5);
return geo;
}

