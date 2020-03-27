#define Pi (4.*atan(1.))
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void set_mcmc_priors(double alpha_min, double alpha_max, double params_low[], double params_high[],int N)
{
int i;
double params_low_in[1000];
double params_high_in[1000];

//P0+P2 fixed sigma NGC+SGC
//params_low_in[0]=alpha_min;params_high_in[0]=alpha_max;

//params_low_in[1]=alpha_min;params_high_in[1]=alpha_max;

params_low_in[0]=0.5;params_high_in[0]=1.5;

params_low_in[1]=0.5;params_high_in[1]=1.5;

params_low_in[2]=0;params_high_in[2]=10;

//params_low_in[2]=0;params_high_in[2]=30;//sigmapara

//params_low_in[3]=0;params_high_in[3]=30;//sigmaperp

params_low_in[3]=0;params_high_in[3]=30;//beta

params_low_in[4]=0; params_high_in[4]=20;//B NGC

params_low_in[5]=-20000; params_high_in[5]=20000;//A1NGC-P0

params_low_in[6]=-20000; params_high_in[6]=20000;//A2NGC-P0

params_low_in[7]=-20000; params_high_in[7]=20000;//A3NGC-P0

params_low_in[8]=-20000; params_high_in[8]=20000;// A1 NGC P2

params_low_in[9]=-20000; params_high_in[9]=20000;//A2 NGC P2

params_low_in[10]=-20000; params_high_in[10]=20000;// A3 NGC P2

params_low_in[11]=0;params_high_in[11]=30;//b1 NGC

params_low_in[12]=-10;params_high_in[12]=10;//b2 NGC

params_low_in[13]=-5; params_high_in[13]=5;//A NGC

params_low_in[14]=0; params_high_in[14]=20;//sigmaP NGC



for(i=0;i<N;i++)
{
params_low[i]=params_low_in[i];
params_high[i]=params_high_in[i];
}

}

void set_proposal_mean(double mean[], int N)
{
//Recall ordering

//For BAO-non-Bispectrum
//for P0,P2,P4 with sigma fixed: alpha024, B, A1, A2,....,Bsgc,A1sgc,...
//for P0,P2,P4 with sigma free: alpha024,sigma,B,A1,A2,....,Bsgc,A1sgc,A2sgc,....
//for P02,P24,P04,P024 with sigma fixed: alpha_para,alpha_perp, B,A1,A2,....Bsgc,A1sgc,....
//for P02,P24,P04,P024 with sigma free: alpha_para,alpha_perp,sigma,B,A1,A2,....Bsgc,A1sgc,....
//
//For BAO-with-Bispectrum
//
//For RSD
//b1,b2,ANoise,sigma_FOG_P,f,alpha_para,alpha_perp,s8 (1chunk,s8-free, only-PS, local lagrangian) 
//b1,b2,ANoise,sigma_FOG_P,b1SGC,b2SGC,AnoiseSGC,sigma_FOG_PSGC,f,alpha_para,alpha_perp,s8 (2chunks, s8-free, onlyPS, local lagrangian)
//b1,b2,ANoise,sigma_FOG_P,b2s2,b3nl,b1SGC,b2SGC,AnoiseSGC,sigma_FOG_PSGC,b2s2SGC,b3nlSGC,f,alpha_para,alpha_perp,s8 (2chunks, s8-free, onlyPS, non-local lagrangian)
///b1,b2,ANoise,sigma_FOG_P,sigma_FOG_B,b2s2,b3nl,b1SGC,b2SGC,AnoiseSGC,sigma_FOG_PSGC,sigma_FOG_BSGC,b2s2SGC,b3nlSGC,f,alpha_para,alpha_perp,s8 (2chunks, s8-free, with Bispectrum and sigmaFogBis allowed, non-local lagrangian)
//No-rsd f=constant...
//Number maximum of free parameters: 18.
int i;
double mean_in[1000];

//9.462094e-01    3.098746e+03    -1.470274e+03   3.126687e+02    1.108556e+00    4.691189e+03    -2.453760e+03   4.397504e+02    1.220000e+01     9.950000e-01    1.010000e+00    3.260518e+01

//P02+NGC+SGC

mean_in[0]= 1.0;
mean_in[1]= 1.0;
mean_in[2]= 0.74;
mean_in[3]= 0.350876;
mean_in[4]= 4.319223;
mean_in[5]= 2791.484000;
mean_in[6]= -290.749500;
mean_in[7]= -129.247100;
mean_in[8]= 3588.262000;
mean_in[9]= -2712.177000;
mean_in[10]= 79.526880;
mean_in[11]= 1.981132;
mean_in[12]= 0.699853;
mean_in[13]= 1.002701;
mean_in[14]= 3.822350;


/*
mean_in[0]= 1.006212;
mean_in[1]= 1.000509;
//mean_in[2]= 12.289970;
//mean_in[3]= 7.352326;
mean_in[2]= 0.333623;
mean_in[3]= 3.725070;
mean_in[4]= 824.228100;
mean_in[5]= -425.596300;
mean_in[6]= 252.139500;
mean_in[7]= -131.571100;
mean_in[8]= -698.828700;
mean_in[9]= 87.067340;
mean_in[10]= 3.664607;
mean_in[11]= 294.036200;
mean_in[12]= -231.411100;
mean_in[13]= 200.949300;
mean_in[14]= -329.709500;
mean_in[15]= -510.818600;
mean_in[16]= 47.709560;
*/
for(i=0;i<N;i++){mean[i]=mean_in[i];}

}


void set_proposal_error(double error[],int N)
{
//same ordering that for set_proposal_mean
int i;
double error_in[1000];
double times;
times=0.1;

//P0+P2 NGC+SGC fixed sigma

error_in[0]= 0.001882*times;
error_in[1]= 0.001025*times;
error_in[2]= 0.006596*times;
error_in[3]= 0.005371*times;
error_in[4]= 0.040218*times;
error_in[5]= 138.160356*times;
error_in[6]= 83.901961*times;
error_in[7]= 19.352462*times;
error_in[8]= 225.216223*times;
error_in[9]=  84.118760*times;
error_in[10]= 6.492040*times;
error_in[11]= 0.006369*times;
error_in[12]= 0.206511*times;
error_in[13]= 0.041540*times;
error_in[14]= 0.045204*times;

/*
error_in[0]=0.001*times;
error_in[1]=0.001*times;

//error_in[2]=0.146496*times;
//error_in[3]=0.118276*times;

error_in[2]=0.007896*times;

error_in[3]=0.024578*times;

error_in[4]=35.621597*times;
error_in[5]=24.712764*times;
error_in[6]= 7.163482*times;

error_in[7]=75.572501*times;
error_in[8]=51.875696*times;
error_in[9]=15.036457*times;

error_in[10]=0.027438*times;
error_in[11]=46.530735*times;
error_in[12]=31.705435*times;
error_in[13]= 9.197304*times;

error_in[14]=101.420516*times;
error_in[15]=58.456889*times;
error_in[16]=15.196578*times;
*/

for(i=0;i<N;i++){error[i]=error_in[i];}

}
