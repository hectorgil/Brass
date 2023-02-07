#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void get_lineMask(FILE *f, double params[])
{
double sav,seff,w0,w2,w4,w6,w8;

fscanf(f,"%lf %lf %lf %lf %lf %lf %lf\n",&sav,&seff,&w0,&w2,&w4,&w6,&w8);

params[0]=sav;
params[1]=seff;
params[2]=w0;
params[3]=w2;
params[4]=w4;
params[5]=w6;
params[6]=w8;

}
void get_lineB(FILE *f, double params[],int type)
{
double k1cen,k1eff,k2cen,k2eff,k3cen,k3eff,b0,q0,NoiseB,NoiseQ;

k1cen=0;k2cen=0;k3cen=0;
k1eff=0;k2eff=0;k3eff=0;
b0=0;q0=0;
NoiseB=0;NoiseQ=0;

if(type==0)//DATA B(k)
{
fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %*f\n",&k1cen,&k1eff,&k2cen,&k2eff,&k3cen,&k3eff,&b0,&NoiseB,&q0,&NoiseQ);
}

if(type==1)//Mocks/Covariance B(k)
{
fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %*f\n",&k1cen,&k1eff,&k2cen,&k2eff,&k3cen,&k3eff,&b0,&NoiseB,&q0,&NoiseQ);
}

params[0]=k1cen;
params[1]=k1eff;

params[2]=k2cen;
params[3]=k2eff;

params[4]=k3cen;
params[5]=k3eff;

params[6]=b0;
params[7]=NoiseB;

params[8]=q0;
params[9]=NoiseQ;
}

void get_lineP(FILE *f, double params[],int type)
{

double kcen,keff,p0,p2,p4,noiseP;

kcen=0;
keff=0;
p0=0;p2=0;p4=0;noiseP=0;
if(type==0)//DATA P(k)
{
fscanf(f,"%lf %lf %lf %lf %lf %*d %lf\n",&kcen,&keff,&p0,&p2,&p4,&noiseP);
}

if(type==1)//Mocks/Covariance P(k)
{
fscanf(f,"%lf %lf %lf %lf %lf %*d %lf\n",&kcen,&keff,&p0,&p2,&p4,&noiseP);
}

params[0]=kcen;
params[1]=keff;
params[2]=p0;
params[3]=p2;
params[4]=p4;
params[5]=noiseP;

}

void get_linealphas(FILE *f, double params[],int type)
{

double p0;

p0=0;
if(type==0)//DATA P(k)
{
fscanf(f,"%lf\n",&p0);
}

if(type==1)//Mocks/Covariance P(k)
{
fscanf(f,"%lf\n",&p0);
}

params[0]=0;
params[1]=0;
params[2]=p0;
params[3]=0;
params[4]=0;
params[5]=0;

}


void getBAOcov(double *paramsBAO, char *path_to_data2_bao)
{
paramsBAO[0]=-2;
paramsBAO[1]=-2;
paramsBAO[2]=-2;

int N=2+8;//number of total fitted parameters
int i;
FILE *f;
f=fopen(path_to_data2_bao,"r");

fscanf(f,"%*s %*s\n");
for(i=0;i<N;i++){fscanf(f,"%*s\n");}
fscanf(f,"%*s %*s\n");

fscanf(f,"%lf %lf ",&paramsBAO[0],&paramsBAO[2]);for(i=2;i<N;i++){fscanf(f,"%*s ");}
fscanf(f,"\n");
fscanf(f,"%*f %lf ",&paramsBAO[1]);
fclose(f);

paramsBAO[2]=paramsBAO[2]/sqrt(paramsBAO[0]*paramsBAO[1]);

//printf("%lf %lf %lf\n",paramsBAO[0],paramsBAO[1],paramsBAO[2]);

if(paramsBAO[2]>=1 || paramsBAO[2]<=-1){printf("Error, cross apara-aperp out of statistical range: %lf\n",paramsBAO[2]);exit(0);}
if(paramsBAO[0]<=0 || paramsBAO[1]<=0){printf("Error, variance of apara,aperp out of statistical range: %lf %lf\n",paramsBAO[0],paramsBAO[1]);exit(0);}

//paramBAO[0] -> sigma_apara^2
//paramBAO[1] -> sigma_aperp^2
//paramBAO[2] -> r_apara_aperp
}
