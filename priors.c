#define Pi (4.*atan(1.))
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//IMPORTANT
//** You should not modify the structure of the functions in this file, just fill in with numerical value for the priors etc you want **/

void set_mcmc_priors(double alpha_min, double alpha_max, double params_low[], double params_high[],int N,char *type_of_analysis, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum, char *Sigma_def_type, char *Sigma_independent, double *Sigma_type, char *local_b2s2, char *local_b3nl, char *fog_free, char *fogmodel_ps, char *fog_bs, char *sigma8_free, char *RSD_fit, int Nchunks, int Npolynomial, char *path, char *id)
{
int i,j;
int bao=0;
int rsd=0;
double params_low_in[1000];
double params_high_in[1000];
int modeP0,modeP2,modeP4,modeB0;
char printout[10];
int Nalphas=0;
FILE *f;
char file_name[200];
sprintf(file_name,"%s/Variable_order_%s.txt",path,id);

//option to print the variables selected in order in the mcmc and log files. 
sprintf(printout,"yes");
//sprintf(printout,"no");


if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0){bao=1;}

if( strcmp( type_of_analysis,"FS") == 0 || strcmp( type_of_analysis,"FSBAOISO") == 0 || strcmp( type_of_analysis,"FSBAOANISO") == 0 || strcmp( type_of_analysis,"FSalphasrecon") == 0){

if( strcmp( RSD_fit,"yes") == 0 || strcmp( RSD_fit,"shape") == 0 || strcmp( RSD_fit,"shape2") == 0  ){rsd=1;}

}

if( strcmp(printout,"yes") == 0){
f=fopen(file_name,"w");
if(f==NULL){printf("Error writting the variable order file. Exiting now...\n");exit(0);}
}

modeP0=0;
modeP2=0;
modeP4=0;
modeB0=0;
if( strcmp(do_power_spectrum,"yes") == 0 && bao==1){
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}
}
if( strcmp(do_bispectrum,"yes") == 0 && bao==1){modeB0=1;}

if( strcmp(type_of_analysis,"BAOISO") == 0)
{
Nalphas=1;
if( strcmp(do_power_spectrum,"yes") == 0 && modeP0+modeP2+modeP4>1){Nalphas=2;}
if( strcmp(do_bispectrum,"yes") == 0 && strcmp(do_power_spectrum,"yes") == 0 && modeP0+modeP2+modeP4==1 && modeP0==0){Nalphas=2;}

}else// FS, FSBAISO, FSBAOANISO, BAOANISO, FSalphasrecon
{
Nalphas=2;
if( strcmp(do_bispectrum,"yes") == 0 && strcmp(do_power_spectrum,"no") == 0 ){Nalphas=1;}//it says baniso but it does iso for bis always
if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0){Nalphas=2;}//corrects above line
if(bao == 0 && rsd ==0){Nalphas=0;}//corrects above line
}

j=0;

if( Nalphas == 1)
{

if(modeP0 == 1 || modeB0 == 1){
params_low_in[j]=alpha_min; params_high_in[j]=alpha_max;j++;//a0 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a0\n",j-1);}
}

if(modeP2 == 1){
params_low_in[j]=alpha_min; params_high_in[j]=alpha_max;j++;//a2 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a2\n",j-1);}
}

if(modeP4 == 1){
params_low_in[j]=alpha_min; params_high_in[j]=alpha_max;j++;//a4 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a4\n",j-1);}
}

}else{

if( bao == 1 || rsd == 1){
params_low_in[j]=alpha_min; params_high_in[j]=alpha_max;j++;//apara
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a_para\n",j-1);}

params_low_in[j]=alpha_min;  params_high_in[j]=alpha_max;j++;//aperp
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a_perp\n",j-1);}
}

}


if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0)
{

if( strcmp(RSD_fit,"shape") == 0 || strcmp(RSD_fit,"shape2") == 0){
params_low_in[j]=-1.0; params_high_in[j]=1.0;j++;//m shape-slope (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d m_BGV\n",j-1);}
}

if( strcmp(RSD_fit,"shape2") == 0){
params_low_in[j]=-1.0; params_high_in[j]=1.0;j++;//m2 shape-slope (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d m2_BGV\n",j-1);}
}

if( strcmp( RSD_fit, "yes") ==0 || strcmp( RSD_fit, "shape") ==0 || strcmp( RSD_fit, "shape2") ==0){
params_low_in[j]=-20; params_high_in[j]=20;j++;//f growth factor (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d f\n",j-1);}
}


if( strcmp(sigma8_free,"yes") == 0){
params_low_in[j]=0; params_high_in[j]=10;j++;//sigma8 (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d s8\n",j-1);}
}


}



if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0){


if( strcmp(Sigma_def_type,"effective") == 0)
{

if( Sigma_type[0]>0 && modeP0==1){
params_low_in[j]=0; params_high_in[j]=20;j++;//sigmaP eff0
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP0_eff\n",j-1);}
}

if( Sigma_type[1]>0 && modeP2==1){
params_low_in[j]=0; params_high_in[j]=20;j++;//sigmaP eff2
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP2_eff\n",j-1);}
}

if( Sigma_type[2]>0 && modeP4==1){
params_low_in[j]=0; params_high_in[j]=20;j++;//sigmaP eff4
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP4_eff\n",j-1);}
}

}
else{//para-perp

if( Sigma_type[3]>0 && modeP0+modeP2+modeP4>0){
params_low_in[j]=0; params_high_in[j]=20;j++;//sigmaP para
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP_para\n",j-1);}
}
//independent yes 
if( Sigma_type[4]>0 && strcmp(Sigma_independent,"yes") ==0 && modeP0+modeP2+modeP4>0)
{
params_low_in[j]=0; params_high_in[j]=20;j++;//sigmaP perp
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP_perp\n",j-1);}
}


}

if( strcmp(do_power_spectrum,"yes") == 0 )
{
if( strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0)
{
params_low_in[j]=0;params_high_in[j]=30;j++;//beta_eff (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d beta_eff\n",j-1);}

params_low_in[j]=0; params_high_in[j]=20;j++;//Baniso NGC (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d B NGC\n",j-1);}

}

}

if(modeP0==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{

params_low_in[j]=0; params_high_in[j]=20;j++;//B P0 NGC  (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp0 NGC\n",j-1);}

}


for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i P0 NGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap0_%d NGC\n",j-1,i);}
}

}

if(modeP2==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
params_low_in[j]=0; params_high_in[j]=20;j++;//Biso P2 NGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp2 NGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i P2 NGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap2_%d NGC\n",j-1,i);}
}

}

if(modeP4==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
params_low_in[j]=0; params_high_in[j]=20;j++;//Biso P4 NGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp4 NGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i P4 NGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap4_%d NGC\n",j-1,i);}
}

}


//if Nchuncks=2
if(Nchunks==2)
{

if( strcmp(do_power_spectrum,"yes") == 0 )
{
if( strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0)
{

params_low_in[j]=0; params_high_in[j]=20;j++;//Baniso SGC (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d B SGC\n",j-1);}

}

}

if(modeP0==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{

params_low_in[j]=0; params_high_in[j]=20;j++;//B P0 SGC  (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp0 SGC\n",j-1);}

}


for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i P0 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap0_%d SGC\n",j-1,i);}
}

}

if(modeP2==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
params_low_in[j]=0; params_high_in[j]=20;j++;//Biso P2 SGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp2 SGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i P2 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap2_%d SGC\n",j-1,i);}
}

}

if(modeP4==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
params_low_in[j]=0; params_high_in[j]=20;j++;//Biso P4 SGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp4 SGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i P4 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap4_%d SGC\n",j-1,i);}
}

}

}//Nchunks2


//if bispectrum yes
if( strcmp(do_bispectrum,"yes") == 0)
{

if( Sigma_type[5]>0)
{
params_low_in[j]=0; params_high_in[j]=20;j++;//sigmaB eff0
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaB0_eff\n",j-1);}
}


params_low_in[j]=-10; params_high_in[j]=10;j++;//betaF NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaF NGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//betaG NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaG NGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//betaS NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaS NGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//betamu NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betamu NGC\n",j-1);}

params_low_in[j]=-5; params_high_in[j]=5;j++;//C1 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C1 NGC\n",j-1);}

params_low_in[j]=-5; params_high_in[j]=5;j++;//C2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C2 NGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//Bbis NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bb0 NGC\n",j-1);}

for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i B0 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ab0_%d NGC\n",j-1,i);}
}

if(Nchunks==2)
{

params_low_in[j]=-10; params_high_in[j]=10;j++;//betaF SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaF SGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//betaG SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaG SGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//betaS SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaS SGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//betamu SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betamu SGC\n",j-1);}

params_low_in[j]=-5; params_high_in[j]=5;j++;//C1 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C1 SGC\n",j-1);}

params_low_in[j]=-5; params_high_in[j]=5;j++;//C2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C2 SGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//Bbis SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bb0 SGC\n",j-1);}

for(i=0;i<Npolynomial;i++){
params_low_in[j]=-100000; params_high_in[j]=100000;j++;//A_i B0 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ab0_%d SGC\n",j-1,i);}
}



}//Nchunks
}//bispectrum

}//bao



if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0)
{

params_low_in[j]=0; params_high_in[j]=10;j++;//b1 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b1 NGC\n",j-1);}

params_low_in[j]=-20; params_high_in[j]=20;j++;//b2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2 NGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//A NGC (amplitude of shot noise)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Anoise NGC\n",j-1);}

if( strcmp( local_b2s2, "no") == 0 ){
params_low_in[j]=-20; params_high_in[j]=20;j++;//b2s2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2s2 NGC\n",j-1);}
}

if( strcmp( local_b3nl, "no") == 0 ){
params_low_in[j]=-20; params_high_in[j]=20;j++;//b3nl NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b3nl NGC\n",j-1);}
}

if( strcmp( fog_free, "yes") == 0 ){
params_low_in[j]=0; params_high_in[j]=10;j++;//sigmafog_P NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_P NGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 )
{
params_low_in[j]=0; params_high_in[j]=10;j++;//sigmafog_B NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_B NGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0 )
{
params_low_in[j]=0; params_high_in[j]=10;j++;//avir NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d avir NGC\n",j-1);}
}

if(Nchunks==2)
{

params_low_in[j]=0; params_high_in[j]=10;j++;//b1 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b1 SGC\n",j-1);}

params_low_in[j]=-20; params_high_in[j]=20;j++;//b2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2 SGC\n",j-1);}

params_low_in[j]=-10; params_high_in[j]=10;j++;//A SGC (amplitude of shot noise)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Anoise SGC\n",j-1);}

if( strcmp( local_b2s2, "no") == 0 ){
params_low_in[j]=-20; params_high_in[j]=20;j++;//b2s2 SGC 
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2s2 SGC\n",j-1);}
}

if( strcmp( local_b3nl, "no") == 0 ){
params_low_in[j]=-20; params_high_in[j]=20;j++;//b3nl SGC 
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b3nl SGC\n",j-1);}
}

if( strcmp( fog_free, "yes") == 0 ){
params_low_in[j]=0; params_high_in[j]=10;j++;//sigmafog_P SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_P SGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 )
{
params_low_in[j]=0; params_high_in[j]=10;j++;//sigmafog_B SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_B SGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0 )
{
params_low_in[j]=0; params_high_in[j]=10;j++;//avir SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d avir SGC\n",j-1);}
}

}

}//FS


if( strcmp(printout,"yes") == 0){fclose(f);}

if(N!=j){printf("Error in set_mcmc_priors, number of free elements reported by do_bao or do_rsd functions (%d) does not match with the number of elements found by prior functions (%d). Exiting now...\n",N,j);exit(0);}

//exit(0);//remove after test


for(i=0;i<N;i++)
{
params_low[i]=params_low_in[i];
params_high[i]=params_high_in[i];
}

}

void set_proposal_mean(double mean[], int N, char *type_of_analysis, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum, char *Sigma_def_type, char *Sigma_independent, double *Sigma_type, char *local_b2s2, char *local_b3nl, char *fog_free, char *fogmodel_ps, char *fog_bs, char *sigma8_free, char *RSD_fit, int Nchunks, int Npolynomial,char *path, char *id)
{
int bao=0;
int rsd=0;
int i,j;
double mean_in[1000];
int modeP0,modeP2,modeP4,modeB0;
int Nalphas=0;
FILE *f;
char file_name[200];
char printout[10];
sprintf(file_name,"%s/Variable_order_%s.txt",path,id);

//option to print the variables selected in order in the mcmc and log files.
sprintf(printout,"yes");
//sprintf(printout,"no");

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0){bao=1;}

if( strcmp( type_of_analysis,"FS") == 0 || strcmp( type_of_analysis,"FSBAOISO") == 0 || strcmp( type_of_analysis,"FSBAOANISO") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0){

if( strcmp( RSD_fit,"yes") == 0 || strcmp( RSD_fit,"shape") == 0 || strcmp( RSD_fit,"shape2") == 0  ){rsd=1;}
}

if( strcmp(printout,"yes") == 0){
f=fopen(file_name,"w");
if(f==NULL){printf("Error writting the variable order file. Exiting now...\n");exit(0);}
}

modeP0=0;
modeP2=0;
modeP4=0;
modeB0=0;

if( strcmp(do_power_spectrum,"yes") == 0 && bao==1){
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}
}
if( strcmp(do_bispectrum,"yes") == 0 && bao==1){modeB0=1;}




if( strcmp(type_of_analysis,"BAOISO") == 0)
{
Nalphas=1;
if( strcmp(do_power_spectrum,"yes") == 0 && modeP0+modeP2+modeP4>1){Nalphas=2;}
if( strcmp(do_bispectrum,"yes") == 0 && strcmp(do_power_spectrum,"yes") == 0 && modeP0+modeP2+modeP4==1 && modeP0==0){Nalphas=2;}

}else// FS, FSBAISO, FSBAOANISO, BAOANISO, FSalphasrecon
{
Nalphas=2;
if( strcmp(do_bispectrum,"yes") == 0 && strcmp(do_power_spectrum,"no") == 0 ){Nalphas=1;}//it says baniso but it does iso for bis always
if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0){Nalphas=2;}//corrects above line
if(bao == 0 && rsd == 0){Nalphas=0;}
}

j=0;

if( Nalphas == 1)
{

if(modeP0 == 1 || modeB0 == 1){
mean_in[j]=0.987230;j++;//a0 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a0\n",j-1);}
}

if(modeP2 == 1){
mean_in[j]=1;j++;//a2 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a2\n",j-1);}
}

if(modeP4 == 1){
mean_in[j]=1;j++;//a4 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a4\n",j-1);}
}

}else{

if(bao == 1 || rsd == 1){

mean_in[j]=1.003877e+00;j++;//apara 
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a_para\n",j-1);}

mean_in[j]=9.987614e-01;j++;//aperp 
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a_perp\n",j-1);}
}

}

if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0 ||  strcmp(type_of_analysis,"FSalphasrecon") == 0)
{

if( strcmp( RSD_fit, "shape") ==0 || strcmp( RSD_fit, "shape2") ==0){
mean_in[j]=-6.331823e-02;j++;//mBGV factor (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d m_BGV\n",j-1);}
}

if( strcmp( RSD_fit, "shape2") ==0){
mean_in[j]=0;j++;//m2BGV factor (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d m2_BGV\n",j-1);}
}

if( strcmp( RSD_fit, "yes") ==0 || strcmp( RSD_fit, "shape") ==0 || strcmp( RSD_fit, "shape2") ==0){
mean_in[j]=8.138821e-01;j++;//f growth factor (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d f\n",j-1);}
}


if( strcmp(sigma8_free,"yes") == 0){
mean_in[j]=0.8;j++;//sigma8 (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d s8\n",j-1);}
}


}

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0){


if( strcmp(Sigma_def_type,"effective") == 0)
{

if( Sigma_type[0]>0 && modeP0==1){
mean_in[j]=7.7;j++;//sigmaP eff0
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP0_eff\n",j-1);}
}

if( Sigma_type[1]>0 && modeP2==1){
mean_in[j]=6.02;j++;//sigmaP eff2
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP2_eff\n",j-1);}
}

if( Sigma_type[2]>0 && modeP4==1){
mean_in[j]=6.03;j++;//sigmaP eff4
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP4_eff\n",j-1);}
}

}
else{//para-perp

if( Sigma_type[3]>0 && modeP0+modeP2+modeP4>0){
mean_in[j]=6.04;j++;//sigmaP para
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP_para\n",j-1);}
}

if( Sigma_type[4]>0 && strcmp(Sigma_independent,"yes") ==0 && modeP0+modeP2+modeP4>0)
{
mean_in[j]=6.05;j++;//sigmaP perp
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP_perp\n",j-1);}
}


}

if( strcmp(do_power_spectrum,"yes") == 0 )
{
if( strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0)
{
mean_in[j]=4.042171e-01;j++;//beta_eff (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d beta_eff\n",j-1);}

mean_in[j]=3.910621e+00;j++;//Baniso NGC (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d B NGC\n",j-1);}

}

}

if(modeP0==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{

mean_in[j]=4.434856;j++;//B P0 NGC  (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp0 NGC\n",j-1);}

}


for(i=0;i<Npolynomial;i++){

if(i==0){mean_in[j]=3.699983e+03;j++;}//A_i P0 NGC (BAOISO & BAOANISO)
if(i==1){mean_in[j]=-9.889215e+02;j++;}//A_i P0 NGC (BAOISO & BAOANISO)
if(i==2){mean_in[j]=6.139992e+01;j++;}//A_i P0 NGC (BAOISO & BAOANISO)
if(i==3){mean_in[j]=-4.295601e-01;j++;}//A_i P0 NGC (BAOISO & BAOANISO)
if(i==4){mean_in[j]=3.608886e-02;j++;}//A_i P0 NGC (BAOISO & BAOANISO)
if(i>4){mean_in[j]=0;j++;}

if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap0_%d NGC\n",j-1,i);}
}

}

if(modeP2==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
mean_in[j]=1.01;j++;//Biso P2 NGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp2 NGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
//mean_in[j]=0.001*(i+1);j++;//A_i P2 NGC (BAOISO & BAOANISO)
if(i==0){mean_in[j]=5.401990e+03;j++;}//A_i P2 NGC (BAOISO & BAOANISO)
if(i==1){mean_in[j]=-3.378355e+03;j++;}//A_i P2 NGC (BAOISO & BAOANISO)
if(i==2){mean_in[j]=1.052580e+02;j++;}//A_i P2 NGC (BAOISO & BAOANISO)
if(i==3){mean_in[j]=-2.334719e+00;j++;}//A_i P2 NGC (BAOISO & BAOANISO)
if(i==4){mean_in[j]=7.426446e-02;j++;}//A_i P2 NGC (BAOISO & BAOANISO)
if(i>4){mean_in[j]=0;j++;}

if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap2_%d NGC\n",j-1,i);}
}

}

if(modeP4==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
mean_in[j]=1.0005;j++;//Biso P4 NGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp4 NGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
mean_in[j]=0.0001*(i+1);j++;//A_i P4 NGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap4_%d NGC\n",j-1,i);}
}

}

if(Nchunks==2)
{

if( strcmp(do_power_spectrum,"yes") == 0 )
{

if( strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0)
{
mean_in[j]=1.123;j++;//Baniso SGC (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d B SGC\n",j-1);}

}

}

if(modeP0==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{

mean_in[j]=2.0;j++;////B P0 SGC  (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp0 SGC\n",j-1);}

}


for(i=0;i<Npolynomial;i++){
mean_in[j]=1.0*(i+1);j++;//A_i P0 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap0_%d SGC\n",j-1,i);}
}

}

if(modeP2==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
mean_in[j]=1.0;j++;//Biso P2 SGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp2 SGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
mean_in[j]=0.00001*(i+1);j++;//A_i P2 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap2_%d SGC\n",j-1,i);}
}

}

if(modeP4==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
mean_in[j]=1.0;j++;//Biso P4 SGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp4 SGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
mean_in[j]=0.000001*(i+1);j++;//A_i P4 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap4_%d SGC\n",j-1,i);}
}

}

}//Nchunks2

if( strcmp(do_bispectrum,"yes") == 0)
{

if( Sigma_type[5]>0)
{
mean_in[j]=7.409052;j++;//sigmaB eff0
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaB0_eff\n",j-1);}
}

mean_in[j]=0.685928;j++;//betaF NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaF NGC\n",j-1);}

mean_in[j]=0.078665;j++;//betaG NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaG NGC\n",j-1);}

mean_in[j]=-0.463318;j++;//betaS NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaS NGC\n",j-1);}

mean_in[j]=-0.346577;j++;//betamu NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betamu NGC\n",j-1);}

mean_in[j]=-0.121952;j++;//C1 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C1 NGC\n",j-1);}

mean_in[j]=1.234003;j++;//C2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C2 NGC\n",j-1);}

mean_in[j]=-1.240867;j++;//Bbis NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bb0 NGC\n",j-1);}

for(i=0;i<Npolynomial;i++){

if(i==0){mean_in[j]=14027.389238;j++;}//A_i B0 NGC
if(i==1){mean_in[j]=-8792.631358;j++;}//A_i B0 NGC
if(i==2){mean_in[j]=2075.363426;j++;}//A_i B0 NGC
if(i==3){mean_in[j]=-6.905314;j++;}//A_i B0 NGC
if(i==4){mean_in[j]=-0.202867;j++;}//A_i B0 NGC
if(i>4){mean_in[j]=0;j++;}

if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ab0_%d NGC\n",j-1,i);}
}

if(Nchunks==2)
{

mean_in[j]=-0.1751;j++;//betaF SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaF SGC\n",j-1);}

mean_in[j]=0.077;j++;//betaG SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaG SGC\n",j-1);}

mean_in[j]=0.357;j++;//betaS SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaS SGC\n",j-1);}

mean_in[j]=-0.15;j++;//betamu SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betamu SGC\n",j-1);}

mean_in[j]=0.226;j++;//C1 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C1 SGC\n",j-1);}

mean_in[j]=0.53;j++;//C2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C2 SGC\n",j-1);}

mean_in[j]=3.28;j++;//Bbis SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bb0 SGC\n",j-1);}

for(i=0;i<Npolynomial;i++){
//mean_in[j]=4.+0.01*(i+1);j++;//A_i B0 SGC
if(i==0){mean_in[j]=-0.052;j++;}//A_i B0 NGC
if(i==1){mean_in[j]=0.271;j++;}//A_i B0 NGC
if(i==2){mean_in[j]=-0.106;j++;}//A_i B0 NGC
if(i==3){mean_in[j]=0.0265;j++;}//A_i B0 NGC
if(i==4){mean_in[j]=0.0176;j++;}//A_i B0 NGC
if(i>4){mean_in[j]=0;j++;}

if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ab0_%d SGC\n",j-1,i);}
}

}//Nchunks
}//bispectrum

}//bao

if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0)
{

mean_in[j]=1.965136e+00;j++;//b1 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b1 NGC\n",j-1);}

mean_in[j]=-6.895879e-01;j++;//b2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2 NGC\n",j-1);}

mean_in[j]=1.405268e+00;j++;// A NGC (amplitude of shot noise)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Anoise NGC\n",j-1);}

if( strcmp( local_b2s2, "no") == 0 ){
mean_in[j]=0.1;j++;// b2s2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2s2 NGC\n",j-1);}
}

if( strcmp( local_b3nl, "no") == 0 ){
mean_in[j]=0.2;j++;// b3nl NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b3nl NGC\n",j-1);}
}

if( strcmp( fog_free, "yes") == 0 ){
mean_in[j]=3.899547e+00;j++;// sigmafog_P NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_P NGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 )
{
mean_in[j]=6.0;j++;// sigmafog_B NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_B NGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0 )
{
mean_in[j]=2.1;j++;// avir NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d avir NGC\n",j-1);}
}

if(Nchunks==2)
{

mean_in[j]=2.372098;j++;//b1 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b1 SGC\n",j-1);}

mean_in[j]=1.655197;j++;//b2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2 SGC\n",j-1);}

mean_in[j]=1.093171;j++;//A SGC (amplitude of shot noise)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Anoise SGC\n",j-1);}

if( strcmp( local_b2s2, "no") == 0 ){
mean_in[j]=0.11;j++;//b2s2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2s2 SGC\n",j-1);}
}

if( strcmp( local_b3nl, "no") == 0 ){
mean_in[j]=0.22;j++;//b3nl SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b3nl SGC\n",j-1);}
}

if( strcmp( fog_free, "yes") == 0 ){
mean_in[j]=5.389390;j++;//sigmafog_P SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_P SGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 )
{
mean_in[j]=6.1;j++;//sigmafog_B SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_B SGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0 )
{
mean_in[j]=2.1;j++;// avir SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d avir SGC\n",j-1);}
}

}

}//FS

if( strcmp(printout,"yes") == 0){fclose(f);}

if(N!=j){printf("Error in set_proposal_mean, number of free elements reported by do_bao or do_rsd functions (%d) does not match with the number of elements found by prior functions (%d). Exiting now...\n",N,j);exit(0);}


for(i=0;i<N;i++){mean[i]=mean_in[i];}

}


void set_proposal_error(double error[],int N, char *type_of_analysis, char *fit_BAO, char *do_power_spectrum, char *do_bispectrum, char *Sigma_def_type, char *Sigma_independent, double *Sigma_type, char *local_b2s2, char *local_b3nl, char *fog_free, char *fogmodel_ps, char *fog_bs, char *sigma8_free, char *RSD_fit, int Nchunks, int Npolynomial,char *path, char *id)
{
int bao=0;
int rsd=0;
int i,j;
double error_in[1000];
double times;
times=0.01;

int modeP0,modeP2,modeP4,modeB0;
int Nalphas=0;
FILE *f;
char file_name[200];
char printout[10];
sprintf(file_name,"%s/Variable_order_%s.txt",path,id);

//option to print the variables selected in order in the mcmc and log files.
sprintf(printout,"yes");
//sprintf(printout,"no");

if( strcmp(printout,"yes") == 0){
f=fopen(file_name,"w");
if(f==NULL){printf("Error writting the variable order file. Exiting now...\n");exit(0);}
}

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0){bao=1;}

if( strcmp( type_of_analysis,"FS") == 0 || strcmp( type_of_analysis,"FSBAOISO") == 0 || strcmp( type_of_analysis,"FSBAOANISO") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0){

if( strcmp( RSD_fit,"yes") == 0 || strcmp( RSD_fit,"shape") == 0 || strcmp( RSD_fit, "shape2") ==0){rsd=1;}
}

modeP0=0;
modeP2=0;
modeP4=0;
modeB0=0;

if( strcmp(do_power_spectrum,"yes") == 0 && bao==1){
if(strcmp(fit_BAO, "P0") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP0=1;}
if(strcmp(fit_BAO, "P2") == 0 || strcmp(fit_BAO, "P02") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP2=1;}
if(strcmp(fit_BAO, "P4") == 0 || strcmp(fit_BAO, "P04") == 0 || strcmp(fit_BAO, "P24") == 0 || strcmp(fit_BAO, "P024") == 0 ){modeP4=1;}
}
if( strcmp(do_bispectrum,"yes") == 0 && bao==1){modeB0=1;}


if( strcmp(type_of_analysis,"BAOISO") == 0)
{
Nalphas=1;
if( strcmp(do_power_spectrum,"yes") == 0 && modeP0+modeP2+modeP4>1){Nalphas=2;}
if( strcmp(do_bispectrum,"yes") == 0 && strcmp(do_power_spectrum,"yes") == 0 && modeP0+modeP2+modeP4==1 && modeP0==0){Nalphas=2;}

}else// FS, FSBAISO, FSBAOANISO, BAOANISO, FSalphas
{
Nalphas=2;
if( strcmp(do_bispectrum,"yes") == 0 && strcmp(do_power_spectrum,"no") == 0 ){Nalphas=1;}//it says baniso but it does iso for bis always
if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0){Nalphas=2;}//corrects above line
if(bao == 0 && rsd ==0){Nalphas=0;}
}

j=0;

if(  Nalphas == 1)
{

if(modeP0 == 1 || modeB0 == 1){
error_in[j]=0.001759*times;j++;//a0 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a0\n",j-1);}
}

if(modeP2 == 1){
error_in[j]=times;j++;//a2 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a2\n",j-1);}
}

if(modeP4 == 1){
error_in[j]=times;j++;//a4 for some BAOISO options
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a4\n",j-1);}
}

}else{

if( bao == 1 || rsd == 1){

error_in[j]=3.862645e-02*times;j++;//apara
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a_para\n",j-1);}

error_in[j]=2.728535e-02*times;j++;//aperp
if( strcmp(printout,"yes") == 0){fprintf(f,"%d a_perp\n",j-1);}
}

}

if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0)
{

if( strcmp( RSD_fit, "shape") ==0 || strcmp( RSD_fit, "shape2") ==0 ){
error_in[j]=3.295162e-02*times;j++;//mBGV  (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d m_BGV\n",j-1);}
}

if( strcmp( RSD_fit, "shape2") ==0 ){
error_in[j]=times;j++;//m2BGV (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d m2_BGV\n",j-1);}
}

if( strcmp( RSD_fit, "yes") ==0 || strcmp( RSD_fit, "shape") ==0 || strcmp( RSD_fit, "shape2") ==0){
error_in[j]=1.050844e-01*times;j++;//f growth factor (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d f\n",j-1);}
}


if( strcmp(sigma8_free,"yes") == 0){
error_in[j]=times;j++;//sigma8 (for FS or FS+BAO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d s8\n",j-1);}
}


}

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0){


if( strcmp(Sigma_def_type,"effective") == 0)
{

if( Sigma_type[0]>0 && modeP0==1){
error_in[j]=times;j++;//sigmaP eff0
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP0_eff\n",j-1);}
}

if( Sigma_type[1]>0 && modeP2==1){
error_in[j]=times;j++;//sigmaP eff2
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP2_eff\n",j-1);}
}

if( Sigma_type[2]>0 && modeP4==1){
error_in[j]=times;j++;//sigmaP eff4
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP4_eff\n",j-1);}
}

}
else{//para-perp

if( Sigma_type[3]>0 && modeP0+modeP2+modeP4>0){
error_in[j]=times;j++;//sigmaP para
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP_para\n",j-1);}
}

if( Sigma_type[4]>0 && strcmp(Sigma_independent,"yes") ==0 && modeP0+modeP2+modeP4>0)
{
error_in[j]=times;j++;//sigmaP perp
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaP_perp\n",j-1);}
}


}

if( strcmp(do_power_spectrum,"yes") == 0 )
{
if( strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0)
{
error_in[j]=times;j++;//beta_eff (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d beta_eff\n",j-1);}

error_in[j]=times;j++;//Baniso NGC (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d B NGC\n",j-1);}

}
}
if(modeP0==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{

error_in[j]=times;j++;//B P0 NGC  (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp0 NGC\n",j-1);}

}


for(i=0;i<Npolynomial;i++){

error_in[j]=10*times;j++;//A_i P0 NGC (BAOISO & BAOANISO)

if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap0_%d NGC\n",j-1,i);}
}

}

if(modeP2==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
error_in[j]=times;j++;//Biso P2 NGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp2 NGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
error_in[j]=10*times;j++;//A_i P2 NGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap2_%d NGC\n",j-1,i);}
}

}

if(modeP4==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
error_in[j]=times;j++;//Biso P4 NGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp4 NGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
error_in[j]=times;j++;//A_i P4 NGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap4_%d NGC\n",j-1,i);}
}

}

if(Nchunks==2)
{

if( strcmp(do_power_spectrum,"yes") == 0 )
{
if( strcmp(type_of_analysis,"BAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0)
{
error_in[j]=times;j++;//Baniso SGC (BAOANISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d B SGC\n",j-1);}

}

}

if(modeP0==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{

error_in[j]=times;j++;////B P0 SGC  (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp0 SGC\n",j-1);}

}


for(i=0;i<Npolynomial;i++){
error_in[j]=times;j++;//A_i P0 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap0_%d SGC\n",j-1,i);}
}

}

if(modeP2==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
error_in[j]=times;j++;//Biso P2 SGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp2 SGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
error_in[j]=times;j++;//A_i P2 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap2_%d SGC\n",j-1,i);}
}

}

if(modeP4==1)
{

if( strcmp(type_of_analysis,"BAOISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0)
{
error_in[j]=times;j++;//Biso P4 SGC (BAOISO only)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bp4 SGC\n",j-1);}
}

for(i=0;i<Npolynomial;i++){
error_in[j]=times;j++;//A_i P4 SGC (BAOISO & BAOANISO)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ap4_%d SGC\n",j-1,i);}
}

}

}//Nchunks2


if( strcmp(do_bispectrum,"yes") == 0)
{

if( Sigma_type[5]>0)
{
error_in[j]=times;j++;//sigmaB eff0
if( strcmp(printout,"yes") == 0){fprintf(f,"%d SigmaB0_eff\n",j-1);}
}

error_in[j]=0.078585*times;j++;//betaF NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaF NGC\n",j-1);}

error_in[j]=0.004860*times;j++;//betaG NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaG NGC\n",j-1);}

error_in[j]=0.460783*times;j++;//betaS NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaS NGC\n",j-1);}

error_in[j]=0.115251*times;j++;//betamu NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betamu NGC\n",j-1);}

error_in[j]=0.153846*times;j++;//C1 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C1 NGC\n",j-1);}

error_in[j]=0.358202*times;j++;//C2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C2 NGC\n",j-1);}

error_in[j]=0.767593*times;j++;//Bbis NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bb0 NGC\n",j-1);}

for(i=0;i<Npolynomial;i++){
//error_in[j]=times;j++;//A_i B0 NGC
if(i==0){error_in[j]=3040.630980*times;j++;}
if(i==1){error_in[j]=1875.452100*times;j++;}
if(i==2){error_in[j]=424.377316*times;j++;}
if(i==3){error_in[j]=3.370433*times;j++;}
if(i==4){error_in[j]=0.048477*times;j++;}
if(i>4){error_in[j]=times;j++;}

if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ab0_%d NGC\n",j-1,i);}
}

if(Nchunks==2)
{

error_in[j]=times;j++;//betaF SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaF SGC\n",j-1);}

error_in[j]=times;j++;//betaG SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaG SGC\n",j-1);}


error_in[j]=times;j++;//betaS SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betaS SGC\n",j-1);}

error_in[j]=times;j++;//betamu SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d betamu SGC\n",j-1);}

error_in[j]=times;j++;//C1 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C1 SGC\n",j-1);}

error_in[j]=times;j++;//C2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d C2 SGC\n",j-1);}

error_in[j]=times;j++;//Bbis SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Bb0 SGC\n",j-1);}

for(i=0;i<Npolynomial;i++){
error_in[j]=times;j++;//A_i B0 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Ab0_%d SGC\n",j-1,i);}
}

}//Nchunks
}//bispectrum

}//bao

if( strcmp(type_of_analysis,"FS") == 0 || strcmp(type_of_analysis,"FSBAOANISO") == 0 || strcmp(type_of_analysis,"FSBAOISO") == 0 || strcmp(type_of_analysis,"FSalphasrecon") == 0)
{

error_in[j]=7.412931e-03*times;j++;//b1 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b1 NGC\n",j-1);}

error_in[j]=2.916427e-01*times;j++;//b2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2 NGC\n",j-1);}

error_in[j]=5.470813e-02*times;j++;// A NGC (amplitude of shot noise)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Anoise NGC\n",j-1);}

if( strcmp( local_b2s2, "no") == 0 ){
error_in[j]=times;j++;// b2s2 NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2s2 NGC\n",j-1);}
}

if( strcmp( local_b3nl, "no") == 0 ){
error_in[j]=times;j++;// b3nl NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b3nl NGC\n",j-1);}
}

if( strcmp( fog_free, "yes") == 0 ){
error_in[j]=5.258373e-02*times;j++;// sigmafog_P NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_P NGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 )
{
error_in[j]=times;j++;// sigmafog_B NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_B NGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0 )
{
error_in[j]=0.01*times;j++;// avir NGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d avir NGC\n",j-1);}
}

if(Nchunks==2)
{

error_in[j]=8.258670e-02*times;j++;//b1 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b1 SGC\n",j-1);}

error_in[j]=1.477681e+00*times;j++;//b2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2 SGC\n",j-1);}

error_in[j]=7.000780e-02*times;j++;//A SGC (amplitude of shot noise)
if( strcmp(printout,"yes") == 0){fprintf(f,"%d Anoise SGC\n",j-1);}

if( strcmp( local_b2s2, "no") == 0 ){
error_in[j]=times;j++;//b2s2 SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b2s2 SGC\n",j-1);}
}

if( strcmp( local_b3nl, "no") == 0 ){
error_in[j]=times;j++;//b3nl SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d b3nl SGC\n",j-1);}
}

if( strcmp( fog_free, "yes") == 0 ){
error_in[j]=4.412008e-01*times;j++;//sigmafog_P SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_P SGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fog_bs, "no") == 0 && strcmp(do_bispectrum, "yes") == 0 && strcmp(do_power_spectrum, "yes") == 0 )
{
error_in[j]=times;j++;//sigmafog_B SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d sigma_fog_B SGC\n",j-1);}
}

if(strcmp(fog_free, "yes") == 0 && strcmp(fogmodel_ps, "Exponential_avir") == 0 )
{
error_in[j]=0.01*times;j++;// avir SGC
if( strcmp(printout,"yes") == 0){fprintf(f,"%d avir SGC\n",j-1);}
}

}

}//FS


if( strcmp(printout,"yes") == 0){fclose(f);}

if(N!=j){printf("Error in set_proposal_error, number of free elements reported by do_bao or do_rsd functions (%d) does not match with the number of elements found by prior functions (%d). Exiting now...\n",N,j);exit(0);}


for(i=0;i<N;i++){error[i]=error_in[i];}

}
