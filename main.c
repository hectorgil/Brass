/*
RUSTICO Rapid foUrier STatIstics COde
RUSTICOX Rapid foUrier STatIstics COde for X-correlations
All rights reserved
Author: Hector Gil Marin
Date: 1st May 2020
email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "mass_assignment.h"
#include "read_positions.h"
#include "functions.h"
#include "fftw_compute.h"
#include "ps_write.h"
//#include "order_algorithm.h"
#include "bispectrum.h"
#include "mask.h"
//#include "structures.h"
/*
typedef struct{
          double OMEGA_M;
}f_params;
*/

int main(int argc, char *argv[])
{
FILE *f,*g;
int i,tid;
char type_of_code[200];//rustico vs rusticoX
char type_of_input[200];//input particles or density
char name_file_ini[2000];//name of inizialization file
//Main Parameters read in
double L1,L2;//Box limit in Mpc/h
double kf,kny;
char type_of_survey[50];//type of survey: Periodic or Cutsky or PeriodicFKP
char type_of_computation[10];//type of computation: DSY /  FFT / DSE
int power_grid;//power in number of grid cells: integer between 6 and 15 
int ngrid;//Number of grid cells: pow(2,power_grid)
char binning_type[10];//Type of binning for the power spectrum: linear or log
double bin_ps;//Size of the bin for the power spectrum in h/Mpc
int mubin;//number of bins for mu
char do_mu_bins[10];//perform mu bining (only for periodic boxes)
char do_anisotropy[10];// compute P2,P4 and maybe P1 P3
char do_odd_multipoles[10];//compute P1,P3
char write_kvectors[10];
char file_for_mu[10];
char header[10];//write header yes/no
char type_of_file[10];//ascii or gadget(only for periodic boxes)
char type_of_fileB[10];//ascii or gadget(only for periodic boxes)
int gadget_files,gadget_filesB;//number of gadget files per realization
int snapshot_num;//individual gadget file
char RSD[10];//RSD distorsion on Gadget boxes?
char RSDB[10];//RSD distorsion on Gadget boxes?
double Rsmoothing;
double I22delta,I33delta,Pnoisedelta,Bnoise1delta,Bnoise2delta;
char output_density[200];//output density in a grid (only for FFT). The output is delta(x) for a periodic box, and F(x) =  n(x)-alpha*n_ran(x), where n(x) and n_ran(x) are the number density of galaxies and randoms. Note that F is not normalized by I2^0.5 or I3^0.333. 

//Bispectrum parameters
char do_bispectrum[10];//Do Bispectrum at all?
char do_bispectrum2[10];//Do bispectrum multipoles?
double Deltakbis;//Triangle Bin size in terms of k-fundamental
double kmin,kmax;//Minimum and maximum k
char triangles_num[2000];//FFT,APR_SUM,APR_EFF,EXA_EFF
char write_triangles[2000];
char path_for_triangles[2000];
char triangles_id[2000];
char do_multigrid[10];
char bispectrum_optimization[2000];
char triangle_shapes[10];

//Read inout parameters
char name_data_in[2000];//Path and name of data source
char name_dataB_in[2000];//Path and name of data source

char name_gadget_file[2000];//Full path for gadget like file
char name_gadget_fileB[2000];//Full path for gadget like file

long int Ndata;//Number of lines of data source
long int Ndata2;//Number of lines of data used
double   Ndata2w;//Number of weighted particles used

long int NdataB;//Number of lines of data source
long int Ndata2B;//Number of lines of data used
double   Ndata2Bw;//Number of weighted particles used
    
char name_randoms_in[2000];//Path and name of random source
char name_randomsB_in[2000];//Path and name of random source

long int Nrand;//Number of lines of random source
long int Nrand2;//Number of lines of random used

    long int NrandB;//Number of lines of random source
    long int Nrand2B;//Number of lines of random used
    
char name_path_out[2000];//Path where to write output
char name_id[2000];//String to identify output
char name_ps_out[2000];//Final name of the output for the power spectrum
char name_ps_kvectors[2000];
char name_psBB_out[2000];//Final name of the output for the power spectrum
char name_psAB_out[2000];//Final name of the output for the power spectrum
    
char name_ps_out2[2000];//Final name of the output for the power spectrum for mu bin filess
char name_psBB_out2[2000];//Final name of the output for the power spectrum
char name_psAB_out2[2000];//Final name of the output for the power spectrum

char name_bs_out[2000];//Final name of the output for the bispectrum
char name_bsBBB_out[2000];//Final name of the output for the bispectrum

char name_bsAAB_out[2000];//Final name of the output for the bispectrum
char name_bsBBA_out[2000];//Final name of the output for the bispectrum

char name_bsBAB_out[2000];//Final name  of the output for the bispectrum
char name_bsABA_out[2000];//Final name of the output for the bispectrum

char name_bsABB_out[2000];//Final name of the output for the bispectrum
char name_bsBAA_out[2000];//Final name of the output for the bispectrum

//Bispectrum quadrupole for x-correlations may be too much. Don't allow for it for now. Also need to work out the shot noise
char name_bs002_out[2000];//Final name of the output for the bispectrum002
char name_bs020_out[2000];//Final name of the output for the bispectrum020
char name_bs200_out[2000];//Final name of the output for the bispectrum200

char name_den_out[2000];//Final name of the output for the density
char name_wink_out[2000];//Final name of the output for the power spectrum
char name_winkBB_out[2000];//Final name of the output for the power spectrum
char name_winkAB_out[2000];//Final name of the output for the power spectrum


//FFT parameters
char type_of_mass_assigment[10];//Type of mass assignment: NGC, CIC, TSC, PCS, P4S, P5S
int Ninterlacing;//Number of interlacing steps (1 for no-interlacing steps)
char grid_correction_string[10];// Do Grid Correction: yes/no
int mode_correction;//Correction factor power
char type_of_yamamoto[20];

//Cutsky parameters
double z_min,z_max;//Minimum and maximum (excluding) redshift cuts
double z_minB,z_maxB;//Minimum and maximum (excluding) redshift cuts
double Omega_m;//Value of Omega matter;
double Omega_L;
double speed_of_light=299792.458;
double Area_survey;//Value of Area of the survey in deg^2
double Area_surveyB;//Value of Area of the survey in deg^2
char Hexadecapole_type[20];//L4 or L2L2 or L1L3
char Octopole_type[20];//L3 or L1L2
char Quadrupole_type[20];//L2 or L1L1
char type_normalization_mode[20];//Normalize according to the area of the survey or the n(z) value: Area / Density
char type_normalization_mode2[20];//Normalize according to the n(z) of data or randoms file
double Shot_noise_factor;// Factor between 0 and 1. 0 Correspond ....
char shuffle[10];//Generate n(z) for the randoms
char window_function[10];//Compute Window Selection function
char write_shuffled_randoms[10];
char name_out_randoms[200];
char name_out_density[200];
char name_out_densityB[200];
int window_norm_bin;
double deltaS_window;
double percentage_randoms_window;
char yamamoto4window[10];
int reverse;

//Positions pointers
  double* pos_x;
  double* pos_y;
  double* pos_z;
  double* weight;
  double* radata;
  double* decdata;
  double* zdata;
  double* wcoldata;
  double* wsysdata;
  double* wfkpdata;
  double* nzdata;
 

  double* pos_x_rand;
  double* pos_y_rand;
  double* pos_z_rand;
  double* weight_rand;
    
    //Positions pointers B
     double* pos_xB;
     double* pos_yB;
     double* pos_zB;
     double* weightB;
     double* radataB;
     double* decdataB;
     double* zdataB;
     double* wcoldataB;
     double* wsysdataB;
     double* wfkpdataB;
     double* nzdataB;
    

     double* pos_x_randB;
     double* pos_y_randB;
     double* pos_z_randB;
     double* weight_randB;
    
//  double DeltaR;  

 //Maximum and minimum values for positions of galaxies and randoms.
  double max,min;
  double maxB,minB;

//Parameters relative to shot noise, effective redshifts, number of particles and normalization
double Psn_1a, Psn_1b, Psn_2a,Psn_2b, z_effective_data,I_norm_data,z_effective_rand,I_norm_rand,alpha_data,alpha_data1,alpha_rand,alpha_rand1,alpha,alpha1,I22,I_norm_data2,I_norm_data3,I_norm_data4, I_norm_rand2,I_norm_rand3,I_norm_rand4;
    double Psn_1aB, Psn_1bB, Psn_2aB,Psn_2bB, z_effective_dataB,I_norm_dataB,z_effective_randB,I_norm_randB,alpha_dataB,alpha_data1B,alpha_randB,alpha_rand1B,alphaB,alpha1B,I22B,I_norm_data2B,I_norm_data3B,I_norm_data4B, I_norm_rand2B,I_norm_rand3B,I_norm_rand4B;
double P_shot_noise1,P_shot_noise2,P_shot_noise;
double P_shot_noise1B,P_shot_noise2B,P_shot_noiseB;

double Bsn1a,Bsn2a,Bsn1b,Bsn2b,Bsn1,Bsn2,IN1,IN2,IN11,IN22,I3_norm_data,I3_norm_data2,I3_norm_data3,I3_norm_data4,I3_norm_rand,I3_norm_rand2,I3_norm_rand3,I3_norm_rand4,I33,Bsn,IN;

double Bsn1aB,Bsn2aB,Bsn1bB,Bsn2bB,Bsn1B,Bsn2B,IN1B,IN2B,IN11B,IN22B,I3_norm_dataB,I3_norm_data2B,I3_norm_data3B,I3_norm_data4B,I3_norm_randB,I3_norm_rand2B,I3_norm_rand3B,I3_norm_rand4B,I33B,BsnB,INB;
    
double num_effective,num_effective_rand;
double num_effective2,num_effective2_rand;
double num_effective3,num_effective3_rand;
    
    double num_effectiveB,num_effective_randB;
    double num_effective2B,num_effective2_randB;
    double num_effective3B,num_effective3_randB;

    int n_lines_parallel;//Number of parallel threads



//Read Inizialization parameters
sprintf(name_file_ini,argv[1]);
f=fopen(name_file_ini,"r");
if(f==NULL){printf("File %s not found...\t exiting now\n",name_file_ini);return 0;}
else{printf("Reading Inizialization file: %s\n\n",name_file_ini);}
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_code);
    if(strcmp(type_of_code, "rustico") != 0 && strcmp(type_of_code, "rusticoX") != 0){printf("Error, type of code has to be either 'rustico' or 'rusticoX'. Input read: %s. Exiting now...\n",type_of_code);exit(0);}
    if(strcmp(type_of_code, "rustico") == 0){printf("Code runninng as auto - delta.\n");}
    if(strcmp(type_of_code, "rusticoX") == 0){printf("Code runninng as cross - delta.\n");}


fscanf(f,"%*s %*s %*s %*s %s\n",type_of_survey);
    if(strcmp(type_of_code, "rustico") == 0){fscanf(f,"%*s %*s %*s %*s %s\n",type_of_file); if( strcmp(type_of_file,"ascii") !=0 && strcmp(type_of_file,"gadget") !=0 ){printf("Error, type of file must be either ascii or gadget. Input read: %s. Exiting now...\n",type_of_file);exit(0);} }
    if(strcmp(type_of_code, "rusticoX") == 0){fscanf(f,"%*s %*s %*s %*s %s %s\n",type_of_file,type_of_fileB);

if( strcmp(type_of_file,"ascii") !=0 && strcmp(type_of_file,"gadget") !=0 ){printf("Error, type of file A must be either ascii or gadget. Input read: %s. Exiting now...\n",type_of_file);exit(0);}
if( strcmp(type_of_fileB,"ascii") !=0 && strcmp(type_of_fileB,"gadget") !=0 ){printf("Error, type of file B must be either ascii or gadget. Input read: %s. Exiting now...\n",type_of_fileB);exit(0);}

}

    fscanf(f,"%*s %*s %*s %*s %s\n",type_of_input);
    if(strcmp(type_of_input, "particles") != 0 && strcmp(type_of_input, "density") != 0){printf("Error, type of input has to be either 'particles' or 'density'. Input read: %s. Exiting now...\n",type_of_input);exit(0);}
    if(strcmp(type_of_input, "density") == 0 && strcmp(type_of_file,"gadget") == 0 ){ printf("Error, a density type of input requires an ascii type of file. Exiting now...\n");exit(0); }
    if(strcmp(type_of_input, "density") == 0 && strcmp(type_of_fileB,"gadget") == 0 && strcmp(type_of_code,"rusticoX") == 0 ){printf("Error, a density type of input requires an ascii type of file. Exiting now...\n");exit(0);}
    if(strcmp(type_of_code, "rustico") == 0){fscanf(f,"%*s %*s %*s %*s %d\n",&gadget_files);}
    if(strcmp(type_of_code, "rusticoX") == 0){fscanf(f,"%*s %*s %*s %*s %d %d\n",&gadget_files,&gadget_filesB);}
    if(strcmp(type_of_code, "rustico") == 0){fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %s\n",RSD);}
    if(strcmp(type_of_code, "rusticoX") == 0){fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %s %s\n",RSD,RSDB);}
printf("== Inizialization Parameters ==\n");
if(strcmp(type_of_code, "rustico") == 0){
printf("Type of Survey: %s\n",type_of_survey);
printf("Type of File: %s\n",type_of_file);
if(strcmp(type_of_file, "gadget") == 0){printf("Number of gadget files: %d\n",gadget_files);}
if(strcmp(type_of_file, "gadget") == 0){printf("RSD distortion: %s\n",RSD);}
}

if(strcmp(type_of_code, "rusticoX") == 0){
     printf("Type of Survey: %s\n",type_of_survey);
     printf("Type of File A: %s\n",type_of_file);
     printf("Type of File B: %s\n",type_of_fileB);
     if(strcmp(type_of_file , "gadget") == 0){printf("Number of gadget files: %d\n",gadget_files);}
     if(strcmp(type_of_file, "gadget") == 0){printf("RSD distortion: %s\n",RSD);}
     if(strcmp(type_of_fileB , "gadget") == 0){printf("Number of gadget files: %d\n",gadget_filesB);}
     if(strcmp(type_of_fileB, "gadget") == 0){printf("RSD distortion: %s\n",RSDB);}
}
     fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&L1,&L2);
     printf("Box edges at %lf Mpc/h and %lf Mpc/h; Size of the Box %lf Mpc/h\n",L1,L2,L2-L1);
     fscanf(f,"%*s %*s %*s %*s %s\n",type_of_computation);
     printf("Type of Computation: %s\n",type_of_computation);
 

fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n",binning_type);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %lf\n",&bin_ps);
fscanf(f,"%*s %*s %*s %*s %lf %lf\n\n",&kmin,&kmax);
printf("Binning for the Power Spectrum: %s; Size of bin: %lf; k-range %lf < k[h/Mpc] < %lf\n\n",binning_type,bin_ps,kmin,kmax);
fscanf(f,"%*s %*s %*s %*s %s\n",do_anisotropy);
fscanf(f,"%*s %*s %*s %*s %s\n",do_odd_multipoles);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",write_kvectors);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",do_mu_bins);
fscanf(f,"%*s %*s %*s %*s %d\n",&mubin);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",file_for_mu);
if( strcmp(do_mu_bins, "yes") == 0){printf("Mu-binned in %d parts\n",mubin);}
if( strcmp(do_mu_bins, "no") == 0){mubin=1;}

fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",do_bispectrum);
fscanf(f,"%*s %*s %*s %*s %s\n",do_bispectrum2);
fscanf(f,"%*s %*s %*s %s\n",do_multigrid);
fscanf(f,"%*s %*s %*s %s\n",bispectrum_optimization);
fscanf(f,"%*s %*s %*s %s\n",triangle_shapes);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %lf\n",&Deltakbis);
fscanf(f,"%*s %*s %*s %s\n",triangles_num);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",write_triangles);
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n\n",path_for_triangles);

printf("== Bispectrum Parameters ==\n");
printf("Do Bispectrum? %s\n",do_bispectrum);
printf("Do Bispectrum multipoles? %s\n",do_bispectrum2);
if(strcmp(do_bispectrum, "yes") == 0){printf("Bispectrum bins %lf\nMultigrid Computation:%s\nBispectrum Optimization:%s\nTriangle Shapes: %s\nTriangle normalization: %s\nWrite individual Triangles:%s \n\n",Deltakbis,do_multigrid,bispectrum_optimization,triangle_shapes,triangles_num,write_triangles);}

fscanf(f,"%*s %*s %*s\n");
if(strcmp(type_of_code, "rustico") == 0){fscanf(f,"%*s %*s %*s %s\n",name_data_in);}
if(strcmp(type_of_code, "rusticoX") == 0){fscanf(f,"%*s %*s %*s %s %s\n",name_data_in,name_dataB_in);}

if(strcmp(type_of_code, "rustico") == 0){fscanf(f,"%*s %*s %*s %s\n",name_randoms_in);}
if(strcmp(type_of_code, "rusticoX") == 0){fscanf(f,"%*s %*s %*s %s %s\n",name_randoms_in,name_randomsB_in);}


fscanf(f,"%*s %*s %*s %s\n",name_path_out);
printf("Output files at %s\n",name_path_out);
fscanf(f,"%*s %*s %*s %s\n",name_id);
printf("Output Id %s\n",name_id);
fscanf(f,"%*s %*s %*s %s\n",header);
printf("Write header? %s\n",header);
fscanf(f,"%*s %*s %*s %s\n",output_density);
printf("Write density? %s\n\n",output_density);
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %d\n",&power_grid);
ngrid=pow(2,power_grid);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",type_of_mass_assigment);
printf("== FFT options ==\n");
printf("Number of k-modes and grid cells per side: %d\n",ngrid);
printf("Type of mass assingment %s\n",type_of_mass_assigment);

if(strcmp(type_of_mass_assigment, "NGC") == 0){mode_correction=1;}
if(strcmp(type_of_mass_assigment, "CIC") == 0){mode_correction=2;}
if(strcmp(type_of_mass_assigment, "TSC") == 0){mode_correction=3;}
if(strcmp(type_of_mass_assigment, "PCS") == 0){mode_correction=4;}
if(strcmp(type_of_mass_assigment, "P4S") == 0){mode_correction=5;}
if(strcmp(type_of_mass_assigment, "P5S") == 0){mode_correction=6;}
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_yamamoto);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Type of Yamamoto: %s\n",type_of_yamamoto);}
fscanf(f,"%*s %*s %*s %*s %*s %d\n",&Ninterlacing);
printf("Number of Interlacing steps %d\n",Ninterlacing);
fscanf(f,"%*s %*s %*s %*s %s\n",grid_correction_string);
printf("Do grid correction? %s\n\n",grid_correction_string);
if(strcmp(grid_correction_string, "no") == 0){mode_correction=0;}
fscanf(f,"%*s %*s\n");
    if(strcmp(type_of_code, "rustico") == 0){fscanf(f,"%*s %*s %*s %lf %lf\n",&z_min,&z_max);}
    if(strcmp(type_of_code, "rusticoX") == 0){fscanf(f,"%*s %*s %*s %lf %lf %lf %lf\n",&z_min,&z_max,&z_minB,&z_maxB);}

if(strcmp(type_of_survey, "cutsky") == 0){printf("== Cutsky options ==\n");}
if(strcmp(type_of_survey, "cutsky") == 0){printf("Redshift cuts: %lf < z < %lf\n",z_min,z_max);}
if(strcmp(type_of_code, "rusticoX") == 0 && strcmp(type_of_survey, "cutsky") == 0){printf("Redshift cuts: %lf < z < %lf\n",z_minB,z_maxB);}
fscanf(f,"%*s %*s %*s %*s %lf %lf\n",&Omega_m,&Omega_L);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Omega_m: %lf; Omega_L: %lf; Omega_k=%lf\n",Omega_m,Omega_L,1-Omega_m-Omega_L);}
if(strcmp(type_of_code, "rustico") == 0){fscanf(f,"%*s %*s %*s %*s %*s %*s %lf\n",&Area_survey);}
if(strcmp(type_of_code, "rusticoX") == 0){fscanf(f,"%*s %*s %*s %*s %*s %*s %lf %lf\n",&Area_survey,&Area_surveyB);}
if(strcmp(type_of_survey, "cutsky") == 0){printf("Area of the survey %lf deg2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_code, "rusticoX") == 0){printf("Area of the survey %lf deg2\n",Area_surveyB);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %lf\n",&Rsmoothing);
fscanf(f,"%*s %*s %*s %s\n",Quadrupole_type);
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_anisotropy, "yes") == 0){printf("Do Quadrupole as %s\n",Quadrupole_type);}
fscanf(f,"%*s %*s %*s %s\n",Octopole_type);
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_odd_multipoles, "yes") == 0 && strcmp(do_anisotropy, "yes") == 0){printf("Do Octopole as %s\n",Octopole_type);}
fscanf(f,"%*s %*s %*s %s\n",Hexadecapole_type);
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_anisotropy, "yes") == 0){printf("Do Hexadecapole as %s\n",Hexadecapole_type);}
fscanf(f,"%*s %*s %*s %*s %s\n",type_normalization_mode);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Compute normalization using %s\n",type_normalization_mode);}
fscanf(f,"%*s %*s %*s %*s %s\n",type_normalization_mode2);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Compute normalization using %s file n(z)\n",type_normalization_mode2);}
fscanf(f,"%*s %*s %*s %*s %*s %lf\n",&Shot_noise_factor);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Shot noise factor set to %lf\n",Shot_noise_factor);}
fscanf(f,"%*s %*s %*s %s\n",shuffle);
printf("\n == Shuffling Options ==\n");
if( strcmp(shuffle, "no") ==0 ){printf("Random catalogue provided used.\n");}
if( strcmp(shuffle, "radec") ==0 ){printf("RA-dec shuffling.\n");}
if( strcmp(shuffle, "redshift") ==0 ){printf("Redshift shuffling\n");}
if( strcmp(shuffle, "both") ==0 ){printf("Both RA-dec & redshift shuffling (random catalogue only used for setting the size of the new random catalogue).\n");}

fscanf(f,"%*s %*s %*s %*s %s\n",write_shuffled_randoms);
if( strcmp(write_shuffled_randoms, "yes") == 0 ||  strcmp(write_shuffled_randoms, "no") ==0 ){printf("Write Shuffled randoms: %s\n",write_shuffled_randoms);}

    fscanf(f,"%*s %*s %*s\n");
    fscanf(f,"%*s %*s %*s %*s %*s %s\n",window_function);
    if( strcmp(window_function, "yes") == 0 ||  strcmp(window_function, "no") ==0 ){printf("Compute window function: %s\n",window_function);}
    fscanf(f,"%*s %*s %*s %*s %*s %d\n",&window_norm_bin);
    fscanf(f,"%*s %*s %*s %lf\n",&deltaS_window);
    fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %lf\n",&percentage_randoms_window);
    fscanf(f,"%*s %*s %*s %s\n",yamamoto4window);

    if( strcmp(window_function, "yes") == 0)
    {
    //print window conditions
    printf("Window normalized to 1 on the bin %d\n",window_norm_bin);
    printf("Window binning %lf Mpc/h\n",deltaS_window);
    printf("Percentage of randoms for window calculation %lf %\n",percentage_randoms_window);
    printf("Yamamoto approximation for the window: %s\n\n",yamamoto4window);
    }

    fscanf(f,"%*s %*s %*s\n");
    fscanf(f,"%*s %*s %*s %*s %lf\n",&Pnoisedelta);
    fscanf(f,"%*s %*s %*s %*s %lf\n",&Bnoise1delta);
    fscanf(f,"%*s %*s %*s %*s %lf\n",&Bnoise2delta);
    fscanf(f,"%*s %*s %*s %lf\n",&I22delta);
    fscanf(f,"%*s %*s %*s %lf\n",&I33delta);
 
   if( strcmp(type_of_input, "density") == 0)
   {
     printf("Shot noise input: %lf\n",Pnoisedelta);
     printf("Normalization input: %lf\n",I22delta);
     if( strcmp(do_bispectrum,"yes") == 0){
     printf("Shot noise input for bispectrum: %lf, %lf\n",Bnoise1delta,Bnoise2delta);
     printf("Normalization input for bispectrum: %lf\n",I33delta);
      }

   }
fclose(f);
    

if(  strcmp(type_of_input,"particles") ==0 ){

  printf("\n== Read in/out options ==\n");
if(strcmp(type_of_file, "gadget") == 0)//periodic + gadget
{
  Ndata=0;
  for(snapshot_num=0;snapshot_num<gadget_files;snapshot_num++)
  {
     sprintf(name_gadget_file, "%s.%d", name_data_in, snapshot_num);
     g=fopen(name_gadget_file,"r");
     if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_data_in);return 0;}
     fclose(g);
     Ndata+=count_particles_gadget(name_gadget_file);
  }
P_shot_noise=pow(L2-L1,3)/Ndata*1.;
printf("Reading data file %s; %ld lines\n",name_data_in,Ndata);

    if(strcmp(type_of_code, "rusticoX") == 0)
    {
       if(strcmp(type_of_fileB, "gadget") == 0)//periodic + gadget
       {
      NdataB=0;
      for(snapshot_num=0;snapshot_num<gadget_filesB;snapshot_num++)
      {
         sprintf(name_gadget_fileB, "%s.%d", name_dataB_in, snapshot_num);
         g=fopen(name_gadget_fileB,"r");
         if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_dataB_in);return 0;}
         fclose(g);
         NdataB+=count_particles_gadget(name_gadget_fileB);
      }
    P_shot_noiseB=pow(L2-L1,3)/NdataB*1.;
    printf("Reading data file %s; %ld lines\n",name_dataB_in,NdataB);
    }
       if(strcmp(type_of_fileB, "ascii") == 0)//periodic + gadget
       {

           g=fopen(name_dataB_in,"r");
           if(g==NULL &&  strcmp(window_function, "no") == 0 ){printf("File %s does not exist. Exiting now...\n",name_dataB_in);return 0;}
           if(g!=NULL){
           fclose(g);
           NdataB=countlines(name_dataB_in);
  }
  if(g==NULL){Ndata=0;}//you can continue but just to get the window from the randoms

  if(NdataB>0){
  printf("Reading data file %s; %ld lines\n",name_dataB_in,NdataB);
  }


       }
   }
}

  if(strcmp(type_of_file, "ascii") == 0)//periodic or cutsky
  {

  g=fopen(name_data_in,"r");
  if(g==NULL &&  strcmp(window_function, "no") == 0 ){printf("File %s does not exist. Exiting now...\n",name_data_in);return 0;}
  if(g!=NULL){
  fclose(g);
  Ndata=countlines(name_data_in);
  }
  if(g==NULL){Ndata=0;}//you can continue but just to get the window from the randoms

  if(Ndata>0){
  printf("Reading data file %s; %ld lines\n",name_data_in,Ndata);
  }

    if(strcmp(type_of_code, "rusticoX") == 0)
    {
       if(strcmp(type_of_fileB, "ascii") == 0)
       {
      g=fopen(name_dataB_in,"r");
      if(g==NULL &&  strcmp(window_function, "no") == 0){printf("File %s does not exist. Exiting now...\n",name_dataB_in);return 0;}
      if(g!=NULL){
      fclose(g);
      NdataB=countlines(name_dataB_in);
    }
    if(g==NULL){NdataB=0;}//you can continue but just to get the window from the randoms

      }
       if(strcmp(type_of_fileB, "gadget") == 0)
       { 
    NdataB=0;
      for(snapshot_num=0;snapshot_num<gadget_filesB;snapshot_num++)
      {
         sprintf(name_gadget_fileB, "%s.%d", name_dataB_in, snapshot_num);
         g=fopen(name_gadget_fileB,"r");
         if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_dataB_in);return 0;}
         fclose(g);
         NdataB+=count_particles_gadget(name_gadget_fileB);
      }
    P_shot_noiseB=pow(L2-L1,3)/NdataB*1.;
    printf("Reading data file %s; %ld lines\n",name_dataB_in,NdataB);


       }

    }

  if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0)
  {
  g=fopen(name_randoms_in,"r");
  if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_randoms_in);return 0;}
  fclose(g);

  Nrand=countlines(name_randoms_in);
  printf("Reading randoms file %s; %ld lines\n",name_randoms_in,Nrand);

    if(strcmp(type_of_code, "rusticoX") == 0){
       g=fopen(name_randomsB_in,"r");
       if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_randomsB_in);return 0;}
       fclose(g);

  NrandB=countlines(name_randomsB_in);
  printf("Reading randoms file %s; %ld lines\n",name_randomsB_in,NrandB);

       }
  }

}

}
//exit(0);
if(  strcmp(type_of_input,"density") ==0 ){
    
     g=fopen(name_data_in,"r");
     if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_data_in);return 0;}
     fclose(g);
           Ndata=countlines(name_data_in);
           if(ngrid*ngrid*ngrid != Ndata){printf("Error, input density vector has different size than expected by the size of cells inputed: Ndata=%ld, Ngrid=%d**3.Exiting now...\n",Ndata,ngrid);exit(0);}

/*           NdataB=countlines(name_dataB_in);
           if(ngrid*ngrid*ngrid != NdataB){printf("Error, input density vector has different size than expected by the size of cells inputed: NdataB=%ld, Ngrid=%d**3.Exiting now...\n",NdataB,ngrid);exit(0);}
*/


}


//Error conditions

if(  strcmp(type_of_input,"particles") ==0 || strcmp(type_of_input,"density") ==0 ){

if(strcmp(type_of_file, "gadget") != 0 && strcmp(type_of_file, "ascii") != 0){printf("File type must be either 'gadget' or 'ascii'. Entry read %s. Exiting now...\n",type_of_file);return 0;}
if(strcmp(type_of_code, "rusticoX") == 0){if(strcmp(type_of_fileB, "gadget") != 0 && strcmp(type_of_fileB, "ascii") != 0){printf("File type must be either 'gadget' or 'ascii'. Entry read %s. Exiting now...\n",type_of_fileB);return 0;}}

if( strcmp(type_of_survey, "cutsky") != 0 && strcmp(type_of_survey, "periodic") != 0  && strcmp(type_of_survey, "periodicFKP") != 0){printf("Survey type must be either 'cutsky', 'periodic' or 'periodicFKP'. Entry read %s. Exiting now...\n",type_of_survey);return 0;}
if( strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_file, "gadget") == 0 ){printf("Warning. Cutsky+gadget option not available. Exiting now...\n");return 0;}
if( strcmp(type_of_survey, "periodicFKP") == 0 && strcmp(type_of_file, "gadget") == 0 ){printf("Warning. PeriodicFKP + gadget option not available. Exiting now...\n");return 0;}
if(strcmp(type_of_code, "rusticoX") == 0){
if( strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_fileB, "gadget") == 0 ){printf("Warning. Cutsky+gadget option not available. Exiting now...\n");return 0;}
if( strcmp(type_of_survey, "periodicFKP") == 0 && strcmp(type_of_fileB, "gadget") == 0 ){printf("Warning. PeriodicFKP + gadget option not available. Exiting now...\n");return 0;}
}
if(gadget_files<1 && strcmp(type_of_file,"gadget") == 0){printf("Warning. gadget files entry must be >0. Entered value %d. Exiting now...\n",gadget_files);return 0;}
if(strcmp(type_of_code, "rusticoX") == 0){if(gadget_filesB<1 &&  strcmp(type_of_fileB,"gadget") == 0){printf("Warning. gadget files entry must be >0. Entered value %d. Exiting now...\n",gadget_files);return 0;}}
if( strcmp(RSD, "yes") != 0 && strcmp(RSD, "no") != 0 ){printf("Warning. RSD option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",RSD);return 0;}
if(strcmp(type_of_code, "rusticoX") == 0){if( strcmp(RSDB, "yes") != 0 && strcmp(RSDB, "no") != 0 ){printf("Warning. RSD option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",RSDB);return 0;}}
if(L2<=L1){printf("Error: L2 parameter has to be large then L1. Exiting now...\n");return 0;}
if( strcmp(type_of_computation, "DSE") !=0 &&  strcmp(type_of_computation, "DSY") !=0 &&  strcmp(type_of_computation, "FFT") !=0){printf("Type of computation only accepts 'DSE', 'DSY' or 'FFT' options. Entry read %s. Exiting now...\n",type_of_computation);return 0;}
if( strcmp(type_of_computation, "FFT") !=0 && strcmp(do_bispectrum, "yes") ==0){printf("Warning. Bispectrum computation only accepts FFT computation option. Entry read %s. Exiting now...\n",do_bispectrum);return 0;}
if( strcmp(binning_type, "log10") !=0 && strcmp(binning_type, "linear") !=0){printf("Warning. Binning type must be either 'log10' or 'linear'. Entry read %s. Exiting now...\n",binning_type);return 0;}
if( strcmp(binning_type, "log10") ==0 && strcmp(do_bispectrum, "yes") ==0){printf("Warning. 'log10' type of binning not available for bispectrum computation on this version. Exiting now...\n"); return 0;}
if(bin_ps<=0){printf("Error: Size of the power spectrum bin has to be greater than 0: Entry read %lf. Exiting now...\n",bin_ps);return 0;}
if(kmin<0 || kmax<=0 || kmin>kmax){printf("Warning: Unusual values for maximum and/or minimum k-values for the bispectrum computation: kmin=%lf, kmax=%lf. Exiting now...\n",kmin, kmax);return 0;}
if(kmin==0 && strcmp(binning_type, "log10") == 0){printf("Cannot set kmin=0 and log-k binning. Exiting now...\n");return 0;}
if( strcmp(do_bispectrum, "yes") != 0 && strcmp(do_bispectrum, "no") != 0){printf("Error. Bispectrum entry must be either 'yes' or 'no'. Read entry %s. Exiting now...\n",do_bispectrum);return 0;}
if( strcmp(bispectrum_optimization,"lowmem") != 0 && strcmp(bispectrum_optimization,"standard") != 0 && strcmp(bispectrum_optimization,"himem") != 0){printf("Error. Bispectrum Optimmization must be either 'lowmem', 'standard' or 'himem'. Read entry %s. Exiting now...\n", bispectrum_optimization);exit(0);}
if( strcmp(do_multigrid, "yes") != 0 && strcmp(do_multigrid, "no") != 0){printf("Error. Multigrid entry must be either 'yes' or 'no'. Read entry %s. Exiting now...\n",do_multigrid);return 0;}
if( strcmp(do_multigrid, "yes") == 0 && Ninterlacing<2 ){printf("Warning. Multigrid option requires a number of interlacing steps >1. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"NGC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"CIC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
//if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"TSC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(triangle_shapes,"ALL") !=0 && strcmp(triangle_shapes,"EQU") !=0 && strcmp(triangle_shapes,"ISO") !=0 && strcmp(triangle_shapes,"SQU") !=0 ){printf("Error. Triangle shapes entry only accepts 'ALL', 'ISO', 'EQU' or 'SQU'. Read entry %s. Exiting now...\n",triangle_shapes);return 0;}
if(Deltakbis<=0){printf("Error: Size of the bispectrum bin has to be greater than 0: %lf kf. Exiting now...\n",Deltakbis);return 0;}
if( strcmp(triangles_num, "FFT") != 0 && strcmp(triangles_num, "APR_SUM") !=0 && strcmp(triangles_num, "EXA_SUM") !=0  && strcmp(triangles_num, "APR_EFF") !=0 && strcmp(triangles_num, "EXA_EFF") !=0){printf("Error. Number of Triangles normalization option only accepts: 'FFT', 'APR_SUM', 'EXA_SUM', 'APR_EFF', 'EXA_EFF'. Entry read %s. Exiting now...\n",triangles_num);return 0;}
if( strcmp(triangles_num, "EXA_SUM") == 0 &&  strcmp(triangle_shapes,"EQU") !=0){printf("Warning. 'EXA_SUM' triangle normalization option it is only available for equilateral triangles. Exiting now...\n");return 0;}//not available at the moment
if( strcmp(triangles_num, "EXA_EFF") == 0 &&  strcmp(triangle_shapes,"EQU") !=0){printf("Warning. 'EXA_EFF' triangle normalization option it is only available for equilateral triangles. Exiting now...\n");return 0;}//not available at the moment
if( strcmp(write_triangles, "yes") != 0 &&  strcmp(write_triangles, "no") !=0){printf("Error. write triangle entry must be either 'yes' or 'now'. Exiting now...\n");return 0;}
if( strcmp(write_triangles, "yes") == 0 && strcmp(triangle_shapes,"SQU") !=0 && strcmp(type_of_code,"rusticoX") == 0){printf("Warning. Write triangle option 'yes' is only recomended for squeezed triangle shapes 'SQU' and in auto-bispectrum mode. Exiting now...\n");return 0;}
if( strcmp(header, "yes") !=0 && strcmp(header, "no") !=0){printf("Warning. Write header option must be either 'yes' or 'no'. Entry read %s. Exiting now...\n",header);return 0;}

if( strcmp(output_density, "corrected") !=0 && strcmp(output_density, "uncorrected") !=0 && strcmp(output_density, "no") !=0){printf("Warning. Write density option must be either 'corrected', 'uncorrected' or 'no'. Entry read %s. Exiting now...\n",output_density);return 0;}

if(power_grid<4 || power_grid>15){printf("Warning: Unusual value for number of grid cells per side: 2^%d=%d. Exiting now...\n",power_grid,ngrid);return 0;}
if(strcmp(type_of_mass_assigment,"NGC") !=0 && strcmp(type_of_mass_assigment,"CIC") !=0 && strcmp(type_of_mass_assigment,"TSC") !=0 && strcmp(type_of_mass_assigment,"PCS") !=0 && strcmp(type_of_mass_assigment,"P4S") !=0 &&  strcmp(type_of_mass_assigment,"P5S") !=0){printf("Error. Type of mass assigment must be either 'NGC', 'CIC', 'TSC', 'PCS', 'P4S' or 'P5S'. Entry read %s. Exiting now...\n",type_of_mass_assigment);return 0;}
if( strcmp(type_of_yamamoto, "GridCenter") != 0 && strcmp(type_of_yamamoto, "GridAverage") != 0){printf("Error. Type of Yamamoto option must be either 'GridCenter' or 'GridAverage'. Entry read %s. Exiting now...\n",type_of_yamamoto);return 0;}
if(Ninterlacing<=0){printf("Error: Number of interglacing steps has to be equal or larger than 1. %d\n",Ninterlacing);return 0;}
if( strcmp(grid_correction_string, "yes") !=0 && strcmp(grid_correction_string, "no") !=0){printf("Grid correction input must be either 'yes' or 'no'. Entry read %s. Exiting now...\n",grid_correction_string);return 0;}

if( strcmp(type_of_computation, "FFT") !=0 && strcmp(output_density, "corrected") == 0){printf("Warning. In order to write out the densities the type of calculation must be FFT. Exiting now...\n");exit(0);}
if( strcmp(type_of_computation, "FFT") !=0 && strcmp(output_density, "uncorrected") == 0){printf("Warning. In order to write out the densities the type of calculation must be FFT. Exiting now...\n");exit(0);}

if( Ninterlacing !=1 && strcmp(output_density, "uncorrected") == 0){printf("Warning. In order to write out the densities (uncorrected) the number of interlacing steps must be 1. Exiting now...\n");exit(0);}


//

if(strcmp(type_of_survey, "cutsky") == 0){

if(Rsmoothing<0 || Rsmoothing>40){printf("Warning. Strange value for Rsmoothing: %lf Mpc/h. Exiting now...\n",Rsmoothing);exit(0);}

if(z_min>=z_max){printf("Error. Minimum value for redshift is larger than the maximum: z_min=%lf; z_max=%lf. Exiting now...\n",z_min,z_max);return 0;}
if(strcmp(type_of_code, "rusticoX") == 0){if(z_minB>=z_maxB){printf("Error. Minimum value for redshift is larger than the maximum: z_min=%lf; z_max=%lf. Exiting now...\n",z_minB,z_maxB);return 0;}}
if(Omega_m<=0 || Omega_m>1){printf("Warning. Unusual value for Omega_m, Omega_m=%lf. Exiting now...\n",Omega_m);return 0;}
if(Omega_L<=0 || Omega_L>1){printf("Warning. Unusual value for Omega_L, Omega_L=%lf. Exiting now...\n",Omega_L);return 0;}
if(Area_survey<=0){printf("Warning. Usual value for the Area of the survey: %lf. Exiting now...\n",Area_survey);return 0;}
if(strcmp(type_of_code, "rusticoX") == 0){if(Area_surveyB<=0){printf("Warning. Usual value for the Area of the survey: %lf. Exiting now...\n",Area_surveyB);return 0;}}
if( strcmp(Hexadecapole_type, "L0L4") !=0 && strcmp(Hexadecapole_type, "L2L2") !=0 && strcmp(Hexadecapole_type, "L1L3") !=0){printf("Hexadecapole option must be either 'L2L2' or 'L0L4' or 'L1L3'. Entry read %s. Exiting now...\n",Hexadecapole_type);return 0;}
if( strcmp(Octopole_type, "L0L3") !=0 && strcmp(Octopole_type, "L1L2") !=0){printf("Octopole option must be either 'L1L2' or 'L0L3'. Entry read %s. Exiting now...\n",Octopole_type);return 0;}
if( strcmp(Quadrupole_type, "L0L2") !=0 && strcmp(Quadrupole_type, "L1L1") !=0){printf("Quadrupole option must be either 'L1L1' or 'L0L2'. Entry read %s. Exiting now...\n",Quadrupole_type);return 0;}

if( strcmp(type_normalization_mode, "area") !=0 && strcmp(type_normalization_mode, "density") !=0){printf("Error. Normalisation type must be either 'area' or 'density'. Entry read %s. Exiting now...\n",type_normalization_mode);return 0;}
if( strcmp(type_normalization_mode2, "data") !=0 && strcmp(type_normalization_mode2, "randoms") !=0){printf("Error. Normalisation type must be either 'data' or 'randoms'. Entry read %s. Exiting now...\n",type_normalization_mode2);return 0;}
if(Shot_noise_factor>1 || Shot_noise_factor<0){printf("Warning. Usual value for the Shot noise factor: %lf. Exiting now...\n",Shot_noise_factor);return 0;}
if( strcmp(shuffle, "radec") != 0 && strcmp(shuffle, "no") != 0 && strcmp(shuffle, "redshift") != 0 && strcmp(shuffle, "both") != 0){printf("Warning. Shuffling option not reconised: %s. Exiting now...\n",shuffle);return 0;}
if( strcmp(window_function, "yes") != 0 && strcmp(window_function, "no") != 0 ){printf("Warning. Window function option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",window_function);return 0;}
if( strcmp(write_shuffled_randoms, "yes") != 0 && strcmp(write_shuffled_randoms, "no") != 0 ){printf("Warning. Write shuffled randoms option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",write_shuffled_randoms);return 0;}

if( strcmp(shuffle, "no") != 0 && Ndata==0){printf("Error, when shuffling option enabled, a valid data file must be provided. Exiting now...\n");exit(0);}
if( strcmp(shuffle, "no") != 0 && NdataB==0 && strcmp(type_of_code,"rusticoX") == 0 ){printf("Error, when shuffling option enabled, a valid data file must be provided. Exiting now...\n");exit(0);}

}//cutsky 
 
    if( strcmp(window_function, "yes") == 0){

    if(window_norm_bin<=0 || window_norm_bin>20){printf("Error, the minimum bin at which the window is normalized must be >0 and <20: %d. Exiting now...\n",window_norm_bin);exit(0);}
    if(deltaS_window<=0){printf("Error, Delta-s for window function can't be negative: %lf. Exiting now...\n",deltaS_window);exit(0);}
    if(deltaS_window>100){printf("Warning, Delta-s for window function seems too large: %lf. Exiting now...\n",deltaS_window);exit(0);}
    if(percentage_randoms_window<=0 || percentage_randoms_window>100){printf("Error, wrong  percentage for the selection of randoms when computing the window: %lf %. Exiting now...\n",percentage_randoms_window);exit(0);}
    if( strcmp(yamamoto4window,"yes") !=0 &&  strcmp(yamamoto4window,"no") != 0){printf("Error, Yamamoto approximation for the window must be yes/no  option: %s. Exiting now...\n",yamamoto4window);exit(0);}

    if( strcmp(name_randoms_in,"none") == 0){printf("Error, no random file selected for the window calculation: %s. Exiting now...\n",name_randoms_in);exit(0);}
    g=fopen(name_randoms_in,"r");
    if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_randoms_in);return 0;}
    fclose(g);


    }


if( strcmp(do_bispectrum, "no") == 0 && strcmp(do_bispectrum2, "yes") == 0 ){printf("Warning. Bispectrum multipoles computation requires a bispectrum calculation. Exiting now...\n");return 0;}
    
    if( strcmp(do_bispectrum2, "yes") == 0 && strcmp(type_of_code, "rusticoX") == 0 && strcmp(do_bispectrum, "yes") == 0){printf("Warning. Bispectrum multipoles computation not available for cross-statistics. Exiting now...\n");return 0;}

if( strcmp(do_bispectrum2, "yes") == 0 && strcmp(type_of_yamamoto,"GridAverage") == 0){printf("Warning, Yamamoto GridAverage not currently available with the bispectrum quadrupole. Exiting now...\n");exit(0);}


if( strcmp(do_mu_bins, "yes") != 0 &&  strcmp(do_mu_bins, "no") != 0){printf("Warning. Mu-binning option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_mu_bins);return 0;}

if( strcmp(do_anisotropy, "yes") != 0 &&  strcmp(do_anisotropy, "no") != 0){printf("Warning. Anisotropy option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_odd_multipoles);return 0;}

if( strcmp(do_odd_multipoles, "yes") != 0 &&  strcmp(do_odd_multipoles, "no") != 0){printf("Warning. Odd multipole option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_odd_multipoles);return 0;}

if( strcmp(write_kvectors, "yes") != 0 &&  strcmp(write_kvectors, "no") != 0){printf("Warning. Write k-vectors option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",write_kvectors);return 0;}

if( strcmp(do_mu_bins, "yes") == 0 &&  strcmp(type_of_survey, "cutsky") == 0){printf("Warning. Mu-binning option only available for non-varying line-of-sight samples. Exiting now...\n");return 0;}

if(mubin<0){printf("Warning. The numer of mu-bins has to be a positive integer. Exiting now...\n");return 0;}

if( strcmp(do_mu_bins, "yes") == 0){
if( strcmp(file_for_mu, "yes") != 0 &&  strcmp(file_for_mu, "no") != 0){printf("Warning. Different files when Mu-binning option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",do_mu_bins);return 0;}
}

if( strcmp(do_mu_bins, "yes") == 0 &&  strcmp(do_odd_multipoles, "yes") == 0){printf("Warning. Cannot do odd multipoles if mu-binning option is also selected. Enable either mu-binning or odd-multipoles. Exiting now...\n");return 0;}

if( strcmp(do_anisotropy, "no") == 0 &&  strcmp(do_odd_multipoles, "yes") == 0){printf("Warning. Cannot do odd multipoles if anisotropy option is not selected. Enable either anisotropy or disable odd-multipoles. Exiting now...\n");return 0;}

if( strcmp(do_mu_bins, "yes") == 0 &&  strcmp(do_anisotropy, "no") == 0){printf("Warning. Cannot do mu-binning if anisotropy option is not selected. Enable either anisotropy or disable mu-binning. Exiting now...\n");return 0;}

if( strcmp(type_of_computation, "DSE") ==0 && strcmp(do_anisotropy, "no") == 0){printf("Warning. There's no point of doing an DSE computation for an isotropic signal, as this is equivalent to a DSY or FFT calculation. Exiting now...\n");exit(0);}

//if( strcmp(type_of_computation, "FFT") !=0 && strcmp(output_density, "corrected") == 0){printf("Warning. In order to write out the densities the type of calculationn must be FFT. Exiting now...\n");exit(0);}
//if( strcmp(type_of_computation, "FFT") !=0 && strcmp(output_density, "uncorrected") == 0){printf("Warning. In order to write out the densities the type of calculationn must be FFT. Exiting now...\n");exit(0);}

//if( Ninterlacing !=1 && strcmp(output_density, "uncorrected) == 0){printf("Warning. In order to write out the densities the number of interlacing steps must be 1. Exiting now...\n");exit(0);}

//}//old cutsky

}//particles
if(  strcmp(type_of_input,"density") ==0 ){

//rustico only
if( strcmp(type_of_code,"rusticox") == 0){printf("Error, density-input option not available for X-power (rusticoX). Exiting now...\n");exit(0);}

//Interlacing must be 1
if( Ninterlacing>1){printf("Error, can't use interlacing technique with a density-input option. Exting now...\n");exit(0);}

if( strcmp(type_of_computation,"FFT") !=0  ){printf("Error, density-input option requires FFT type of computation. Exting now...\n");exit(0);}


//I22,I33 >0
if(I22delta<=0 || I33delta<=0){printf("Error, input I22 and I33 normalizations must be positive: I22=%lf, I33=%lf. Exiting  now...\n",I22delta,I33delta);exit(0);}

//random file must be none
if( strcmp(name_randoms_in,"none") !=0){printf("Error, For the density-input option, randoms must be set to be none. Exting now...\n");exit(0);}
//Grid center option only

//can't ouput density if density is the input
if( strcmp(output_density, "corrected") == 0){printf("Warning, printing density option enabled when density is also selected as an input. Exiting now...\n");exit(0);}
if( strcmp(output_density, "uncorrected") == 0){printf("Warning, printing density option enabled when density is also selected as an input. Exiting now...\n");exit(0);}

if( strcmp(output_density, "corrected") == 0 && strcmp(grid_correction_string, "no") == 0){printf("Warning, can't produce a corrected-delta field if grid-correction option is off. Exiting now...\n");exit(0);}

if( strcmp( type_of_yamamoto ,"GridAverage") ==0 ){printf("Error, density-input option requires GridCenter Yamamoto. Extiting now...\n");exit(0);}

//shuffle set to no
if( strcmp( shuffle,"no") !=0){printf("Warning, no shuffle option available when density-input option is selected. Exiting now...\n");exit(0);}

//window computation set to no
if( strcmp( window_function, "yes") ==0){printf("Warning, no window computation option available when density-input option is selected. Exiting now...\n");exit(0);}
}

if( strcmp(type_of_survey,"periodicFKP") == 0 &&  strcmp( shuffle,"no") !=0  ){printf("Warning, no shuffle option for periodic boxes. Exiting now...\n");exit(0);}

//etc strings.....

//Determine number of processors available for openmp  
        #pragma omp parallel for private(i,tid) shared(n_lines_parallel,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){n_lines_parallel=omp_get_num_threads();}
        }
        i=fftw_init_threads();
        fftw_plan_with_nthreads(n_lines_parallel);
        printf("Number of processors used: %d\n",n_lines_parallel);


if( strcmp(output_density,"no") != 0)//either corrected or uncorrected
{

if(strcmp(output_density,"corrected") == 0){
if(strcmp(type_of_code, "rustico") == 0){sprintf(name_out_density,"%s/delta_%s_%s.txt",name_path_out,type_of_mass_assigment,name_id);}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_out_density,"%s/deltaA_%s_%s.txt",name_path_out,type_of_mass_assigment,name_id);}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_out_densityB,"%s/deltaB_%s_%s.txt",name_path_out,type_of_mass_assigment,name_id);}
}
if(strcmp(output_density,"uncorrected") == 0){
if(strcmp(type_of_code, "rustico") == 0){sprintf(name_out_density,"%s/deltaraw_%s_%s.txt",name_path_out,type_of_mass_assigment,name_id);}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_out_density,"%s/deltarawA_%s_%s.txt",name_path_out,type_of_mass_assigment,name_id);}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_out_densityB,"%s/deltarawB_%s_%s.txt",name_path_out,type_of_mass_assigment,name_id);}
}

}
else
{
if(strcmp(type_of_code, "rustico") == 0){sprintf(name_out_density,"no");}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_out_density,"no");}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_out_densityB,"no");}
}


//Reading files.
double parameter_value[31];
double parameter_valueB[31];

if( strcmp(type_of_input,"particles") == 0){

        parameter_value[0]=Omega_m;
        parameter_value[29]=Omega_L;
        parameter_value[30]=speed_of_light;
        parameter_value[1]=z_min;
        parameter_value[2]=z_max;
        parameter_value[3]=Ndata;
        parameter_value[13]=Area_survey;

    if(strcmp(type_of_code, "rusticoX") == 0){
        parameter_valueB[0]=Omega_m;
        parameter_valueB[29]=Omega_L;
        parameter_valueB[30]=speed_of_light;
         parameter_valueB[1]=z_minB;
         parameter_valueB[2]=z_maxB;
         parameter_valueB[3]=NdataB;  
         parameter_valueB[13]=Area_surveyB;
    }
        

Ndata2=0;
if(strcmp(type_of_survey, "cutsky") == 0 && Ndata>0)
{
Ndata2=get_number_used_lines_data(name_data_in,parameter_value);if(Ndata2<1000){printf("Warning, unusual low value for Ndata=%ld\n",Ndata2);if(Ndata2==0){exit(0);}}
alpha_data1=parameter_value[12];

Ndata2B=0;
if(strcmp(type_of_code, "rusticoX") == 0 && NdataB>0){
Ndata2B=get_number_used_lines_data(name_dataB_in,parameter_valueB);if(Ndata2B<1000){printf("Warning, unusual low value for NdataB=%ld\n",Ndata2B);if(Ndata2B==0){exit(0);}}
alpha_data1B=parameter_valueB[12];
}
    
    
}
    
if(strcmp(type_of_survey, "periodic") == 0 || strcmp(type_of_survey, "periodicFKP") == 0)
{

parameter_value[3]=Ndata;
if( strcmp(type_of_file,"ascii") == 0){
Ndata2=get_number_used_lines_periodic(name_data_in, parameter_value,0);
Ndata2w=get_number_used_lines_weighted_periodic(name_data_in, parameter_value,0);}
else{//gadget (no-weight by default)
Ndata2=Ndata;
Ndata2w=Ndata;
}

if(Ndata2==0 || Ndata2w == 0){printf("Error, no particles to be used found in %s. Ndata=%ld, Ndata_w=%lf. Exiting now...\n",name_data_in,Ndata2,Ndata2w);exit(0);}

if(strcmp(type_of_survey, "periodicFKP") == 0)
{
parameter_value[3]=Nrand;
Nrand2=get_number_used_lines_periodic(name_randoms_in, parameter_value,1);
}
if(strcmp(type_of_code, "rusticoX") == 0){

parameter_valueB[3]=NdataB;
if( strcmp(type_of_file,"ascii") == 0){
Ndata2B=get_number_used_lines_periodic(name_dataB_in, parameter_valueB,0);
Ndata2Bw=get_number_used_lines_weighted_periodic(name_dataB_in, parameter_valueB,0);}
else{//gadget (no-weight by default)
Ndata2B=NdataB;
Ndata2Bw=NdataB;
}

if(Ndata2B==0 || Ndata2Bw == 0){printf("Error, no particles to be used found in %s. Ndata=%ld, Ndata_w=%lf. Exiting now...\n",name_dataB_in,Ndata2B,Ndata2Bw);exit(0);}
if(strcmp(type_of_survey, "periodicFKP") == 0)
{
parameter_valueB[3]=NrandB;
Nrand2B=get_number_used_lines_periodic(name_randomsB_in, parameter_valueB,1);
}

}
}

if(strcmp(type_of_file, "ascii") == 0 && Ndata2>0)//these are only kept stored during all the process for ascii files. Gadget files keep the name of the file and read it each time they need
{
        pos_x = (double*) calloc(Ndata2, sizeof(double));
        pos_y = (double*) calloc(Ndata2, sizeof(double));
        pos_z = (double*) calloc(Ndata2, sizeof(double));
	weight = (double*) calloc(Ndata2, sizeof(double));
}
    
if(strcmp(type_of_fileB, "ascii") == 0 && Ndata2B>0)//these are only kept stored during all the process for ascii files. Gadget files keep the name of the file and read it each time they need
{
            pos_xB = (double*) calloc(Ndata2B, sizeof(double));
            pos_yB = (double*) calloc(Ndata2B, sizeof(double));
            pos_zB = (double*) calloc(Ndata2B, sizeof(double));
            weightB = (double*) calloc(Ndata2B, sizeof(double));
}

if(strcmp(type_of_survey, "cutsky") == 0 && Ndata2>0)
{
//pos_x,pos_y,pos_z,weight are loaded with the position of particles
//Ndata is uploaded to the number of particles used
//Psn_1a,Psn_1b,Psn_2a,Psn_2b are uploaded with information on the shot noise
//z_efffective is uploaded with the effective redshift of the sample
//num_effective is uploaded with the effective number of particles
//I_norm is uploaded with information relative to the normalization based on density
//alpha is uploeaded with information relative to the effective ratio between data and randoms
printf("Reading %s...",name_data_in);

    I33=0;I22=0;IN=0;Bsn=0;alpha=0;

                radata = (double*) calloc(Ndata2, sizeof(double));
                decdata = (double*) calloc(Ndata2, sizeof(double));
                zdata = (double*) calloc(Ndata2, sizeof(double));
                wcoldata = (double*) calloc(Ndata2, sizeof(double));
                wsysdata = (double*) calloc(Ndata2, sizeof(double));
                wfkpdata = (double*) calloc(Ndata2, sizeof(double));
                nzdata = (double*) calloc(Ndata2, sizeof(double));
                
parameter_value[3]=Ndata;
get_skycuts_data(name_data_in, pos_x, pos_y, pos_z, weight, parameter_value,type_normalization_mode,radata,decdata,zdata,wcoldata,wsysdata,wfkpdata,nzdata,shuffle,Rsmoothing);


Ndata2=parameter_value[3];
Psn_1a=parameter_value[4];
Psn_1b=parameter_value[5];
I3_norm_data=parameter_value[6];
z_effective_data=parameter_value[7];
num_effective=parameter_value[8];
I_norm_data=parameter_value[9];
min=parameter_value[10];
max=parameter_value[11];
alpha_data=parameter_value[12];

I_norm_data2=parameter_value[14];
I_norm_data3=parameter_value[15];
I_norm_data4=parameter_value[16];
num_effective2=parameter_value[17];
num_effective3=parameter_value[28];

I3_norm_data2=parameter_value[18];
I3_norm_data3=parameter_value[19];
I3_norm_data4=parameter_value[20];

Bsn1a=parameter_value[21];
Bsn1b=parameter_value[22];
IN1=parameter_value[23];
IN2=parameter_value[24];
IN11=parameter_value[26];
IN22=parameter_value[27];

parameter_value[3]=Ndata;

if(strcmp(type_of_code, "rustico") == 0){sprintf(name_den_out,"%s/Density_galaxies_%s.txt",name_path_out,name_id); }
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_den_out,"%s/Density_galaxiesA_%s.txt",name_path_out,name_id); }
get_skycuts_write_density_data(name_data_in, parameter_value,name_den_out);

printf("Ok!\n");

parameter_value[0]=Omega_m;
parameter_value[29]=Omega_L;
parameter_value[30]=speed_of_light;
parameter_value[1]=z_min;
parameter_value[2]=z_max;
parameter_value[3]=Nrand;


if( strcmp(shuffle, "no") == 0 )
{
Nrand2=get_number_used_lines_randoms(name_randoms_in,parameter_value);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2);if(Nrand2==0){exit(0);}}
}
else//shuffle == yes
{
Nrand2=get_number_used_lines_randoms(name_randoms_in,parameter_value);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2);if(Nrand2==0){exit(0);}}
Nrand2=Nrand;//if shuffle is enabled number of used randoms is the total number of randoms, no matter of the z-cuts (those are applied later if necessary)
}
alpha_rand1=parameter_value[12];
alpha1=alpha_data1/alpha_rand1;


        pos_x_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_y_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_z_rand = (double*) calloc(Nrand2, sizeof(double));
        weight_rand = (double*) calloc(Nrand2, sizeof(double));

								
printf("Reading %s...",name_randoms_in);
parameter_value[12]=alpha_data;
if(strcmp(type_of_code, "rustico") == 0){sprintf(name_den_out,"%s/Density_randoms_%s.txt",name_path_out,name_id);}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_den_out,"%s/Density_randomsA_%s.txt",name_path_out,name_id);}

if(strcmp(type_of_code, "rustico") == 0){sprintf(name_out_randoms,"%s/Randoms_%s.txt",name_path_out,name_id);}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_out_randoms,"%s/RandomsA_%s.txt",name_path_out,name_id);}

get_skycuts_randoms(name_path_out,name_id,name_randoms_in, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, parameter_value,type_normalization_mode, type_normalization_mode2,radata,decdata,zdata,wcoldata,wsysdata,wfkpdata,nzdata,Ndata2,alpha1,shuffle,write_shuffled_randoms,name_out_randoms,Rsmoothing);

Nrand2=parameter_value[3];//now Nrand2 is the number of randoms used (for shuffle == no this is always the case, but for shuffle == yes not necesseraly) 
Psn_2a=parameter_value[4];
Psn_2b=parameter_value[5];
z_effective_rand=parameter_value[7];
num_effective_rand=parameter_value[8];
num_effective2_rand=parameter_value[20];
num_effective3_rand=parameter_value[28];

if(min>parameter_value[10]){min=parameter_value[10];}
if(max<parameter_value[11]){max=parameter_value[11];}
alpha_rand=parameter_value[12];

I_norm_rand=parameter_value[9];
I_norm_rand2=parameter_value[14];
I_norm_rand3=parameter_value[15];
I_norm_rand4=parameter_value[16];

I3_norm_rand=parameter_value[6];
I3_norm_rand2=parameter_value[17];
I3_norm_rand3=parameter_value[18];
I3_norm_rand4=parameter_value[19];
    
    Bsn2a=parameter_value[21];
    Bsn2b=parameter_value[22];

alpha=alpha_data/alpha_rand;
I_norm_rand=I_norm_rand*alpha;
I3_norm_rand=I3_norm_rand*alpha;

    if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "density") == 0 ){I22=I_norm_rand;I33=I3_norm_rand;}
    if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "density") == 0){I22=I_norm_data;I33=I3_norm_data;}

if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "area") == 0 ){I22=I_norm_rand2;I33=I3_norm_rand2;}
if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "area") == 0){I22=I_norm_data2;I33=I3_norm_data2;}


P_shot_noise1=(Psn_1a+alpha*alpha*Psn_2a)/I22;
P_shot_noise2=(Psn_1b+alpha*alpha*Psn_2b)/I22;
P_shot_noise=P_shot_noise1*Shot_noise_factor+P_shot_noise2*(1.-Shot_noise_factor);

Bsn1=(Bsn1a-alpha*alpha*alpha*Bsn2a)/I33;
Bsn2=(Bsn1b-alpha*alpha*alpha*Bsn2b)/I33;
    
Bsn=Bsn1*Shot_noise_factor+(1.-Shot_noise_factor)*Bsn2;
if( strcmp(type_normalization_mode, "area") == 0 ){IN=(IN1*Shot_noise_factor+IN2*(1.-Shot_noise_factor))/I33;}
if( strcmp(type_normalization_mode, "density") == 0 ){IN=(IN11*Shot_noise_factor+IN22*(1.-Shot_noise_factor))/I33;}

parameter_value[3]=Nrand;

if(strcmp(type_of_code, "rustico") == 0){sprintf(name_den_out,"%s/Density_randoms_%s.txt",name_path_out,name_id);}
if(strcmp(type_of_code, "rusticoX") == 0){sprintf(name_den_out,"%s/Density_randomsA_%s.txt",name_path_out,name_id);}
get_skycuts_write_density_randoms(name_randoms_in, parameter_value,alpha,name_den_out,radata,decdata,zdata,wcoldata,wsysdata,wfkpdata,nzdata,Ndata2,alpha1,shuffle);

free(radata);
free(decdata);
free(zdata);
free(wcoldata);
free(wsysdata);
free(wfkpdata);
free(nzdata);
printf("Ok!\n\n");

    I33B=0;I22B=0;INB=0;BsnB=0;alphaB=0;
    if(strcmp(type_of_code, "rusticoX") == 0 && Ndata2B>0){
        
    printf("Reading %s...",name_dataB_in);

                    radataB = (double*) calloc(Ndata2B, sizeof(double));
                    decdataB = (double*) calloc(Ndata2B, sizeof(double));
                    zdataB = (double*) calloc(Ndata2B, sizeof(double));
                    wcoldataB = (double*) calloc(Ndata2B, sizeof(double));
                    wsysdataB = (double*) calloc(Ndata2B, sizeof(double));
                    wfkpdataB = (double*) calloc(Ndata2B, sizeof(double));
                    nzdataB = (double*) calloc(Ndata2B, sizeof(double));
                    
parameter_valueB[3]=NdataB;
    get_skycuts_data(name_dataB_in, pos_xB, pos_yB, pos_zB, weightB, parameter_valueB,type_normalization_mode,radataB,decdataB,zdataB,wcoldataB,wsysdataB,wfkpdataB,nzdataB,shuffle,Rsmoothing);

    Ndata2B=parameter_valueB[3];
    Psn_1aB=parameter_valueB[4];
    Psn_1bB=parameter_valueB[5];
    I3_norm_dataB=parameter_valueB[6];
    z_effective_dataB=parameter_valueB[7];
    num_effectiveB=parameter_valueB[8];
    I_norm_dataB=parameter_valueB[9];
    minB=parameter_valueB[10];
    maxB=parameter_valueB[11];
    alpha_dataB=parameter_valueB[12];

    I_norm_data2B=parameter_valueB[14];
    I_norm_data3B=parameter_valueB[15];
    I_norm_data4B=parameter_valueB[16];
    num_effective2B=parameter_valueB[17];
    num_effective3B=parameter_valueB[28];

    I3_norm_data2B=parameter_valueB[18];
    I3_norm_data3B=parameter_valueB[19];
    I3_norm_data4B=parameter_valueB[20];

    Bsn1aB=parameter_valueB[21];
    Bsn1bB=parameter_valueB[22];
    IN1B=parameter_valueB[23];
    IN2B=parameter_valueB[24];
    IN11B=parameter_valueB[26];
    IN22B=parameter_valueB[27];

    parameter_valueB[3]=NdataB;

    sprintf(name_den_out,"%s/Density_galaxiesB_%s.txt",name_path_out,name_id);
    get_skycuts_write_density_data(name_dataB_in, parameter_valueB,name_den_out);
    parameter_valueB[3]=Ndata2B;

    printf("Ok!\n");

    parameter_valueB[0]=Omega_m;
    parameter_valueB[29]=Omega_L;
    parameter_valueB[30]=speed_of_light;
    parameter_valueB[1]=z_minB;
    parameter_valueB[2]=z_maxB;
    parameter_valueB[3]=NrandB;


    if( strcmp(shuffle, "no") == 0 )
    {
    Nrand2B=get_number_used_lines_randoms(name_randomsB_in,parameter_valueB);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2B);if(Nrand2B==0){exit(0);}}
    }
    else
    {
    Nrand2B=get_number_used_lines_randoms(name_randomsB_in,parameter_valueB);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2B);if(Nrand2B==0){exit(0);}}
    Nrand2B=NrandB;//if shuffle is enabled number of used randoms is the total number of randoms, no matter of the z-cuts (those are applied later if necessary)
    }
    alpha_rand1B=parameter_valueB[12];
    alpha1B=alpha_data1B/alpha_rand1B;


            pos_x_randB = (double*) calloc(Nrand2B, sizeof(double));
            pos_y_randB = (double*) calloc(Nrand2B, sizeof(double));
            pos_z_randB = (double*) calloc(Nrand2B, sizeof(double));
            weight_randB = (double*) calloc(Nrand2B, sizeof(double));

                                    
    printf("Reading %s...",name_randomsB_in);
    parameter_valueB[12]=alpha_dataB;
    sprintf(name_den_out,"%s/Density_randomsB_%s.txt",name_path_out,name_id);
    sprintf(name_out_randoms,"%s/RandomsB_%s.txt",name_path_out,name_id);

    get_skycuts_randoms(name_path_out,name_id,name_randomsB_in, pos_x_randB, pos_y_randB, pos_z_randB, weight_randB, parameter_valueB,type_normalization_mode, type_normalization_mode2,radataB,decdataB,zdataB,wcoldataB,wsysdataB,wfkpdataB,nzdataB,Ndata2B,alpha1B,shuffle,write_shuffled_randoms,name_out_randoms,Rsmoothing);

    Nrand2B=parameter_valueB[3];//now Nrand2 is the number of randoms used (for shuffle == no this is always the case, but for shuffle == yes not necesseraly)
    Psn_2aB=parameter_valueB[4];
    Psn_2bB=parameter_valueB[5];
    z_effective_randB=parameter_valueB[7];
    num_effective_randB=parameter_valueB[8];
    num_effective2_randB=parameter_valueB[20];
    num_effective3_randB=parameter_valueB[28];

    if(minB>parameter_valueB[10]){minB=parameter_valueB[10];}
    if(maxB<parameter_valueB[11]){maxB=parameter_valueB[11];}
    alpha_randB=parameter_valueB[12];

    I_norm_randB=parameter_valueB[9];
    I_norm_rand2B=parameter_valueB[14];
    I_norm_rand3B=parameter_valueB[15];
    I_norm_rand4B=parameter_valueB[16];

    I3_norm_randB=parameter_valueB[6];
    I3_norm_rand2B=parameter_valueB[17];
    I3_norm_rand3B=parameter_valueB[18];
    I3_norm_rand4B=parameter_valueB[19];
        
        Bsn2aB=parameter_valueB[21];
        Bsn2bB=parameter_valueB[22];

    alphaB=alpha_dataB/alpha_randB;
    I_norm_randB=I_norm_randB*alphaB;
    I3_norm_randB=I3_norm_randB*alphaB;

        if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "density") == 0 ){I22B=I_norm_randB;I33B=I3_norm_randB;}
        if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "density") == 0){I22B=I_norm_dataB;I33B=I3_norm_dataB;}

    if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "area") == 0 ){I22B=I_norm_rand2B;I33B=I3_norm_rand2B;}
    if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "area") == 0){I22B=I_norm_data2B;I33B=I3_norm_data2B;}


    P_shot_noise1B=(Psn_1aB+alphaB*alphaB*Psn_2aB)/I22B;
    P_shot_noise2B=(Psn_1bB+alphaB*alphaB*Psn_2bB)/I22B;
    P_shot_noiseB=P_shot_noise1B*Shot_noise_factor+P_shot_noise2B*(1.-Shot_noise_factor);

    Bsn1B=(Bsn1aB-alphaB*alphaB*alphaB*Bsn2aB)/I33B;
    Bsn2B=(Bsn1bB-alphaB*alphaB*alphaB*Bsn2bB)/I33B;
        
    BsnB=Bsn1B*Shot_noise_factor+(1.-Shot_noise_factor)*Bsn2B;
    if( strcmp(type_normalization_mode, "area") == 0 ){INB=(IN1B*Shot_noise_factor+IN2B*(1.-Shot_noise_factor))/I33B;}
    if( strcmp(type_normalization_mode, "density") == 0 ){INB=(IN11B*Shot_noise_factor+IN22B*(1.-Shot_noise_factor))/I33B;}

    parameter_valueB[3]=NrandB;

   sprintf(name_den_out,"%s/Density_randomsB_%s.txt",name_path_out,name_id);
    get_skycuts_write_density_randoms(name_randomsB_in, parameter_valueB,alphaB,name_den_out,radataB,decdataB,zdataB,wcoldataB,wsysdataB,wfkpdataB,nzdataB,Ndata2B,alpha1B,shuffle);

    free(radataB);
    free(decdataB);
    free(zdataB);
    free(wcoldataB);
    free(wsysdataB);
    free(wfkpdataB);
    free(nzdataB);
    printf("Ok!\n\n");
        
        
    }


}


//Compute the widow
    if(strcmp(window_function, "yes") == 0){

        parameter_value[0]=Omega_m;
        parameter_value[29]=Omega_L;
        parameter_value[30]=speed_of_light;
        parameter_value[1]=z_min;
        parameter_value[2]=z_max;
        if(strcmp(shuffle, "no") == 0){parameter_value[3]=Nrand;}//total number of randoms
        if(strcmp(shuffle, "yes") == 0){parameter_value[3]=Nrand2;}//total number of randoms
        parameter_value[13]=Area_survey;

    if(strcmp(type_of_code, "rusticoX") == 0){
        parameter_valueB[0]=Omega_m;
        parameter_valueB[29]=Omega_L;
        parameter_valueB[30]=speed_of_light;
         parameter_valueB[1]=z_minB;
         parameter_valueB[2]=z_maxB;
         if(strcmp(shuffle, "no") == 0){parameter_valueB[3]=NrandB;}//total number of randoms
         if(strcmp(shuffle, "yes") == 0){parameter_valueB[3]=Nrand2B;}//total number of randoms
         parameter_valueB[13]=Area_surveyB;
    }
    
    if(strcmp(type_of_code, "rustico") == 0){sprintf(name_wink_out,"%s/Window_%s.txt",name_path_out,name_id);}
    if(strcmp(type_of_code, "rusticoX") == 0){
    sprintf(name_wink_out,"%s/WindowAA_%s.txt",name_path_out,name_id);
    sprintf(name_winkBB_out,"%s/WindowBB_%s.txt",name_path_out,name_id);
    sprintf(name_winkAB_out,"%s/WindowAB_%s.txt",name_path_out,name_id);
    }
        

    f=fopen(name_wink_out,"w");
    if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_wink_out);return 0;}
    fclose(f);
if(strcmp(type_of_code, "rusticoX") == 0){
    f=fopen(name_winkBB_out,"w");
    if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_wink_out);return 0;}
    fclose(f);

    f=fopen(name_winkAB_out,"w");
    if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_wink_out);return 0;}
    fclose(f);
}

    printf("Computing the RR counts using %d processors. This can take a while...\n",n_lines_parallel);

 window_mask_function_RRcount(name_wink_out,name_winkBB_out,name_winkAB_out,name_randoms_in,name_randomsB_in,window_norm_bin,deltaS_window,percentage_randoms_window,yamamoto4window,parameter_value,parameter_valueB,header,L2-L1,n_lines_parallel,type_of_code,shuffle,pos_x_rand,pos_y_rand,pos_z_rand,weight_rand,pos_x_randB,pos_y_randB,pos_z_randB,weight_randB,Quadrupole_type,Hexadecapole_type,type_of_survey);


        if(strcmp(type_of_code, "rustico") == 0){
    if(Ndata==0){printf("Window function counts completed. No data file provided for power spectrum computation. Exiting now...\n");return 0;}
    if(Ndata>0){printf("Window function counts completed.\n\n");}
        }
        if(strcmp(type_of_code, "rusticoX") == 0){
            if(Ndata==0 || NdataB==0){printf("Window function counts completed. No data file provided for power spectrum computation. Exiting now...\n");return 0;}
            if(Ndata>0 && NdataB>0){printf("Window function counts completed.\n\n");}

        }

    }//skycut option





//Periodic
if(strcmp(type_of_survey, "periodic") == 0 || strcmp(type_of_survey, "periodicFKP") == 0 )
{
if(strcmp(type_of_file, "ascii") == 0)
{
printf("Reading %s...",name_data_in);
parameter_value[3]=Ndata;
get_periodic_data(name_data_in, pos_x, pos_y, pos_z, weight, parameter_value,0);
min=parameter_value[10];
max=parameter_value[11];
printf("Ok!\n\n");
if(strcmp(type_of_survey, "periodic") == 0)
{
P_shot_noise=pow(L2-L1,3)/Ndata2*1.*((parameter_value[4]/Ndata2)/((parameter_value[12]/Ndata2)*(parameter_value[12]/Ndata2)));// P_shot_noise = 1/nbar * <w^2> / <w>^2
}
if(strcmp(type_of_survey, "periodicFKP") == 0)
{
alpha_data=parameter_value[12];
Psn_1a=parameter_value[4];
Bsn1a=parameter_value[21];

//randoms

        pos_x_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_y_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_z_rand = (double*) calloc(Nrand2, sizeof(double));
        weight_rand = (double*) calloc(Nrand2, sizeof(double));

printf("Reading %s...",name_randoms_in);
parameter_value[3]=Nrand;
get_periodic_data(name_randoms_in, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, parameter_value,1);
if(parameter_value[10]<min){min=parameter_value[10];}
if(parameter_value[11]>max){max=parameter_value[11];}

alpha_rand=parameter_value[12];
alpha=alpha_data/alpha_rand;

Psn_2a=parameter_value[4];
Bsn2a=parameter_value[21];

I22=alpha_data*alpha_data/pow(L2-L1,3);
I33=alpha_data*alpha_data*alpha_data/pow(L2-L1,6);
IN=I22/I33;
P_shot_noise=(Psn_1a+alpha*alpha*Psn_2a)/I22;
Bsn=(Bsn1a-alpha*alpha*alpha*Bsn2a)/I33;

printf("Ok!\n\n");

}
}
}
    
  if(strcmp(type_of_code, "rusticoX") == 0){
  if(strcmp(type_of_survey, "periodic") == 0 || strcmp(type_of_survey,"periodicFKP") == 0){ 
  if( strcmp(type_of_fileB, "ascii") == 0){
   
      printf("Reading %s...",name_dataB_in);
parameter_valueB[3]=NdataB;
      get_periodic_data(name_dataB_in, pos_xB, pos_yB, pos_zB, weightB, parameter_valueB,0);
      minB=parameter_valueB[10];
      maxB=parameter_valueB[11];
      printf("Ok!\n\n");

if(strcmp(type_of_survey, "periodic") == 0)
{
      P_shot_noiseB=pow(L2-L1,3)/Ndata2B*1.;
}
if(strcmp(type_of_survey, "periodicFKP") == 0)
{
alpha_data1B=parameter_value[12];
Psn_1aB=parameter_value[4];
Bsn1aB=parameter_value[21];

//randoms
        pos_x_randB = (double*) calloc(Nrand2B, sizeof(double));
        pos_y_randB = (double*) calloc(Nrand2B, sizeof(double));
        pos_z_randB = (double*) calloc(Nrand2B, sizeof(double));
        weight_randB = (double*) calloc(Nrand2B, sizeof(double));

printf("Reading %s...",name_randomsB_in);
parameter_value[3]=NrandB;
get_periodic_data(name_randomsB_in, pos_x_randB, pos_y_randB, pos_z_randB, weight_randB, parameter_valueB,1);
if(parameter_valueB[10]<minB){minB=parameter_valueB[10];}
if(parameter_valueB[11]>maxB){maxB=parameter_valueB[11];}

alpha_randB=parameter_valueB[12]; 
alphaB=alpha_dataB/alpha_randB;
    
Psn_2aB=parameter_valueB[4];
Bsn2aB=parameter_valueB[21];

I22B=alpha_dataB*alpha_dataB/pow(L2-L1,3);
I33B=alpha_dataB*alpha_dataB*alpha_dataB/pow(L2-L1,6);
INB=I22B/I33B;

P_shot_noiseB=(Psn_1aB+alphaB*alphaB*Psn_2aB)/I22B;
BsnB=(Bsn1aB-alphaB*alphaB*alphaB*Bsn2aB)/I33B;

printf("Ok!\n\n");



}

    
}
}
}
    



if(max>L2 || min<L1){printf("Warning: Limits of the box are exceeded by the data or random galaxies: data particles found at the limits %lf and %lf. Exiting now...\n",min,max);return 0;}
    
    if(strcmp(type_of_code, "rusticoX") == 0){

        if(maxB>L2 || minB<L1){printf("Warning: Limits of the box-B are exceeded by the data or random galaxies: data particles found at the limits %lf and %lf. Exiting now...\n",minB,maxB);return 0;}
        
    }


/*
if( strcmp(shuffle, "no") == 0 )
{
Nrand2=get_number_used_lines_randoms(name_randoms_in,parameter_value);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2);if(Nrand2==0){exit(0);}}
}
else//shuffle == yes
{
Nrand2=get_number_used_lines_randoms(name_randoms_in,parameter_value);if(Nrand2<1000){printf("Warning, unusual low value for Nrandoms=%ld\n",Nrand2);if(Nrand2==0){exit(0);}}
Nrand2=Nrand;//if shuffle is enabled number of used randoms is the total number of randoms, no matter of the z-cuts (those are applied later if necessary)
}
alpha_rand1=parameter_value[12];
alpha1=alpha_data1/alpha_rand1;


        pos_x_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_y_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_z_rand = (double*) calloc(Nrand2, sizeof(double));
        weight_rand = (double*) calloc(Nrand2, sizeof(double));
*/


}//end of particle if
if( strcmp(type_of_input,"density") == 0){

I22=I22delta;
I33=I33delta;
Ndata2=0;
num_effective=0;
num_effective_rand=0;
num_effective2=0;
num_effective2_rand=0;
num_effective3=0;
num_effective3_rand=0;
z_effective_data=0;
z_effective_rand=0;
Area_survey=0;
alpha=0;
P_shot_noise=Pnoisedelta;
IN=Bnoise1delta;
Bsn=Bnoise2delta;

}
//exit(0);

if(strcmp(write_kvectors,"yes")==0){sprintf(name_ps_kvectors,"%s/Power_Spectrum_kvectors_%s",name_path_out,name_id);}

if( strcmp(do_mu_bins,"yes") == 0  && strcmp(file_for_mu,"yes") == 0 ){
    
    
    if(strcmp(type_of_code, "rustico") == 0){sprintf(name_ps_out,"%s/Power_Spectrum_%s",name_path_out,name_id);}
    if(strcmp(type_of_code, "rusticoX") == 0){
        sprintf(name_ps_out,"%s/Power_SpectrumAA_%s",name_path_out,name_id);
        sprintf(name_psAB_out,"%s/Power_SpectrumAB_%s",name_path_out,name_id);
        sprintf(name_psBB_out,"%s/Power_SpectrumBB_%s",name_path_out,name_id);
      }
         
}
else{
    if(strcmp(type_of_code, "rustico") == 0){sprintf(name_ps_out,"%s/Power_Spectrum_%s.txt",name_path_out,name_id);}
    if(strcmp(type_of_code, "rusticoX") == 0){
        sprintf(name_ps_out,"%s/Power_SpectrumAA_%s.txt",name_path_out,name_id);
        sprintf(name_psAB_out,"%s/Power_SpectrumAB_%s.txt",name_path_out,name_id);
        sprintf(name_psBB_out,"%s/Power_SpectrumBB_%s.txt",name_path_out,name_id);
    }

}

if(strcmp(type_of_code, "rustico") == 0){
sprintf(name_bs_out,"%s/Bispectrum_%s.txt",name_path_out,name_id);

sprintf(name_bs002_out,"%s/Bispectrum_Quadrupole002_%s.txt",name_path_out,name_id);
sprintf(name_bs020_out,"%s/Bispectrum_Quadrupole020_%s.txt",name_path_out,name_id);
sprintf(name_bs200_out,"%s/Bispectrum_Quadrupole200_%s.txt",name_path_out,name_id);
}
    if(strcmp(type_of_code, "rusticoX") == 0){

        sprintf(name_bs_out,"%s/BispectrumAAA_%s.txt",name_path_out,name_id);
        sprintf(name_bsBBB_out,"%s/BispectrumBBB_%s.txt",name_path_out,name_id);
    
        sprintf(name_bsABB_out,"%s/BispectrumABB_%s.txt",name_path_out,name_id);
        sprintf(name_bsBAA_out,"%s/BispectrumBAA_%s.txt",name_path_out,name_id);
        
        sprintf(name_bsBAB_out,"%s/BispectrumBAB_%s.txt",name_path_out,name_id);
        sprintf(name_bsABA_out,"%s/BispectrumABA_%s.txt",name_path_out,name_id);

        sprintf(name_bsAAB_out,"%s/BispectrumAAB_%s.txt",name_path_out,name_id);
        sprintf(name_bsBBA_out,"%s/BispectrumBBA_%s.txt",name_path_out,name_id);

        //bispectrum multipoles not available for ruticoX
        //sprintf(name_bs002_out,"%s/Bispectrum_Quadrupole002_%s.txt",name_path_out,name_id);
        //sprintf(name_bs020_out,"%s/Bispectrum_Quadrupole020_%s.txt",name_path_out,name_id);
        //sprintf(name_bs200_out,"%s/Bispectrum_Quadrupole200_%s.txt",name_path_out,name_id);
    }

sprintf(triangles_id,"%s/Triangles_%s",path_for_triangles,name_id);


for(i=0;i<mubin;i++){

sprintf(name_ps_out2,"%s_%d.txt",name_ps_out,i);
if( strcmp(do_mu_bins,"yes") == 0  && strcmp(file_for_mu,"yes") == 0 ){f=fopen(name_ps_out2,"w");}
else{f=fopen(name_ps_out,"w");}

if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_ps_out);return 0;}

//write header for the power spectrum file
if( strcmp(header, "yes") == 0)
{
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys*wfkp: %lf\n",num_effective3);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys*wfkp: %lf\n",num_effective3_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s\n", type_normalization_mode );}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",I22);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Shot noise factor %lf\n",Shot_noise_factor);}
fprintf(f,"#Shot noise value %lf\n",P_shot_noise);
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0 ){fprintf(f,"#Quadrupole as %s\n",Quadrupole_type);}
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_odd_multipoles,"yes") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Octopole as %s\n",Octopole_type);}
if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Hexadecapole as %s\n",Hexadecapole_type);}

if(strcmp(do_anisotropy,"yes") == 0){
if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"no") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Quadrupole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Dipole\t Quadrupole\t Octopole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
if(strcmp(do_mu_bins, "yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t mu-centerbin\t mu-eff\t P(k,mu)-Pshotnoise\t number of modes\t Pshotnoise\n");}
}else{
fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t number of modes\t Pshotnoise\n");
}

}
fclose(f);

    if(strcmp(type_of_code, "rusticoX") == 0){

        sprintf(name_psBB_out2,"%s_%d.txt",name_psBB_out,i);
        if( strcmp(do_mu_bins,"yes") == 0  && strcmp(file_for_mu,"yes") == 0 ){f=fopen(name_psBB_out2,"w");}
        else{f=fopen(name_psBB_out,"w");}

        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_psBB_out);return 0;}

        //write header for the power spectrum file
        if( strcmp(header, "yes") == 0)
        {
        fprintf(f,"#Data file %s\n",name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s\n",name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld\n",Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol: %lf\n",num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys: %lf\n",num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys*wfkp: %lf\n",num_effective3B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys*wfkp: %lf\n",num_effective3_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf\n",alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s\n", type_normalization_mode );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",I22B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Shot noise factor %lf\n",Shot_noise_factor);}
        fprintf(f,"#Shot noise value %lf\n",P_shot_noiseB);
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSDB); }
        if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0 ){fprintf(f,"#Quadrupole as %s\n",Quadrupole_type);}
        if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_odd_multipoles,"yes") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Octopole as %s\n",Octopole_type);}
        if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Hexadecapole as %s\n",Hexadecapole_type);}

        if(strcmp(do_anisotropy,"yes") == 0){
        if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"no") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Quadrupole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
        if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Dipole\t Quadrupole\t Octopole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
        if(strcmp(do_mu_bins, "yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t mu-centerbin\t mu-eff\t P(k,mu)-Pshotnoise\t number of modes\t Pshotnoise\n");}
        }else{
        fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t number of modes\t Pshotnoise\n");
        }

        }
        fclose(f);
        
        sprintf(name_psAB_out2,"%s_%d.txt",name_psAB_out,i);
        if( strcmp(do_mu_bins,"yes") == 0  && strcmp(file_for_mu,"yes") == 0 ){f=fopen(name_psAB_out2,"w");}
        else{f=fopen(name_psAB_out,"w");}

        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_psAB_out);return 0;}

        //write header for the power spectrum file
        if( strcmp(header, "yes") == 0)
        {
        fprintf(f,"#Data file %s %s\n",name_data_in,name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s %s\n",name_randoms_in,name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld %ld\n",Ndata2,Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld %ld\n",Nrand2,Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf %lf\n",num_effective,num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol: %lf %lf\n",num_effective_rand,num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf %lf\n",num_effective2,num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys: %lf %lf\n",num_effective2_rand,num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys*wfkp: %lf %lf\n",num_effective3,num_effective3B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of random elements weighted by wcol*wsys*wfkp: %lf %lf\n",num_effective3_rand,num_effective3_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf %lf\n",z_effective_data,z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf %lf\n",z_effective_rand,z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf %lf deg^2\n",Area_survey,Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf %lf\n",alpha,alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s\n", type_normalization_mode );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",sqrt(I22*I22B));}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Shot noise factor %lf\n",Shot_noise_factor);}
        fprintf(f,"#Shot noise value 0\n");
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s %s\n",RSD,RSDB); }
        if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0 ){fprintf(f,"#Quadrupole as %s\n",Quadrupole_type);}
        if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(do_odd_multipoles,"yes") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Octopole as %s\n",Octopole_type);}
        if(strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_computation, "DSE") != 0){fprintf(f,"#Hexadecapole as %s\n",Hexadecapole_type);}

        if(strcmp(do_anisotropy,"yes") == 0){
        if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"no") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Quadrupole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
        if(strcmp(do_mu_bins, "no") == 0 && strcmp(do_odd_multipoles,"yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t Dipole\t Quadrupole\t Octopole\t Hexadecapole\t number of modes\t Pshotnoise\n");}
        if(strcmp(do_mu_bins, "yes") == 0){fprintf(f,"#k-centerbin\t k-eff\t mu-centerbin\t mu-eff\t P(k,mu)-Pshotnoise\t number of modes\t Pshotnoise\n");}
        }else{
        fprintf(f,"#k-centerbin\t k-eff\t Monopole-Pshotnoise\t number of modes\t Pshotnoise\n");
        }

        }
        fclose(f);
        
        
    }
        
}

//Write header for the Bispectrum file
if( strcmp(header, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0)
{
f=fopen(name_bs_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");
fclose(f);

    if(strcmp(type_of_code, "rusticoX") == 0){

        ///BBB
        f=fopen(name_bsBBB_out,"w");
        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bsBBB_out);return 0;}
        fprintf(f,"#Data file %s\n",name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s\n",name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld\n",Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf\n",alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",I33B);}
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSDB); }
        fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
        if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
        if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
        if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
        if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
        fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

        fclose(f);
        
        ///AAB
        f=fopen(name_bsAAB_out,"w");
        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bsAAB_out);return 0;}
        fprintf(f,"#Data file %s %s\n",name_data_in,name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s %s\n",name_randoms_in,name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld %ld\n",Ndata2,Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld %ld\n",Nrand2,Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf %lf\n",num_effective,num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf %lf\n",num_effective_rand,num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf %lf\n",num_effective2,num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf %lf\n",num_effective2_rand,num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf %lf\n",z_effective_data,z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf %lf\n",z_effective_rand,z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf %lf deg^2\n",Area_survey,Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf %lf \n",alpha,alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",cbrt(I33*I33*I33B));}
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s %s\n",RSD,RSDB); }
        fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
        if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
        if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
        if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
        if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
        fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

        fclose(f);
        
        ///ABB
        f=fopen(name_bsABB_out,"w");
        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bsABB_out);return 0;}
        fprintf(f,"#Data file %s %s\n",name_data_in,name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s %s\n",name_randoms_in,name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld %ld\n",Ndata2,Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld %ld\n",Nrand2,Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf %lf\n",num_effective,num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf %lf\n",num_effective_rand,num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf %lf\n",num_effective2,num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf %lf\n",num_effective2_rand,num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf %lf\n",z_effective_data,z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf %lf\n",z_effective_rand,z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf %lf deg^2\n",Area_survey,Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf %lf \n",alpha,alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",cbrt(I33*I33B*I33B));}
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s %s\n",RSD,RSDB); }
        fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
        if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
        if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
        if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
        if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
        fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

        fclose(f);
        
        ///BBA
        f=fopen(name_bsBBA_out,"w");
        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bsBBA_out);return 0;}
        fprintf(f,"#Data file %s %s\n",name_data_in,name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s %s\n",name_randoms_in,name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld %ld\n",Ndata2,Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld %ld\n",Nrand2,Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf %lf\n",num_effective,num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf %lf\n",num_effective_rand,num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf %lf\n",num_effective2,num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf %lf\n",num_effective2_rand,num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf %lf\n",z_effective_data,z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf %lf\n",z_effective_rand,z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf %lf deg^2\n",Area_survey,Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf %lf \n",alpha,alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",cbrt(I33B*I33B*I33));}
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s %s\n",RSD,RSDB); }
        fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
        if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
        if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
        if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
        if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
        fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

        fclose(f);
        
        ///ABA
        f=fopen(name_bsABA_out,"w");
        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bsABA_out);return 0;}
        fprintf(f,"#Data file %s %s\n",name_data_in,name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s %s\n",name_randoms_in,name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld %ld\n",Ndata2,Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld %ld\n",Nrand2,Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf %lf\n",num_effective,num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf %lf\n",num_effective_rand,num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf %lf\n",num_effective2,num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf %lf\n",num_effective2_rand,num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf %lf\n",z_effective_data,z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf %lf\n",z_effective_rand,z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf %lf deg^2\n",Area_survey,Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf %lf \n",alpha,alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",cbrt(I33*I33*I33B));}
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s %s\n",RSD,RSDB); }
        fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
        if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
        if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
        if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
        if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
        fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

        fclose(f);
        
        ///BAB
        f=fopen(name_bsBAB_out,"w");
        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bsBAB_out);return 0;}
        fprintf(f,"#Data file %s %s\n",name_data_in,name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s %s\n",name_randoms_in,name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld %ld\n",Ndata2,Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld %ld\n",Nrand2,Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf %lf\n",num_effective,num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf %lf\n",num_effective_rand,num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf %lf\n",num_effective2,num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf %lf\n",num_effective2_rand,num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf %lf\n",z_effective_data,z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf %lf\n",z_effective_rand,z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf %lf deg^2\n",Area_survey,Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf %lf \n",alpha,alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",cbrt(I33B*I33*I33B));}
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s %s\n",RSD,RSDB); }
        fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
        if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
        if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
        if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
        if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
        fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

        fclose(f);
        
        ///BAA
        f=fopen(name_bsBAA_out,"w");
        if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bsBAA_out);return 0;}
        fprintf(f,"#Data file %s %s\n",name_data_in,name_dataB_in);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s %s\n",name_randoms_in,name_randomsB_in);}
        fprintf(f,"#Number of data elements used: %ld %ld\n",Ndata2,Ndata2B);
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld %ld\n",Nrand2,Nrand2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf %lf\n",num_effective,num_effectiveB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf %lf\n",num_effective_rand,num_effective_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf %lf\n",num_effective2,num_effective2B);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf %lf\n",num_effective2_rand,num_effective2_randB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf %lf\n",z_effective_data,z_effective_dataB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf %lf\n",z_effective_rand,z_effective_randB);}
        fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf %lf deg^2\n",Area_survey,Area_surveyB);}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf %lf \n",alpha,alphaB);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
        if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",cbrt(I33*I33*I33B));}
        fprintf(f,"#Type of Computation: %s\n",type_of_computation);
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
        if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
        if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
        if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
        if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s %s\n",RSD,RSDB); }
        fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
        if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
        if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
        if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
        if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
        fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");

        fclose(f);
        
        

    }
}
    
//Bispectrum quadrupole
if( strcmp(header, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0 &&  strcmp(do_bispectrum2, "yes") == 0)
{
f=fopen(name_bs002_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs002_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum002-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");
fclose(f);
///
f=fopen(name_bs020_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs020_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum020-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");
fclose(f);
///
f=fopen(name_bs200_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs200_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %ld\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Number of random elements used: %ld\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %lf\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol: %lf\n",num_effective_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of randoms elements weighted by wcol*wsys: %lf\n",num_effective2_rand);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){fprintf(f,"#Normalization %e\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %d\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf OL=%lf\n",Omega_m,Omega_L);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fprintf(f,"#k1-centerbin\t k1-eff\t k2-centerbin\t k2-eff\t k3-centerbin\t k3-eff\t Bispectrum200-Bshotnoise\t Bsn\t Q-Qsn\t Qsn\t Number of Triangles\n");
fclose(f);

}
//exit(0);
//Start Power Spectrum Computation for Cutsky
printf("== Computing the Power Spectrum ==\n");

                                                          
//Compute and write the Power Spectrum for FFT+Skycut type of survey.
if(strcmp(type_of_computation, "FFT") == 0){

if(strcmp(type_of_survey, "cutsky") == 0 || strcmp(type_of_survey, "periodicFKP") == 0){
parameter_value[0]=L1;
parameter_value[1]=L2;
if(strcmp(type_of_survey, "cutsky") == 0){check_box_for_yamamoto(parameter_value,ngrid);}

L1=parameter_value[0];
L2=parameter_value[1];

if(strcmp(type_of_code, "rusticoX") == 0){
parameter_valueB[0]=L1;
parameter_valueB[1]=L2;
}

    
//AA AB BB
    kf=2.*(4.*atan(1.))/(L2-L1);
    kny=2.*(4.*atan(1.))/(L2-L1)*ngrid/2.;
    
if(strcmp(type_of_yamamoto, "GridCenter") == 0){

    //if(kf>kmin){printf("Warning, computing power spectrum for k>%lf, as kf>kmin.",kf);}
    //else{kf=kmin;}
    //if(kmax>kny){printf("Warning, computing power spectrum for k<%lf, as kny<kmax.",kny);}
    //else{kny=kmax;}
    kf=kmin;kny=kmax;//overwrites previous cuts

if(strcmp(type_of_input,"density") == 0){strcpy (name_out_density,name_data_in);}

loop_interlacing_skycut(kf,kny,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2,pos_xB, pos_yB, pos_zB, weightB,Ndata2B, pos_x_randB, pos_y_randB, pos_z_randB, weight_randB,Nrand2B, L1, L2, ngrid, P_shot_noise,P_shot_noiseB, bin_ps, I22,I22B, alpha,alphaB, mode_correction, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_ps_out,name_psAB_out,name_psBB_out, type_of_mass_assigment,do_bispectrum,type_of_code,output_density,name_out_density,name_out_densityB,type_of_input,type_of_survey,file_for_mu,mubin,write_kvectors, name_ps_kvectors);
        

}

if(strcmp(type_of_yamamoto, "GridAverage") == 0){

//if(strcmp(window_function, "yes") == 0){

    //if(kf>kmin){printf("Warning, computing power spectrum for k>%lf, as kf>kmin.",kf);}
    //else{kf=kmin;}
    //if(kmax>kny){printf("Warning, computing power spectrum for k<%lf, as kny<kmax.",kny);}
    //else{kny=kmax;}
    kf=kmin;kny=kmax;//overwrites previous cuts

loop_interlacing_skycut2(kf,kny,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2, pos_xB, pos_yB, pos_zB, weightB,Ndata2B, pos_x_randB, pos_y_randB, pos_z_randB, weight_randB,Nrand2B, L1, L2, ngrid, P_shot_noise,P_shot_noiseB, bin_ps, I22,I22B, alpha,alphaB, mode_correction, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_ps_out,name_psAB_out,name_psBB_out, type_of_mass_assigment,do_bispectrum,type_of_code,output_density,name_out_density,name_out_densityB,type_of_survey,file_for_mu,mubin,write_kvectors, name_ps_kvectors);


}

}
}
                                                         
                                                          
if(strcmp(type_of_code, "rustico") == 0){

                                                          
if(strcmp(type_of_computation, "DSY") == 0 || strcmp(type_of_computation, "DSE") == 0)
{

     if( strcmp(type_of_survey, "cutsky") == 0)
     {
         loop_directsum_yamamoto_skycut_caller(kmin,kmax,pos_x, pos_y, pos_z, weight, Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2, L1, L2, P_shot_noise, bin_ps, I22, alpha, n_lines_parallel, binning_type, Quadrupole_type,Octopole_type,Hexadecapole_type,do_odd_multipoles,do_anisotropy, name_ps_out,type_of_computation,write_kvectors, name_ps_kvectors);         
     }
     if(strcmp(type_of_survey, "periodic") == 0 || strcmp(type_of_survey, "periodicFKP") == 0)
     {
          printf("No Direct Sum for periodic box at the moment. Exiting now...\n");return 0;
     }

}
            
}
if(strcmp(type_of_code, "rusticoX") == 0){
            
            if(strcmp(type_of_computation, "DSY") == 0 || strcmp(type_of_computation, "DSE") == 0)
            {

                 if( strcmp(type_of_survey, "cutsky") == 0)
                 {
                      printf("No Direct Sum for X-analysis at the moment. Exiting now...\n");return 0;
                               
                 }
                 if(strcmp(type_of_survey, "periodic") == 0 || strcmp(type_of_survey, "periodicFKP") == 0)
                 {
                      printf("No Direct Sum for periodic box at the moment. Exiting now...\n");return 0;
                 }

            }

}
                                                          

//Compute and write the Power Spectrum for FFT+Box with constant line of sight along z.
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "periodic") == 0)
{
    //if(kf>kmin){printf("Warning, computing power spectrum for k>%lf, as kf>kmin.",kf);}
    //else{kf=kmin;}
    //if(kmax>kny){printf("Warning, computing power spectrum for k<%lf, as kny<kmax.",kny);}
    //else{kny=kmax;}
    kf=kmin;kny=kmax;//overwrites previous cuts

if(strcmp(type_of_code, "rustico") == 0){

if(strcmp(type_of_file, "ascii") == 0){

if(strcmp(type_of_input,"density") == 0){strcpy (name_out_density,name_data_in);}

        loop_interlacing_periodic(kf,kny,Ninterlacing, pos_x, pos_y, pos_z, weight, Ndata2,Ndata2w,NULL, NULL, NULL, NULL, 0,0, L1, L2, ngrid, P_shot_noise,0, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out,NULL,NULL, type_of_mass_assigment,do_odd_multipoles,do_anisotropy,do_bispectrum,mubin,file_for_mu,type_of_code,output_density,name_out_density,name_out_densityB,type_of_input,write_kvectors, name_ps_kvectors);
}
if(strcmp(type_of_file, "gadget") == 0){
        loop_interlacing_periodic_gadget(kf,kny,Ninterlacing, name_data_in,NULL ,gadget_files,0, L1, L2, ngrid, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out,NULL,NULL,  type_of_mass_assigment,Shot_noise_factor,grid_correction_string,RSD,NULL,do_odd_multipoles,do_anisotropy,mubin,file_for_mu,type_of_code,output_density,name_out_density,name_out_densityB,type_of_input,write_kvectors, name_ps_kvectors);
}

}


if(strcmp(type_of_code, "rusticoX") == 0){

    if(strcmp(type_of_file, "ascii") == 0 && strcmp(type_of_fileB, "ascii") == 0)
    {
        loop_interlacing_periodic(kf,kny,Ninterlacing, pos_x, pos_y, pos_z, weight, Ndata2,Ndata2w,pos_xB, pos_yB, pos_zB, weightB, Ndata2B,Ndata2Bw, L1, L2, ngrid, P_shot_noise,P_shot_noiseB, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out,name_psAB_out,name_psBB_out, type_of_mass_assigment,do_odd_multipoles,do_anisotropy,do_bispectrum,mubin,file_for_mu,type_of_code,output_density,name_out_density,name_out_densityB,type_of_input,write_kvectors, name_ps_kvectors);


    }
    if(strcmp(type_of_file, "gadget") == 0 && strcmp(type_of_fileB, "gadget") == 0)
    {
        loop_interlacing_periodic_gadget(kf,kny,Ninterlacing, name_data_in,name_dataB_in ,gadget_files,gadget_filesB, L1, L2, ngrid, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out,name_psAB_out,name_psBB_out,  type_of_mass_assigment,Shot_noise_factor,grid_correction_string,RSD,RSDB,do_odd_multipoles,do_anisotropy,mubin,file_for_mu,type_of_code,output_density,name_out_densityB,name_out_density,type_of_input,write_kvectors, name_ps_kvectors);

    }
    if(strcmp(type_of_file, "gadget") == 0 && strcmp(type_of_fileB, "ascii") == 0)
    {
        loop_interlacing_periodic_gadget_x_ascii(kf,kny,Ninterlacing, name_data_in, gadget_files, pos_xB, pos_yB, pos_zB, weightB, Ndata2B,Ndata2Bw, L1, L2, ngrid,P_shot_noiseB, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out,name_psAB_out,name_psBB_out, type_of_mass_assigment,Shot_noise_factor,do_bispectrum,RSD,do_odd_multipoles,do_anisotropy,mubin,file_for_mu,type_of_code,output_density,name_out_densityB,name_out_density,write_kvectors, name_ps_kvectors);

    }
    if(strcmp(type_of_file, "ascii") == 0 && strcmp(type_of_fileB, "gadget") == 0)
    {
        loop_interlacing_periodic_gadget_x_ascii(kf,kny,Ninterlacing, name_dataB_in, gadget_filesB, pos_x, pos_y, pos_z, weight, Ndata2,Ndata2w, L1, L2, ngrid,P_shot_noise, bin_ps, mode_correction, n_lines_parallel, binning_type, name_psBB_out,name_psAB_out,name_ps_out, type_of_mass_assigment,Shot_noise_factor,do_bispectrum,RSDB,do_odd_multipoles,do_anisotropy,mubin,file_for_mu,type_of_code,output_density,name_out_densityB,name_out_density,write_kvectors, name_ps_kvectors);

    }

}
    
}

if(strcmp(do_bispectrum, "no") == 0){
printf("Computation of Power Spectrum finished sucessfully!\n\n");
return 0;
}

printf("== Computing the Bispectrum ==\n");

if(strcmp(type_of_computation, "FFT") == 0){
if(strcmp(type_of_survey,"periodicFKP") == 0 || strcmp(type_of_survey, "cutsky") == 0){

//printf("%e %e %e\n",I33,IN,Bsn);
    
  loop_bispectrum_skycut_caller(kf, kny, Ninterlacing,  pos_x, pos_y, pos_z, weight, Ndata2, 0, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2,pos_xB, pos_yB, pos_zB, weightB, Ndata2B,0, pos_x_randB, pos_y_randB, pos_z_randB, weight_randB, Nrand2B , L1, L2, ngrid, P_shot_noise, P_shot_noiseB, Deltakbis, I33,I33B,I22,I22B, IN, INB, Bsn,BsnB, alpha,alphaB, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bsAAB_out,name_bsABA_out,name_bsBAA_out,name_bsABB_out,name_bsBAB_out,name_bsBBA_out,name_bsBBB_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment,triangles_num,write_triangles,triangles_id, do_multigrid,bispectrum_optimization, triangle_shapes,do_bispectrum2,type_of_code,name_data_in,type_of_input,type_of_survey);


}   
}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "periodic") == 0)
{
if(strcmp(type_of_file, "ascii") == 0)
{
    
    if(strcmp(type_of_code, "rustico") == 0){

   
    loop_bispectrum_skycut_caller(kf, kny, Ninterlacing,  pos_x, pos_y, pos_z, weight, Ndata2,Ndata2w, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2,pos_xB, pos_yB, pos_zB, weightB, Ndata2B,Ndata2Bw, pos_x_randB, pos_y_randB, pos_z_randB, weight_randB, Nrand2B , L1, L2, ngrid, P_shot_noise, P_shot_noiseB, Deltakbis, 0,0,0,0, 0, 0, 0,0, 0,0, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bsAAB_out,name_bsABA_out,name_bsBAA_out,name_bsABB_out,name_bsBAB_out,name_bsBBA_out,name_bsBBB_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment,triangles_num,write_triangles,triangles_id, do_multigrid, bispectrum_optimization, triangle_shapes,do_bispectrum2,type_of_code,name_data_in,type_of_input,type_of_survey);
    }
    if(strcmp(type_of_code, "rusticoX") == 0){
    
        if(strcmp(type_of_fileB, "ascii") == 0){
            
            loop_bispectrum_skycut_caller(kf, kny, Ninterlacing,  pos_x, pos_y, pos_z, weight, Ndata2,Ndata2w, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2,pos_xB, pos_yB, pos_zB, weightB, Ndata2B,Ndata2Bw, pos_x_randB, pos_y_randB, pos_z_randB, weight_randB, Nrand2B , L1, L2, ngrid, P_shot_noise, P_shot_noiseB, Deltakbis, 0,0,0,0, 0, 0, 0,0, 0,0, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bsAAB_out,name_bsABA_out,name_bsBAA_out,name_bsABB_out,name_bsBAB_out,name_bsBBA_out,name_bsBBB_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment,triangles_num,write_triangles,triangles_id, do_multigrid, bispectrum_optimization, triangle_shapes,do_bispectrum2,type_of_code,NULL,type_of_input,type_of_survey);

            
        }
        if(strcmp(type_of_fileB, "gadget") == 0){

            //cross ascii x gadget
            reverse=1;
            loop_bispectrum_periodic_for_gadget_x_ascii_caller(kf,kny,Ninterlacing, L1, L2, ngrid, Deltakbis, mode_correction, n_lines_parallel, binning_type, name_bs_out, name_bsAAB_out, name_bsABA_out, name_bsBAA_out,name_bsABB_out, name_bsBAB_out,name_bsBBA_out, name_bsBBB_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment, triangles_num, write_triangles, triangles_id, name_dataB_in, gadget_filesB,  pos_x, pos_y, pos_z, weight, Ndata2,Ndata2w, P_shot_noise,do_multigrid, bispectrum_optimization, triangle_shapes,RSDB,do_bispectrum2,type_of_code,reverse,type_of_survey);
                       
        }

    
    
    }
        

    

}
if(strcmp(type_of_file, "gadget") == 0)
{
    
    if(strcmp(type_of_code, "rustico") == 0){

loop_bispectrum_periodic_for_gadget_caller(kf,kny,Ninterlacing, L1, L2, ngrid, Deltakbis, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bsAAB_out,name_bsABA_out,name_bsBAA_out,name_bsABB_out,name_bsBAB_out,name_bsBBA_out,name_bsBBB_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment, triangles_num, write_triangles, triangles_id, name_data_in,gadget_files,name_dataB_in,gadget_filesB, do_multigrid, bispectrum_optimization, triangle_shapes,RSD,RSDB,do_bispectrum2,type_of_code,type_of_survey);


    }
    if(strcmp(type_of_code, "rusticoX") == 0){
        
        if(strcmp(type_of_fileB, "gadget") == 0)
        {
            //gadget x gadget
            loop_bispectrum_periodic_for_gadget_caller(kf,kny,Ninterlacing, L1, L2, ngrid, Deltakbis, mode_correction, n_lines_parallel, binning_type, name_bs_out,name_bsAAB_out,name_bsABA_out,name_bsBAA_out,name_bsABB_out,name_bsBAB_out,name_bsBBA_out,name_bsBBB_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment, triangles_num, write_triangles, triangles_id, name_data_in,gadget_files,name_dataB_in,gadget_filesB, do_multigrid, bispectrum_optimization, triangle_shapes,RSD,RSDB,do_bispectrum2,type_of_code,type_of_survey);

        }
        if(strcmp(type_of_fileB, "ascii") == 0){
         
            //gadget x ascii
            reverse=0;
            loop_bispectrum_periodic_for_gadget_x_ascii_caller(kf,kny,Ninterlacing, L1, L2, ngrid, Deltakbis, mode_correction, n_lines_parallel, binning_type, name_bs_out, name_bsAAB_out, name_bsABA_out, name_bsBAA_out,name_bsABB_out, name_bsBAB_out,name_bsBBA_out, name_bsBBB_out,name_bs002_out,name_bs020_out,name_bs200_out, type_of_mass_assigment, triangles_num, write_triangles, triangles_id, name_data_in, gadget_files,  pos_xB, pos_yB, pos_zB, weightB, Ndata2B,Ndata2Bw, P_shot_noiseB,do_multigrid, bispectrum_optimization, triangle_shapes,RSD,do_bispectrum2,type_of_code,reverse,type_of_survey);

            
        }


    }
    
}

}
if(strcmp(type_of_computation, "FFT") != 0)
{
printf("No direct Sum for the Bispectrum yet!\n");
return 0;
}


//Improve: sorting algorithm. Avoid system calls

printf("Computation of the Bispectrum finished sucessfully!\n\n");
return 0;
}
