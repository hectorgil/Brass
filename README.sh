
===== Quick Introduction on how to use BRASS =====

DATE: 7th Feb. 2023
Contact: Hector Gil Marin hectorgil@icc.ub.edu
Please cite 2007.08994 if you are employing this code (even partially as a validation pipeline) for a publication. 

#Be aware you will need to edit get_line.c to read the data format BEFORE you compile your executable file. You might also have to edit (BEFORE compiling) the header of read.c if you want to re-scale the covariance from the mocks (for eg., when fitting the mean of X mocks).

#define SCALING_COVARIANCE (X)


The file priors.c contains information on the priors for each parameters in the function set_mcmc_priors(). If you do not use a proposal covariance, you may also want to edit the functions set_proposal_mean and set_proposal_error for a first point and size of the step for the first mcmc trial run.

#Please read the readme_params for a quick reference on the inputs of the params.ini file

#Examples of compilation with icc and gcc compiler linking to the GSL libraries

icc main.c bao.c bispectrum.c cubature.c fftlog.c functions.c get_line.c integrals_rsd.c mcmc_bao.c priors.c read.c rsd.c structures.c -O3 -lgsl -lm -lgslcblas   -openmp  -lpthread  -lfftw3   -I/users/hectorgm/fftw3_threads/include/ -L/users/hectorgm/fftw3_threads/lib/  -I/OPT/gsl-1.16/include/ -L//OPT/gsl-1.16/lib/ -o file_icc.out

gcc main.c bao.c bispectrum.c cubature.c fftlog.c functions.c get_line.c integrals_rsd.c mcmc_bao.c priors.c read.c rsd.c structures.c -O3 -lgsl -lm -lgslcblas   -fopenmp  -lpthread  -lfftw3 -o file_gcc.out


For running the file you need to define the number of threads you are going to use through,

$export OMP_NUM_THREADS=20
and then run the code as,

$time ./file_gcc.out params.c


You can find the necessary files for a test run in,

https://drive.google.com/drive/folders/1nS1f9Lx8zIdsEXGnKuVW-pOglvAHlkuz?usp=share_link

