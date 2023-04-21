#Examples of compilation. You will need to edit get_line.c to read the data format. You might have to edit the header of read.c if you want to re-scale the covariance from the mocks (for eg., when fitting the mean of X mocks). The file priors.c contains information on the priors for each parameters in the function set_mcmc_priors(). If you do not use a proposal covariance, you may also want to edit the functions set_proposal_mean and set_proposal_error for a first point and size of the step for the first mcmc trial run


icc main.c bao.c bispectrum.c cubature.c fftlog.c functions.c get_line.c integrals_rsd.c mcmc_bao.c priors.c read.c rsd.c structures.c -O3 -lgsl -lm -lgslcblas   -openmp  -lpthread  -lfftw3   -I/users/hectorgm/fftw3_threads/include/ -L/users/hectorgm/fftw3_threads/lib/  -I/OPT/gsl-1.16/include/ -L//OPT/gsl-1.16/lib/ -o file_icc.out

gcc main.c bao.c bispectrum.c cubature.c fftlog.c functions.c get_line.c integrals_rsd.c mcmc_bao.c priors.c read.c rsd.c structures.c -O3 -lgsl -lm -lgslcblas   -fopenmp  -lpthread  -lfftw3 -o file_gcc.out
