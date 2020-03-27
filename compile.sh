export OMP_NUM_THREADS=10

gcc main.c functions.c read.c bao.c mcmc_bao.c fftlog.c rsd.c integrals_rsd.c structures.c priors.c cubature.c -O3 -lgsl -lm -lgslcblas   -fopenmp  -lpthread  -lfftw3 -o file.out

