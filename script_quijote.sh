#!/bin/bash
#SBATCH -p normal #partition (queue) ##Don't change this
#SBATCH -N 1 #number of nodes
#SBATCH -n 20 #number of cores
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu 1000 #memory pool ##Mb of RAM
#SBATCH -t 20-20:00 # time limit (D-HH:MM)
#SBATCH -o /home/DATA/hector/codes/code_brassv15_example/log.o #path file to store output
#SBATCH -e /home/DATA/hector/codes/code_brassv15_example/log.e #path file to store error messages
#SBATCH --job-name shapefit #job name

export OMP_NUM_THREADS=20

time ./file_gcc_quijote.out params_quijote.c
