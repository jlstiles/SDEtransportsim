#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M jlstiles@berkeley.edu 
#$ -m beas
export OMP_NUM_THREADS=1

module load gcc-4.9.3

R --vanilla < sim_new_YZ.R > sim_new_YZ.Rout


