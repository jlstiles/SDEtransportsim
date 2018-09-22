#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
#$ -M jlstiles@berkeley.edu 
#$ -m beas

R --vanilla < sim_kara_crazy.R > sim_kara_crazy.Rout


