#!/bin/bash
#
#SBATCH -J sim
#SBATCH -p long 
#SBATCH -t 70-10:00:00
#SBATCH -o sim.%j.out
#SBATCH -e sim.%j.err
#SBATCH -n 16
#SBATCH --mem=12000
module purge 
module load gnu8 R
Rscript /home/gdeliorm/high_ica_1000_newsim.R