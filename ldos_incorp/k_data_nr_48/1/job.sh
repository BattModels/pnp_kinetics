#!/bin/bash
#SBATCH -J etarun_1_nr48
#SBATCH -n 64 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH -p RM-shared # Partition to submit to
#SBATCH --mem=126000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o out.out # File to which STDOUT will be written
#SBATCH -e out.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=mbabar@andrew.cmu.edu # Email to which notifications will be sent

echo "Job started on `hostname` at `date`"

julia eta_run_script.jl > out.out

echo " "
echo "Job Ended at `date`"
