#!/bin/sh

#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:ompthreads=1:mem=50gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
WAVE=1
MODEL="hsgp-eq-cd-new-hh-2-symmetric-poisson"

# Post-processing
MIXING=TRUE
PPC=TRUE
PLOT=TRUE

module load anaconda3/personal
source activate Renv

MODEL=${MODEL}-${WAVE}
Rscript $REPO_PATH/scripts/postprocess-single-new-hh-2-poisson.R --model $MODEL --mixing $MIXING --ppc $PPC --plot $PLOT
