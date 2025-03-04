#!/bin/sh

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=500gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
WAVE=1
MODEL="hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson-multiple-betas-1"

# Post-processing
MIXING=TRUE
PPC=TRUE
PLOT=TRUE

module load anaconda3/personal
source activate Renv

Rscript $REPO_PATH/scripts/postprocess-single-new-hh-poisson-multiple-betas.R --model $MODEL --mixing $MIXING --ppc $PPC --plot $PLOT
