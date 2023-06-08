#!/bin/sh

#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:ompthreads=1:mem=500gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
WAVE=1
MODEL="hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson-1-hhgreater5"

# Post-processing
MIXING=FALSE
PPC=TRUE
PLOT=TRUE

module load anaconda3/personal
source activate Renv

Rscript $REPO_PATH/scripts/postprocess-single-new-hh-poisson-hh5.R --model $MODEL --mixing $MIXING --ppc $PPC --plot $PLOT
