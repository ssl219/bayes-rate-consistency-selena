#!/bin/sh

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=100gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
WAVE=1
MODEL="hsgp-eq-cd-new-hh-dropping-all-zeros-symmetric-poisson"
HSGP_C=1.5
HSGP_M=20

# HMC Sampler params
CHAINS=1
WARMUP=50
SAMPLING=100

# Post-processing
MIXING=TRUE

module load anaconda3/personal
source activate Renv

Rscript $REPO_PATH/scripts/sim-run-stan-single-new-hh-boarding.R --wave $WAVE --model $MODEL --hsgp_c $HSGP_C --hsgp_m $HSGP_M --chains $CHAINS --iter_warmup $WARMUP --iter_sampling $SAMPLING

MODEL=${MODEL}-${WAVE}
Rscript $REPO_PATH/scripts/sim-postprocess-single-new-hh-boarding.R --model $MODEL --mixing $MIXING
