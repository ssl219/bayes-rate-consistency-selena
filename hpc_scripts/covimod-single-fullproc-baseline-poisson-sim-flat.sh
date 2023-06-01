#!/bin/sh

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=200gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
WAVE=1
MODEL="hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson"
HSGP_C=1.5
HSGP_M=20
SAMPLE_SIZE=450

# HMC Sampler params
CHAINS=4
WARMUP=500
SAMPLING=1000

# Post-processing
MIXING=TRUE

module load anaconda3/personal
source activate Renv

Rscript $REPO_PATH/scripts/sim-run-stan-single-baseline-flat.R --wave $WAVE --model $MODEL --hsgp_c $HSGP_C --hsgp_m $HSGP_M --chains $CHAINS --iter_warmup $WARMUP --iter_sampling $SAMPLING

MODEL=${MODEL}-${WAVE}-"sim-flat-baseline"-${SAMPLE_SIZE}
Rscript $REPO_PATH/scripts/sim-postprocess-single-baseline-flat.R --model $MODEL --mixing $MIXING
