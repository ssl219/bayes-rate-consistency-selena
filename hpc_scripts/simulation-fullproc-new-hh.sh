#!/bin/sh

#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:ompthreads=1:mem=50gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
WAVE=1
MODEL="hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson"
HSGP_C=1.5
HSGP_M=20
SIZE=55
HHSIZE=4
SCENARIO="flat"

# HMC Sampler params
CHAINS=4
WARMUP=500
SAMPLING=1000

# Post-processing
MIXING=TRUE

module load anaconda3/personal
source activate Renv

Rscript $REPO_PATH/scripts/sim-run-stan-single-new-hh.R --wave $WAVE --model $MODEL --hsgp_c $HSGP_C --hsgp_m $HSGP_M --chains $CHAINS --size $SIZE --hhsize $HHSIZE --scenario $SCENARIO --iter_warmup $WARMUP --iter_sampling $SAMPLING

MODEL=${MODEL}-${WAVE}-"sim-hh"${HHSIZE}-${SCENARIO}-${SIZE}
Rscript $REPO_PATH/scripts/sim-postprocess-single-new-hh.R --model $MODEL --size $SIZE --hhsize $HHSIZE --scenario $SCENARIO --mixing $MIXING
