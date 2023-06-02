#!/bin/sh

#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:ompthreads=1:mem=50gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
DATA_PATH=/rds/general/user/ssl219/home
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

Rscript $REPO_PATH/scripts/sim-run-stan-single-new-hh.R --iter_warmup $WARMUP --iter_sampling $SAMPLING --chains $CHAINS --hhsize $HHSIZE --size $SIZE--scenario $SCENARIO --model $MODEL --repo_path $REPO_PATH --data_path $DATA_PATH --hsgp_c $HSGP_C --hsgp_m $HSGP_M --wave $WAVE

MODEL=${MODEL}-${WAVE}-"sim-hh"${HHSIZE}-${SCENARIO}-${SIZE}
Rscript $REPO_PATH/scripts/sim-postprocess-single-new-hh.R --model $MODEL --size $SIZE --hhsize $HHSIZE --scenario $SCENARIO --mixing $MIXING
