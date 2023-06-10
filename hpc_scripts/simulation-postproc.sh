#!/bin/sh

#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:ompthreads=1:mem=200gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
DATA_PATH=/rds/general/user/ssl219/home
WAVE=1
MODEL="hsgp-eq-rd-new-hh-dropping-all-zeros-symmetric-poisson"
HSGP_C=1.5
HSGP_M=20
N=500
HHSIZE=2
SCENARIO="flat"
SIM_NO=1

# HMC Sampler params
CHAINS=4
WARMUP=500
SAMPLING=1000

# Post-processing
MIXING=TRUE

module load anaconda3/personal
source activate Renv

MODEL=${MODEL}-${WAVE}-"sim-hh"${HHSIZE}-${SCENARIO}-${N}-${SIM_NO}
Rscript $REPO_PATH/scripts/sim-postprocess-single-new-hh.R --model $MODEL --sample_size $N --hhsize $HHSIZE --scenario $SCENARIO --mixing $MIXING --sim.no $SIM_NO
