#!/bin/sh

#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:ompthreads=1:mem=50gb

REPO_PATH=/rds/general/user/ssl219/home/bayes-rate-consistency-selena
DATA_PATH=/rds/general/user/ssl219/home
WAVE=1
MODEL="hsgp-eq-rd"
HSGP_C=1.5
HSGP_M=20
N=500
HHSIZE=2
SCENARIO="boarding_school"
SIM_NO=1

# HMC Sampler params
CHAINS=4
WARMUP=500
SAMPLING=1000

# Post-processing
MIXING=TRUE

module load anaconda3/personal
source activate Renv

Rscript $REPO_PATH/scripts/sim-run-stan-single-baseline.R --sim.no $SIM_NO --iter_warmup $WARMUP --iter_sampling $SAMPLING --chains $CHAINS --hhsize $HHSIZE --sample_size $N --scenario $SCENARIO --model $MODEL --repo_path $REPO_PATH --data_path $DATA_PATH --hsgp_c $HSGP_C --hsgp_m $HSGP_M --wave $WAVE

MODEL=${MODEL}-${WAVE}-"sim-hh"${HHSIZE}-${SCENARIO}-${N}-${SIM_NO}
Rscript $REPO_PATH/scripts/sim-postprocess-single-baseline.R --model $MODEL --sample_size $N --hhsize $HHSIZE --scenario $SCENARIO --mixing $MIXING --sim.no $SIM_NO
