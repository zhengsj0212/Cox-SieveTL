#!/bin/bash
#$ -N sievetl_cov1
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_rt=08:00:00,h_data=6G
#$ -pe shared 1
#$ -t 1-500

mkdir -p logs

source /u/local/Modules/default/init/bash

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export NSIM=${NSIM:-500}
export REPS_PER_TASK=${REPS_PER_TASK:-1}
export OUT_BASE=${OUT_BASE:-$PWD/coverage_plugin_wald_results}
export BETA2=${BETA2:-0.5}
export KAPPA=${KAPPA:-2}
export C1=${C1:-1}
export SHIFT=${SHIFT:-0}
export LVAL=${LVAL:-5}
export N_TARGET=${N_TARGET:-100}
export N_SOURCE=${N_SOURCE:-1000}
export MULTIPLIER_BOOT=${MULTIPLIER_BOOT:-1000}
export CONF_LEVEL=${CONF_LEVEL:-0.95}
export SEED_BASE=${SEED_BASE:-100000}

module load R/4.2.2

Rscript run_coverage_plugin_wald_hoffman2.R
