#!/bin/bash
#$ -N sievetl_invh_col
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_rt=08:00:00,h_data=4G
#$ -pe shared 1
#$ -t 1-20

mkdir -p logs

source /u/local/Modules/default/init/bash

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export NPILOT=${NPILOT:-20}
export PILOT_SEED_BASE=${PILOT_SEED_BASE:-1}
export PILOT_OUT_BASE=${PILOT_OUT_BASE:-$PWD}

export BETA2=${BETA2:-0.5}
export KAPPA=${KAPPA:-2}
export C1=${C1:-1}
export SHIFT=${SHIFT:-0}
export LVAL=${LVAL:-5}
export N_TARGET=${N_TARGET:-100}
export N_SOURCE=${N_SOURCE:-1000}

module load R/4.2.2

Rscript run_PLC_CR_spline_hoffman2_pilot_invh_tuning_by_column.R
