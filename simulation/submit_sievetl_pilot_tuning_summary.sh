#!/bin/bash
#$ -N sievetl_pilot_sum
#$ -cwd
#$ -j y
#$ -o logs/$JOB_NAME.$JOB_ID.log
#$ -l h_rt=01:00:00,h_data=4G
#$ -pe shared 1

mkdir -p logs

source /u/local/Modules/default/init/bash

export PILOT_OUT_BASE=${PILOT_OUT_BASE:-$PWD}

module load R/4.2.2

Rscript summarize_spline_sievetl_pilot_tuning.R
