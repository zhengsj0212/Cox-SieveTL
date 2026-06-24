# Simulation Code for Cox-SieveTL

This folder contains the cleaned simulation workflow for the paper:

**Cox-SieveTL: Semiparametric Transfer Learning for Cox Models via Sieve Maximum Likelihood**

The code here keeps only the main simulation pipeline for the SieveTL method. It is designed to:
- tune the main SieveTL penalties from pilot replications
- tune the by-column inverse-Hessian penalties from pilot replications
- run the final coverage / inference jobs
- save raw `.rds` outputs

This folder does **not** contain:
- comparison methods
- plotting code
- summary CSV generation
- LaTeX table generation
- filtered postprocessing

## Workflow Overview

Run the code in this order:

1. stage-1 pilot tuning for the main SieveTL penalties
2. stage-1 pilot summary
3. stage-2 pilot tuning for the by-column inverse-Hessian penalties
4. stage-2 pilot summary
5. final coverage / inference jobs

The final job reads the tuning files produced in steps 2 and 4 and writes raw `.rds` outputs.

## Step-by-Step Procedure

### Step 1. Stage-1 pilot tuning for the main SieveTL penalties

Run:
- [submit_sievetl_pilot_tuning_array.sh](submit_sievetl_pilot_tuning_array.sh)

This script launches Hoffman2 array jobs for the stage-1 pilot replications.

By default it uses:
- `NPILOT = 20`

The underlying R code is:
- [run_PLC_CR_spline_hoffman2_pilot_sievetl_tuning.R](run_PLC_CR_spline_hoffman2_pilot_sievetl_tuning.R)

This pilot step searches over:
- `c_sieve_grid = c(0.5, 1)`
- `c_zeta_grid = c(0, 0.05, 0.1, 0.2)`
- `c_eta_grid = c(0, 0.05, 0.1, 0.2)`

For each pilot seed, it writes files such as:
- `pilot_bic_tables_sievetl_seedXXXX.csv`
- `pilot_selected_sievetl_tuning_seedXXXX.csv`

### Step 2. Stage-1 pilot summary

Run:
- [submit_sievetl_pilot_tuning_summary.sh](submit_sievetl_pilot_tuning_summary.sh)

The underlying R code is:
- [summarize_spline_sievetl_pilot_tuning.R](summarize_spline_sievetl_pilot_tuning.R)

This step reads the pilot outputs from step 1 and aggregates them into:
- `final_selected_tuning.csv`

This is the main SieveTL tuning file used by the later steps.

### Step 3. Stage-2 pilot tuning for the by-column inverse-Hessian penalties

Run:
- [submit_spline_invh_pilot_by_column_array.sh](submit_spline_invh_pilot_by_column_array.sh)

This script launches Hoffman2 array jobs for the by-column inverse-Hessian pilot replications.

By default it also uses:
- `NPILOT = 20`

The underlying R code is:
- [run_PLC_CR_spline_hoffman2_pilot_invh_tuning_by_column.R](run_PLC_CR_spline_hoffman2_pilot_invh_tuning_by_column.R)

This step:
- reads `final_selected_tuning.csv` from step 2
- fits SieveTL on each pilot dataset using that stage-1 tuning
- tunes the column-specific inverse-Hessian penalties

For each column, it searches over:
- `hessian_c_grid = c(0, 0.1, 0.5, 1, 2, 4, 8)`

For each pilot seed, it writes files such as:
- `pilot_bic_tables_by_column_seedXXXX.csv`
- `pilot_bic_tables_by_column_seedXXXX.rds`
- `pilot_selected_c_by_column_seedXXXX.csv`

### Step 4. Stage-2 pilot summary

Run:
- [submit_spline_invh_pilot_by_column_summary.sh](submit_spline_invh_pilot_by_column_summary.sh)

The underlying R code is:
- [summarize_spline_invh_pilot_by_column.R](summarize_spline_invh_pilot_by_column.R)

This step reads the pilot outputs from step 3 and aggregates them into:
- `final_selected_c_by_column.csv`

This is the by-column inverse-Hessian tuning file used by the final jobs.

### Step 5. Final coverage / inference jobs

Run:
- [submit_coverage_plugin_wald.sh](submit_coverage_plugin_wald.sh)

The underlying R code is:
- [run_coverage_plugin_wald_hoffman2.R](run_coverage_plugin_wald_hoffman2.R)

This is the final job in the kept workflow.

It reads:
- `final_selected_tuning.csv`
- `final_selected_c_by_column.csv`

and writes raw `.rds` outputs for the coverage / inference results.

## What Each File Does

### Core R files

- [run_PLC_CR_spline_hoffman2.R](run_PLC_CR_spline_hoffman2.R)
  - shared core implementation used by the pilot scripts
  - defines the simulation model
  - defines SieveTL fitting code
  - defines spline construction helpers
  - defines score and Hessian helpers
  - defines inverse-Hessian helpers
  - defines tuning-file path helpers

- [coverage_plugin_wald_sievetl.R](coverage_plugin_wald_sievetl.R)
  - core inference engine for the final coverage jobs
  - computes the SieveTL estimator
  - computes the one-step estimator
  - computes pointwise intervals
  - computes SCB using the multiplier bootstrap
  - packages replication-level output objects

- [run_coverage_plugin_wald_hoffman2.R](run_coverage_plugin_wald_hoffman2.R)
  - final Hoffman2 runner for the coverage jobs
  - maps array tasks to design rows and replication chunks
  - calls `coverage_plugin_wald_sievetl.R`
  - writes raw `.rds` files

### Stage-1 pilot files

- [run_PLC_CR_spline_hoffman2_pilot_sievetl_tuning.R](run_PLC_CR_spline_hoffman2_pilot_sievetl_tuning.R)
  - stage-1 pilot runner
  - simulates one pilot dataset per task
  - searches the SieveTL tuning grid

- [summarize_spline_sievetl_pilot_tuning.R](summarize_spline_sievetl_pilot_tuning.R)
  - stage-1 pilot summary
  - combines the pilot outputs
  - writes `final_selected_tuning.csv`

### Stage-2 pilot files

- [run_PLC_CR_spline_hoffman2_pilot_invh_tuning_by_column.R](run_PLC_CR_spline_hoffman2_pilot_invh_tuning_by_column.R)
  - stage-2 pilot runner
  - reads the stage-1 tuning file
  - fits SieveTL on pilot data
  - tunes the by-column inverse-Hessian penalties

- [summarize_spline_invh_pilot_by_column.R](summarize_spline_invh_pilot_by_column.R)
  - stage-2 pilot summary
  - combines the by-column pilot outputs
  - writes `final_selected_c_by_column.csv`

### Submit scripts

- [submit_sievetl_pilot_tuning_array.sh](submit_sievetl_pilot_tuning_array.sh)
  - Hoffman2 array script for stage-1 pilot tuning

- [submit_sievetl_pilot_tuning_summary.sh](submit_sievetl_pilot_tuning_summary.sh)
  - Hoffman2 script for stage-1 pilot summary

- [submit_spline_invh_pilot_by_column_array.sh](submit_spline_invh_pilot_by_column_array.sh)
  - Hoffman2 array script for stage-2 pilot tuning

- [submit_spline_invh_pilot_by_column_summary.sh](submit_spline_invh_pilot_by_column_summary.sh)
  - Hoffman2 script for stage-2 pilot summary

- [submit_coverage_plugin_wald.sh](submit_coverage_plugin_wald.sh)
  - Hoffman2 array script for the final coverage / inference jobs

## Main Outputs

This workflow produces three types of outputs.

### 1. Stage-1 pilot tuning output
- `final_selected_tuning.csv`

### 2. Stage-2 pilot tuning output
- `final_selected_c_by_column.csv`

### 3. Final coverage / inference output
- raw per-task `.rds` files written by the final coverage jobs

## Notes

- This folder is focused on raw method outputs rather than postprocessing.
- The scripts assume a Hoffman2-style batch environment.
- The final kept workflow is:
  1. stage-1 pilot
  2. stage-1 summary
  3. stage-2 pilot
  4. stage-2 summary
  5. final coverage `.rds` jobs
