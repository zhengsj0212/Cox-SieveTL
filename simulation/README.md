# SieveTL Simulation Workflow

This folder contains the cleaned main simulation workflow for the SieveTL method.

It keeps only the code needed to:
- tune the main SieveTL penalties from pilot replications
- tune the by-column inverse-Hessian penalties from pilot replications
- run the final coverage / inference jobs
- save raw `.rds` outputs

It does **not** include:
- comparison methods
- plotting code
- summary CSV generation
- LaTeX table generation
- filtered postprocessing

## Procedure

Run the code in this order.

### 1. Tune the main SieveTL penalties
Run:
- [submit_sievetl_pilot_tuning_array.sh](submit_sievetl_pilot_tuning_array.sh)

This launches `NPILOT=20` pilot replications by default and searches over:
- `c_sieve_grid = c(0.5, 1)`
- `c_zeta_grid = c(0, 0.05, 0.1, 0.2)`
- `c_eta_grid = c(0, 0.05, 0.1, 0.2)`

Per-pilot files written:
- `pilot_bic_tables_sievetl_seedXXXX.csv`
- `pilot_selected_sievetl_tuning_seedXXXX.csv`

### 2. Summarize the main SieveTL pilot tuning
Run:
- [submit_sievetl_pilot_tuning_summary.sh](submit_sievetl_pilot_tuning_summary.sh)

This aggregates the pilot outputs and writes:
- `final_selected_tuning.csv`

This file is the stage-1 tuning file used by later steps.

### 3. Tune the by-column inverse-Hessian penalties
Run:
- [submit_spline_invh_pilot_by_column_array.sh](submit_spline_invh_pilot_by_column_array.sh)

This also uses `NPILOT=20` by default.

For each column of the transfer inverse-Hessian correction, it searches over:
- `hessian_c_grid = c(0, 0.1, 0.5, 1, 2, 4, 8)`

Per-pilot files written:
- `pilot_bic_tables_by_column_seedXXXX.csv`
- `pilot_bic_tables_by_column_seedXXXX.rds`
- `pilot_selected_c_by_column_seedXXXX.csv`

### 4. Summarize the by-column inverse-Hessian pilot tuning
Run:
- [submit_spline_invh_pilot_by_column_summary.sh](submit_spline_invh_pilot_by_column_summary.sh)

This aggregates the pilot outputs and writes:
- `final_selected_c_by_column.csv`

This file is the stage-2 tuning file used by the final jobs.

### 5. Run the main coverage / inference jobs
Run:
- [submit_coverage_plugin_wald.sh](submit_coverage_plugin_wald.sh)

This final step reads:
- `final_selected_tuning.csv`
- `final_selected_c_by_column.csv`

and writes raw `.rds` outputs for the coverage / inference results.

## What Each File Does

### R files

- [run_PLC_CR_spline_hoffman2.R](run_PLC_CR_spline_hoffman2.R)
  - shared core implementation
  - defines the simulation model, SieveTL fitting code, spline construction, score/Hessian helpers, inverse-Hessian helpers, and tuning-file path helpers
  - the pilot scripts source this file

- [run_PLC_CR_spline_hoffman2_pilot_sievetl_tuning.R](run_PLC_CR_spline_hoffman2_pilot_sievetl_tuning.R)
  - stage-1 pilot runner
  - simulates one pilot dataset per task and selects the main SieveTL tuning parameters

- [summarize_spline_sievetl_pilot_tuning.R](summarize_spline_sievetl_pilot_tuning.R)
  - stage-1 pilot summary
  - combines the pilot outputs and writes `final_selected_tuning.csv`

- [run_PLC_CR_spline_hoffman2_pilot_invh_tuning_by_column.R](run_PLC_CR_spline_hoffman2_pilot_invh_tuning_by_column.R)
  - stage-2 pilot runner
  - uses the stage-1 tuning file, fits SieveTL on pilot data, and tunes the by-column inverse-Hessian penalties

- [summarize_spline_invh_pilot_by_column.R](summarize_spline_invh_pilot_by_column.R)
  - stage-2 pilot summary
  - combines the columnwise pilot outputs and writes `final_selected_c_by_column.csv`

- [coverage_plugin_wald_sievetl.R](coverage_plugin_wald_sievetl.R)
  - core inference engine for the final coverage jobs
  - computes the SieveTL and one-step estimators
  - computes pointwise intervals
  - computes SCB using the multiplier bootstrap
  - packages replication-level records and output objects

- [run_coverage_plugin_wald_hoffman2.R](run_coverage_plugin_wald_hoffman2.R)
  - final Hoffman2 runner for the coverage jobs
  - maps array tasks to design rows and replication chunks
  - calls the inference engine in `coverage_plugin_wald_sievetl.R`
  - writes raw `.rds` outputs

### Bash files

- [submit_sievetl_pilot_tuning_array.sh](submit_sievetl_pilot_tuning_array.sh)
  - Hoffman2 array submission script for stage-1 pilot tuning

- [submit_sievetl_pilot_tuning_summary.sh](submit_sievetl_pilot_tuning_summary.sh)
  - Hoffman2 submission script for stage-1 pilot summary

- [submit_spline_invh_pilot_by_column_array.sh](submit_spline_invh_pilot_by_column_array.sh)
  - Hoffman2 array submission script for stage-2 inverse-Hessian pilot tuning

- [submit_spline_invh_pilot_by_column_summary.sh](submit_spline_invh_pilot_by_column_summary.sh)
  - Hoffman2 submission script for stage-2 pilot summary

- [submit_coverage_plugin_wald.sh](submit_coverage_plugin_wald.sh)
  - Hoffman2 array submission script for the final coverage / inference jobs

## Main Outputs

The workflow produces three kinds of outputs.

- Stage-1 pilot tuning:
  - `final_selected_tuning.csv`

- Stage-2 pilot tuning:
  - `final_selected_c_by_column.csv`

- Final coverage / inference results:
  - raw per-task `.rds` files

## Notes

- This export is focused on raw method outputs rather than postprocessing.
- The remaining scripts assume a Hoffman2-style batch environment.
