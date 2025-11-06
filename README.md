# Cox-SieveTL: Semiparametric Transfer Learning for Cox Models via Sieve Maximum Likelihood

This repository contains simulation code for the method proposed in

> **"Cox-SieveTL: Semiparametric Transfer Learning for Cox Models via Sieve Maximum Likelihood"**

The code implements a sieve-based transfer learning procedure (`sievetl`) for Cox proportional hazards models, together with data generation and simple utilities for reconstructing the baseline hazard.

------------------------------------------------------------------------

## Repository structure

``` text
GitHub_proj/
  R/
    simulate.R       # simulation of target + source survival data
    sieveTL.R        # core Cox-SieveTL implementation
  run_simulaton.R    # main script to run simulation studies
```

-   **`R/simulate.R`**
    -   Provides `simulate_once_all()` to generate a target dataset and one source dataset under various covariate shift scenarios.
-   **`R/sieveTL.R`**
    -   Core implementation of the Cox-SieveTL method:
        -   `sievetl()` -- fits the sieve-based transfer learning Cox model and returns estimated coefficients for covariates and baseline hazard basis.
        -   `make_baseline_from_gamma()` -- given the estimated baseline coefficients and basis functions, reconstructs the baseline hazard and cumulative baseline hazard on a time grid.
    -   Also contains internal helpers:
        -   `transfer_loss_fast()`, `transfer_grad_fast()`
        -   `debias_loss_fast()`, `debias_grad_fast()`
-   **`run_simulaton.R`**
    -   Example simulation driver:
        -   loops over combinations of `(beta2, kappa, c1, L, covariate_shift)`
        -   runs `nsim` Monte Carlo replications in parallel
        -   saves results (estimated coefficients and baseline hazard paths) as `.rds` files into `results/spline/`.

------------------------------------------------------------------------

## Requirements

This code is written in **R** and relies on the following packages:

-   `MASS`
-   `pracma`
-   `splines`
-   `foreach`
-   `doParallel`
-   `stats` (base R)

You can install missing packages via:

``` r
install.packages(c("MASS", "pracma", "splines", "foreach", "doParallel"))
```

------------------------------------------------------------------------

## How to run the simulations

From R:

``` r
setwd("path/to/GitHub_proj")  # set to the root of this repository
source("run_simulaton.R")
```

The script will:

1.  Set up a parallel cluster.
2.  Repeatedly simulate target and source data via `simulate_once_all()`.
3.  Fit the Cox-SieveTL model using `sievetl()`.
4.  Reconstruct the baseline cumulative hazard on a grid using `make_baseline_from_gamma()`.
5.  Save all results into `results/spline/*.rds`.

By default:

-   `nsim = 100` Monte Carlo replications,
-   target sample size `n_target = 500`,
-   source sample size `n_source = 2000`,
-   spline basis dimension `p_fixed = 4`,
-   number of covariates `L = 5`.

You can edit these values directly in `run_simulaton.R`.

------------------------------------------------------------------------

## Example: basic usage of `sievetl()`

Below is a minimal example showing how to:

1.  Simulate one target--source pair,
2.  Fit the Cox-SieveTL model,
3.  Reconstruct the baseline cumulative hazard on a grid.

``` r
library(MASS)
library(pracma)
library(splines)

source("R/simulate.R")
source("R/sieveTL.R")

set.seed(1)

D <- simulate_once_all(
  n_target       = 200,
  n_source       = 500,
  beta2_source   = 0.5,
  kappa_source   = 2,
  c1             = 4.55,
  covariate_shift = TRUE,
  L              = 5
)

target <- D$target
source_list <- D$source

fit <- sievetl(
  target      = target,
  source      = source_list,
  p           = 4,
  L           = 5,
  lambda_zeta = 0.05,
  lambda_eta  = 0.05
)

fit$beta.hat
fit$gamma.hat
```

### Reconstructing the baseline hazard and cumulative hazard

``` r
all_stime <- c(target$stime, source_list[[1]]$stime)
B0 <- splines::bs(all_stime, df = 4, degree = 3, intercept = TRUE)
at    <- attributes(B0)
knots <- at$knots
bkn   <- at$Boundary.knots
deg   <- at$degree

g_funcs <- lapply(seq_len(ncol(B0)), function(j) {
  function(t) {
    Bt <- splines::bs(t, degree = deg, knots = knots, Boundary.knots = bkn, intercept = TRUE)
    as.numeric(Bt[, j])
  }
})

t_grid <- seq(max(bkn[1], 0.01), min(bkn[2], max(target$stime, na.rm = TRUE)), length.out = 100)

baseline_df <- make_baseline_from_gamma(fit$gamma.hat, g_funcs, t_grid)

plot(baseline_df$time, baseline_df$H0, type = "l", xlab = "Time", ylab = "Cumulative baseline hazard")
```
