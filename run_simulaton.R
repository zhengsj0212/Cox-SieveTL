# run_simulaton.R
#
# Run nsim simulations with spline basis.
# For each replicate, save:
#   - beta_hat  (target coefficient vector)
#   - gamma_hat (baseline coefficient vector)
#   - baseline  (baseline hazard and cumulative baseline hazard on a grid)

library(MASS)
library(foreach)
library(doParallel)
library(pracma)
library(splines)

source("R/simulate.R")
source("R/sieveTL.R")

nsim   <- 100
ncores <- 6
set.seed(123)

beta_kappa_list <- list(
  c(0.5, 2)
  # additional combinations can be added here
)

c1_vals    <- c(1.50, 4.55)
L_vals     <- c(5)
shift_vals <- c(TRUE)

n_target <- 500
n_source <- 2000

cl <- makeCluster(ncores)
registerDoParallel(cl)

dir.create("results", showWarnings = FALSE)
dir.create("results/spline", showWarnings = FALSE)

for (bk in beta_kappa_list) {
  beta2 <- bk[1]
  kappa <- bk[2]
  
  for (c1 in c1_vals) {
    for (L in L_vals) {
      for (covariate_shift in shift_vals) {
        
        filename_base <- sprintf(
          "beta%.1f_kappa%d_L%d_c%.2f_shift%s",
          beta2, kappa, L, c1,
          ifelse(covariate_shift, "TRUE", "FALSE")
        )
        
        file_res <- file.path(
          "results", "spline",
          paste0(filename_base, ".rds")
        )
        
        message(
          "Running: beta2=", beta2,
          ", kappa=", kappa,
          ", c1=", c1,
          ", L=", L,
          ", shift=", covariate_shift
        )
        
        # fixed basis size and penalties (you can tune these)
        p_fixed <- 4
        if (kappa == 2) {
          if (beta2 == 0.5) {
            lambda_zeta <- 0.05
            lambda_eta  <- 0.05
          } else {
            lambda_zeta <- 0.02
            lambda_eta  <- 0
          }
        } else {
          if (kappa %in% c(1, 3)) {
            lambda_eta  <- 0.02
          } else {
            lambda_eta  <- 0
          }
          lambda_zeta <- 0
        }
        
        results <- foreach(
          i = 1:nsim,
          .packages = c("MASS", "pracma", "foreach", "splines", "stats"),
          .export   = c(
            "simulate_once_all", "sievetl",
            "transfer_loss_fast", "debias_loss_fast",
            "transfer_grad_fast", "debias_grad_fast",
            "make_baseline_from_gamma"
          )
        ) %dopar% {
          
          tryCatch({
            
            # 1) simulate data
            D <- simulate_once_all(
              n_target = n_target, n_source = n_source,
              beta2_source = beta2, kappa_source = kappa,
              c1 = c1, covariate_shift = covariate_shift, L = L
            )
            
            # 2) fit spline-based sieveTL model
            fit <- sievetl(
              target = D$target, source = D$source,
              p = p_fixed, L = L,
              lambda_zeta = lambda_zeta, lambda_eta = lambda_eta
            )
            
            # 3) reconstruct spline basis (using pooled times) for baseline path
            combined_stimes <- c(D$target$stime, D$source[[1]]$stime)
            B0 <- splines::bs(
              combined_stimes,
              df = p_fixed,
              degree = 3,
              intercept = TRUE
            )
            at    <- attributes(B0)
            knots <- at$knots
            bkn   <- at$Boundary.knots
            deg   <- at$degree
            
            g_funcs <- lapply(seq_len(ncol(B0)), function(j) {
              function(t) {
                Bt <- splines::bs(
                  t,
                  degree = deg,
                  knots = knots,
                  Boundary.knots = bkn,
                  intercept = TRUE
                )
                as.numeric(Bt[, j])
              }
            })
            
            # 4) baseline cumulative hazard path on a fine grid
            t_cumh <- seq(
              max(bkn[1], 0.01),
              min(bkn[2], max(D$target$stime, na.rm = TRUE)),
              length.out = 100
            )
            
            bh_sieveTL <- make_baseline_from_gamma(
              fit$gamma.hat, g_funcs, t_cumh
            )
            
            list(
              beta_hat  = fit$beta.hat,
              gamma_hat = fit$gamma.hat,
              baseline  = bh_sieveTL
            )
            
          }, error = function(e) {
            message(paste("Task", i, "failed:", e$message))
            NULL
          })
        }
        
        saveRDS(results, file = file_res)
      }
    }
  }
}

stopCluster(cl)

