# sim-run.R -- reanalysis of Graf, Todd & Baksh (2026) for the letter.
# Same reference trials, Weibull PH event times, U[0, b] censoring, beta_trt = -0.4,
# training N in {100, 200, 400}, independent test sets of N = 500.
# Adds: (i) a correctly specified Cox model (stratified by arm, with the true
# treatment-covariate interactions), (ii) calibration against the TRUE S(t|x),
# (iii) the IPCW Brier score and C computed on one common prediction scale.
# Run from the repo root: Rscript sim-run.R

suppressMessages({ library(survival); library(randomForestSRC); library(ranger); library(mcstatsim) })
source("R/dgp-rct.R")

BETA_TRT <- -0.4
# one thread per learner: parallelism comes from the runsim forks only
options(rf.cores = 1L, mc.cores = 1L)
Sys.setenv(OMP_NUM_THREADS = "1")
N_TEST <- 500L
REFS <- list(pbc = load_reference("pbc"), prostate = load_reference("prostate"))

grid <- expand.grid(dataset = c("pbc", "prostate"), hazard = c("dec", "const", "inc", "nph"),
                    n = c(100L, 200L, 400L), cens = c(0.3, 0.6), stringsAsFactors = FALSE)
grid$cens_b <- mapply(function(ds, hz, cs) calibrate_rct_censoring(REFS[[ds]], BETA_TRT, hz, cs),
                      grid$dataset, grid$hazard, grid$cens)

# evaluation times: 10th..90th percentiles of the marginal true event-time distribution
eval_times <- function(ref, hazard) {
  set.seed(7)
  idx <- sample.int(nrow(ref$X), 1e5, TRUE); trt <- rbinom(1e5, 1, 0.5)
  unname(quantile(draw_rct_times(ref$X[idx, ], trt, ref, BETA_TRT, hazard), seq(0.1, 0.9, 0.1)))
}
TIMES <- list()
for (ds in names(REFS)) for (hz in names(hazard_shapes)) TIMES[[paste(ds, hz)]] <- eval_times(REFS[[ds]], hz)

step_at <- function(tgrid, S, times) {
  # S: n x length(tgrid) step function, evaluated at `times` (right-continuous)
  j <- findInterval(times, tgrid)
  out <- matrix(1, nrow(S), length(times))
  out[, j > 0] <- S[, j[j > 0], drop = FALSE]
  out
}

fit_predict <- function(model, train, test, times, ref) {
  covs <- names(ref$X)
  main <- paste(c("trt", covs), collapse = " + ")

  switch(model,
    cox_main = {   # Graf-style: main effects, one PH model
      f <- coxph(as.formula(paste("Surv(time, status) ~", main)), train, x = TRUE)
      sf <- survfit(f, newdata = test)
      t(summary(sf, times = times, extend = TRUE)$surv)
    },
    cox_spec = {   # correctly specified: strata(trt) frees the arm-specific baseline, true interactions
      # interactions as explicit product columns: survfit.coxph refuses trt:x
      # when trt itself only enters through strata()
      ix <- names(ref$inter)
      for (v in ix) { train[[paste0("trtx_", v)]] <- train$trt * train[[v]]
                      test[[paste0("trtx_", v)]] <- test$trt * test[[v]] }
      rhs <- paste(c(covs, if (length(ix)) paste0("trtx_", ix)), collapse = " + ")
      f <- coxph(as.formula(paste("Surv(time, status) ~", rhs, "+ strata(trt)")), train)
      S <- matrix(NA_real_, nrow(test), length(times))
      for (a in 0:1) {
        i <- which(test$trt == a)
        if (length(i)) {
          sf <- survfit(f, newdata = test[i, ])
          # with strata, summary()$surv is a flat vector, curve after curve
          S[i, ] <- t(matrix(summary(sf, times = times, extend = TRUE)$surv, nrow = length(times)))
        }
      }
      S
    },
    rsf_logrank = {
      f <- rfsrc(Surv(time, status) ~ ., data = train, ntree = 500, splitrule = "logrank")
      p <- predict(f, newdata = test)
      step_at(p$time.interest, p$survival, times)
    },
    rsf_extratrees = {
      f <- ranger(Surv(time, status) ~ ., data = train, num.trees = 500, splitrule = "extratrees", num.threads = 1)
      p <- predict(f, data = test)
      step_at(p$unique.death.times, p$survival, times)
    })
}

# IPCW Brier score at each time, KM estimate of the censoring distribution from the test set
ipcw_brier <- function(time, status, S, times) {
  G <- survfit(Surv(time, 1 - status) ~ 1)
  Gf <- stepfun(G$time, c(1, G$surv))
  Gt <- pmax(Gf(time - 1e-8), 1e-3)
  sapply(seq_along(times), function(k) {
    t <- times[k]
    w_ev <- (time <= t & status == 1) / Gt
    w_al <- (time > t) / max(Gf(t), 1e-3)
    mean(w_ev * S[, k]^2 + w_al * (1 - S[, k])^2)
  })
}

models <- c("cox_main", "cox_spec", "rsf_logrank", "rsf_extratrees")

sim_one <- function(dataset, hazard, n, cens, cens_b) {
  ref <- REFS[[dataset]]; times <- TIMES[[paste(dataset, hazard)]]
  train <- simulate_rct(ref, n, BETA_TRT, hazard, cens_b)
  test <- simulate_rct(ref, N_TEST, BETA_TRT, hazard, cens_b)
  tr <- attr(test, "truth")
  S_true <- rct_true_surv(tr$X, tr$trt, ref, BETA_TRT, hazard, times)
  y <- Surv(test$time, test$status)
  do.call(rbind, lapply(models, function(m) {
    S <- tryCatch(fit_predict(m, train, test, times, ref), error = function(e) NULL)
    if (is.null(S) || anyNA(S)) return(data.frame(model = m, metric = c("harrell_c", "uno_c", "ibs", "mae_true", "mae_true_med"), value = NA_real_))
    risk <- 1 - S[, 5]   # predicted risk at the median event time, common to all models
    data.frame(model = m,
      metric = c("harrell_c", "uno_c", "ibs", "mae_true", "mae_true_med"),
      value = c(
        concordance(y ~ risk, reverse = TRUE)$concordance,
        concordance(y ~ risk, reverse = TRUE, timewt = "n/G2", ymax = times[9])$concordance,
        mean(ipcw_brier(test$time, test$status, S, times)),
        mean(abs(S - S_true)),           # mean |S_hat - S_true| over patients and the 9 times
        mean(abs(S[, 5] - S_true[, 5]))  # at the median event time
      ))
  }))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  reps <- if (length(args)) as.integer(args[1]) else 200L
  dir.create("results", showWarnings = FALSE)
  saveRDS(grid, "results/grid.rds")
  raw <- runsim(n = reps, grid_params = grid, sim_func = sim_one, seed = 2026,
                num_cores = as.integer(Sys.getenv("NCORES", 14)), backend = "fork", on_error = "warn")
  write.csv(raw, gzfile(sprintf("results/sim_raw_%dreps.csv.gz", reps)), row.names = FALSE)
}
