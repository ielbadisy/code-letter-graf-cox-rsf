# sim-run.R -- simulation for Sections 3-4 of the letter on Graf, Todd & Baksh (2026).
# Same reference trials, Weibull event times, U[0, b] censoring, beta_trt = -0.4, training
# N in {100, 200, 400}, independent test sets of N = 500, 48 scenarios.
# Models: Graf-style Cox, Cox stratified by arm with the true interactions, two RSFs.
# Metrics: Harrell C and IBS on one prediction scale, mean |S_hat - S_true| (true S known).
# Output: results/sim_raw_<reps>reps.csv.gz. Run from the repo root: Rscript sim-run.R 200
suppressMessages({ library(randomForestSRC); library(ranger); library(mcstatsim) })
source("R/dgp-rct.R")
# one thread per learner: parallelism comes from the runsim forks only
options(rf.cores = 1L, mc.cores = 1L); Sys.setenv(OMP_NUM_THREADS = "1")

BETA_TRT <- -0.4; N_TEST <- 500L
MODELS <- c("cox_main", "cox_spec", "rsf_logrank", "rsf_extratrees")
REFS <- list(pbc = load_reference("pbc"), prostate = load_reference("prostate"))

grid <- expand.grid(dataset = c("pbc", "prostate"), hazard = c("dec", "const", "inc", "nph"),
                    n = c(100L, 200L, 400L), cens = c(0.3, 0.6), stringsAsFactors = FALSE)
grid$cens_b <- mapply(function(ds, hz, cs) calibrate_rct_censoring(REFS[[ds]], BETA_TRT, hz, cs),
                      grid$dataset, grid$hazard, grid$cens)

# evaluation times: 10th..90th percentiles of the marginal true event-time distribution
TIMES <- list()
for (ds in names(REFS)) for (hz in names(hazard_shapes)) {
  set.seed(7); ref <- REFS[[ds]]
  idx <- sample.int(nrow(ref$X), 1e5, TRUE); trt <- rbinom(1e5, 1, 0.5)
  TIMES[[paste(ds, hz)]] <- unname(quantile(draw_rct_times(ref$X[idx, ], trt, ref, BETA_TRT, hz),
                                            seq(0.1, 0.9, 0.1)))
}

# right-continuous step function S (rows x tgrid) read at `times`
step_at <- function(tgrid, S, times) cbind(1, S)[, findInterval(times, tgrid) + 1, drop = FALSE]
# with strata, summary()$surv is a flat vector, curve after curve
cox_surv <- function(f, newdata, times)
  t(matrix(summary(survfit(f, newdata = newdata), times = times, extend = TRUE)$surv, length(times)))

fit_predict <- function(model, train, test, times, ref) {
  covs <- names(ref$X)
  switch(model,
    cox_main = cox_surv(coxph(reformulate(c("trt", covs), "Surv(time, status)"), train), test, times),
    cox_spec = {
      # interactions as product columns: survfit.coxph refuses trt:x when trt is only a stratum
      ix <- sprintf("trtx_%s", names(ref$inter))
      for (v in names(ref$inter)) { train[[paste0("trtx_", v)]] <- train$trt * train[[v]]
                                    test[[paste0("trtx_", v)]] <- test$trt * test[[v]] }
      f <- coxph(reformulate(c(covs, ix, "strata(trt)"), "Surv(time, status)"), train)
      S <- matrix(NA_real_, nrow(test), length(times))
      for (a in 0:1) { i <- test$trt == a; S[i, ] <- cox_surv(f, test[i, ], times) }
      S
    },
    rsf_logrank = {
      p <- predict(rfsrc(Surv(time, status) ~ ., train, ntree = 500, splitrule = "logrank"), test)
      step_at(p$time.interest, p$survival, times)
    },
    rsf_extratrees = {
      p <- predict(ranger(Surv(time, status) ~ ., train, num.trees = 500, splitrule = "extratrees",
                          num.threads = 1), test)
      step_at(p$unique.death.times, p$survival, times)
    })
}

# IPCW Brier score averaged over `times`, KM censoring distribution from the test set
ibs <- function(time, status, S, times) {
  G <- survfit(Surv(time, 1 - status) ~ 1); Gf <- stepfun(G$time, c(1, G$surv))
  Gt <- pmax(Gf(time - 1e-8), 1e-3)
  mean(vapply(seq_along(times), function(k) mean((time <= times[k] & status == 1) / Gt * S[, k]^2 +
    (time > times[k]) / max(Gf(times[k]), 1e-3) * (1 - S[, k])^2), 0))
}

sim_one <- function(dataset, hazard, n, cens, cens_b) {
  ref <- REFS[[dataset]]; times <- TIMES[[paste(dataset, hazard)]]
  train <- simulate_rct(ref, n, BETA_TRT, hazard, cens_b)
  test <- simulate_rct(ref, N_TEST, BETA_TRT, hazard, cens_b)
  tr <- attr(test, "truth"); y <- Surv(test$time, test$status)
  S_true <- rct_true_surv(tr$X, tr$trt, ref, BETA_TRT, hazard, times)
  do.call(rbind, lapply(MODELS, function(m) {
    S <- fit_predict(m, train, test, times, ref)
    risk <- 1 - S[, 5]   # predicted risk at the median event time, common to all models
    data.frame(model = m, metric = c("harrell_c", "ibs", "mae_true"),
               value = c(concordance(y ~ risk, reverse = TRUE)$concordance,
                         ibs(test$time, test$status, S, times), mean(abs(S - S_true))))
  }))
}

if (sys.nframe() == 0L) {
  reps <- as.integer(commandArgs(TRUE)[1])
  raw <- runsim(n = reps, grid_params = grid, sim_func = sim_one, seed = 2026,
                num_cores = as.integer(Sys.getenv("NCORES", 14)), backend = "fork")
  write.csv(raw, gzfile(sprintf("results/sim_raw_%dreps.csv.gz", reps)), row.names = FALSE)
}
