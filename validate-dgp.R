# validate-dgp.R -- cheap full-data checks of R/dgp-rct.R before any sim.
source("R/dgp-rct.R")
set.seed(42)
out <- list()
for (ds in c("pbc", "prostate")) for (hz in names(hazard_shapes)) for (cen in c(TRUE, FALSE)) {
  ref <- load_reference(ds)
  b <- calibrate_rct_censoring(ref, -0.4, hz, 0.3, centre = cen)
  d <- simulate_rct(ref, 20000, -0.4, hz, b, centre = cen)
  tr <- attr(d, "truth")
  tt <- quantile(d$time[d$status == 1], c(.25, .5, .75))
  km <- summary(survfit(Surv(time, status) ~ 1, d), times = tt)$surv
  true <- colMeans(rct_true_surv(tr$X, tr$trt, ref, -0.4, hz, tt, cen))
  cox_trt <- if (hz != "nph") unname(coef(coxph(Surv(time, status) ~ ., d))["trt"]) else NA
  out[[length(out) + 1]] <- data.frame(ds, hz, centre = cen, cens = round(1 - mean(d$status), 3),
    max_km_vs_true = round(max(abs(km - true)), 4), cox_trt_logHR = round(cox_trt, 2),
    median_T = signif(median(d$time), 3))
}
print(do.call(rbind, out), row.names = FALSE)
