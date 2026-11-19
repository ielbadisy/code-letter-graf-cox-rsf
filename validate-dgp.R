# validate-dgp.R -- cheap check of R/dgp-rct.R on 20 000 patients before the simulation:
# censoring rate, max |KM - mean true S| at the event-time quartiles, Cox treatment log-HR
# (PH hazards), median time; centred and uncentred. Output: console.
# Run from the repo root: Rscript validate-dgp.R
source("R/dgp-rct.R")
set.seed(42)
chk <- function(ds, hz, cen) {
  ref <- load_reference(ds)
  d <- simulate_rct(ref, 20000, -0.4, hz, calibrate_rct_censoring(ref, -0.4, hz, 0.3, centre = cen), cen)
  tr <- attr(d, "truth"); tt <- quantile(d$time[d$status == 1], c(.25, .5, .75))
  km <- summary(survfit(Surv(time, status) ~ 1, d), times = tt)$surv
  data.frame(ds, hz, centre = cen, cens = round(1 - mean(d$status), 3),
    max_km_vs_true = round(max(abs(km - colMeans(rct_true_surv(tr$X, tr$trt, ref, -0.4, hz, tt, cen)))), 4),
    cox_trt_logHR = if (hz == "nph") NA else round(coef(coxph(Surv(time, status) ~ ., d))[["trt"]], 2),
    median_T = signif(median(d$time), 3))
}
g <- expand.grid(cen = c(TRUE, FALSE), hz = names(hazard_shapes), ds = c("pbc", "prostate"),
                 stringsAsFactors = FALSE)
print(do.call(rbind, Map(chk, g$ds, g$hz, g$cen)), row.names = FALSE)
