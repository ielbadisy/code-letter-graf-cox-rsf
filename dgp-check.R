# dgp-check.R -- Section 2 of the letter: the published prostate DGP of Graf et al. (2026).
# Input: data-ref/adv_prostate_ca/adv_prostate_ca.csv. Output: results/dgp_check.csv.
# Uncentred Cox refit, treatment log-HR shift at the covariate means, published scale (39.2)
# against Weibull fits, median simulated time under the uncentred DGP against the trial median.
# Run from the repo root: Rscript dgp-check.R
source("R/dgp-rct.R")
ref <- load_reference("prostate")
pc <- read.csv("data-ref/adv_prostate_ca/adv_prostate_ca.csv")
f <- Surv(survtime, cens) ~ rx + age + wt + sbp + dbp + sz + ap + hg + sg + pf + hx + bm + ekg +
  rx:age + rx:bm + rx:ap
cx <- coxph(f, pc); w0 <- survreg(Surv(survtime, cens) ~ 1, pc); w1 <- survreg(f, pc)
shift <- sum(ref$inter * colMeans(ref$X)[names(ref$inter)])
set.seed(42)
med <- sapply(c(dec = "dec", const = "const", inc = "inc"), function(hz) {
  X <- ref$X[sample.int(nrow(ref$X), 1e6, TRUE), ]
  median(draw_rct_times(X, rbinom(1e6, 1, 0.5), ref, -0.4, hz, centre = FALSE))
})
v <- c(n_patients = nrow(pc), max_abs_coef_diff = max(abs(coef(cx)[names(ref$beta)] - ref$beta)),
       refit_coef_rx = coef(cx)[["rx"]], trt_shift_at_means = shift, trt_logHR_at_means = -0.4 + shift,
       published_scale = ref$scale, weibull_aft_scale_null = exp(coef(w0))[[1]],
       weibull_ph_scale_null = exp(coef(w0))[[1]]^(1 / w0$scale), weibull_aft_scale_full = exp(coef(w1))[[1]],
       setNames(med, paste0("median_T_", names(med))),
       median_T_trial = summary(survfit(Surv(survtime, cens) ~ 1, pc))$table[["median"]])
write.csv(data.frame(quantity = names(v), value = v), "results/dgp_check.csv", row.names = FALSE)
