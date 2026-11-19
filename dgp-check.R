# dgp-check.R -- Section 2 of the letter: the published prostate DGP of Graf et al. (2026).
# Input: data-ref/adv_prostate_ca/adv_prostate_ca.csv. Output: results/dgp_check.csv.
# Refits the uncentred Cox model, computes the treatment log-HR shift at the covariate
# means, compares the published scale (39.2) with Weibull fits, and gives the median
# simulated survival time under the published (uncentred) DGP against the trial median.
source("R/dgp-rct.R")
ref <- load_reference("prostate")
pc <- read.csv("data-ref/adv_prostate_ca/adv_prostate_ca.csv")
f <- Surv(survtime, cens) ~ rx + age + wt + sbp + dbp + sz + ap + hg + sg + pf + hx + bm + ekg +
  rx:age + rx:bm + rx:ap
cx <- coxph(f, pc)
shift <- sum(ref$inter * colMeans(ref$X)[names(ref$inter)])
# Weibull fits: AFT scale exp(intercept), and the PH-form scale of the null model
w0 <- survreg(Surv(survtime, cens) ~ 1, pc); w1 <- survreg(f, pc)
set.seed(42)
med <- sapply(c(dec = "dec", const = "const", inc = "inc"), function(hz) {
  n <- 1e6; X <- ref$X[sample.int(nrow(ref$X), n, TRUE), ]
  median(draw_rct_times(X, rbinom(n, 1, 0.5), ref, -0.4, hz, centre = FALSE))
})
km <- unname(summary(survfit(Surv(survtime, cens) ~ 1, pc))$table["median"])
res <- data.frame(
  quantity = c("n_patients", "max_abs_coef_diff", "refit_coef_rx", "trt_shift_at_means",
               "trt_logHR_at_means", "published_scale", "weibull_aft_scale_null",
               "weibull_ph_scale_null", "weibull_aft_scale_full", paste0("median_T_", names(med)),
               "median_T_trial"),
  value = c(nrow(pc), max(abs(coef(cx)[names(ref$beta)] - ref$beta)), coef(cx)[["rx"]], shift,
            -0.4 + shift, ref$scale, exp(coef(w0)), exp(coef(w0))^(1 / w0$scale),
            exp(coef(w1)[[1]]), med, km))
print(res, row.names = FALSE)
write.csv(res, "results/dgp_check.csv", row.names = FALSE)
