# dgp-check.R -- Section 2 of the letter: the published prostate DGP of Graf et al. (2026).
# Input: data-ref/adv_prostate_ca/adv_prostate_ca.csv. Output: results/dgp_check.csv.
# Refits the uncentred Cox model, computes the treatment log-HR shift at the covariate
# means, and the median simulated survival time with the published scale (39.2).
source("R/dgp-rct.R")
ref <- load_reference("prostate")
pc <- read.csv("data-ref/adv_prostate_ca/adv_prostate_ca.csv")
f <- coxph(Surv(survtime, cens) ~ rx + age + wt + sbp + dbp + sz + ap + hg + sg + pf + hx + bm + ekg +
  rx:age + rx:bm + rx:ap, pc)
shift <- sum(ref$inter * colMeans(ref$X)[names(ref$inter)])
set.seed(42)
med <- sapply(c(dec = "dec", const = "const", inc = "inc"), function(hz) {
  n <- 2e5; X <- ref$X[sample.int(nrow(ref$X), n, TRUE), ]
  median(draw_rct_times(X, rbinom(n, 1, 0.5), ref, -0.4, hz, centre = FALSE))
})
km <- unname(summary(survfit(Surv(survtime, cens) ~ 1, pc))$table["median"])
res <- data.frame(quantity = c("refit_coef_rx", "trt_shift_at_means", "trt_logHR_at_means",
                               paste0("median_T_", names(med)), "median_T_trial"),
                  value = round(c(coef(f)[["rx"]], shift, -0.4 + shift, med, km), 3))
print(res, row.names = FALSE)
print(round(cbind(refit = coef(f)[names(ref$beta)], published = ref$beta), 4))
write.csv(res, "results/dgp_check.csv", row.names = FALSE)
