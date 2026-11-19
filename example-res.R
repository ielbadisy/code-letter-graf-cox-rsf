# example-res.R -- Section 1 of the letter: Cox C index on the two real trials, against
# Graf, Todd & Baksh (2026) Table 3 (0.632+ bootstrap: PBC 0.776, prostate 0.521).
# Same full models (main effects + treatment; prostate adds rx:age, rx:bm, rx:ap), Harrell C,
# apparent / out-of-bag / 0.632+ (Efron & Tibshirani 1997, no-information C = 0.5), then
# backward AIC selection repeated inside every bootstrap sample.
# Outputs: results/realdata_cox_cindex.csv, results/realdata_cox_cindex_aicselection.csv.
# Run from the repo root: Rscript example-res.R
suppressMessages(library(survival))

pc <- transform(read.csv("data-ref/adv_prostate_ca/adv_prostate_ca.csv"), time = survtime, status = cens)
f_pc <- Surv(time, status) ~ rx + age + wt + sbp + dbp + sz + ap + hg + sg + pf + hx + bm + ekg +
  rx:age + rx:bm + rx:ap
pb <- transform(subset(survival::pbc, !is.na(trt)), status = as.integer(status > 0))
for (v in c("chol", "copper", "trig", "platelet")) pb[[v]][is.na(pb[[v]])] <- mean(pb[[v]], na.rm = TRUE)
f_pb <- Surv(time, status) ~ trt + age + sex + ascites + hepato + spiders + edema + bili + chol +
  albumin + copper + alk.phos + ast + trig + platelet + protime + stage

cidx <- function(fit, d) concordance(Surv(d$time, d$status) ~ predict(fit, newdata = d),
                                     reverse = TRUE)$concordance
full <- function(f, d) do.call(coxph, list(f, data = d))   # data embedded so step() can refit
# backward AIC, treatment (and its interactions) always kept
aic <- function(lower) function(f, d) step(full(f, d), scope = list(lower = lower, upper = f), trace = 0)

boot632p <- function(f, d, fitter, B) {
  set.seed(1)
  app <- cidx(fitter(f, d), d)
  oob <- mean(replicate(B, { i <- sample.int(nrow(d), replace = TRUE)
                             cidx(fitter(f, d[i, ]), d[-unique(i), ]) }))
  R <- if (app > oob) min(1, (app - oob) / (app - 0.5)) else 0
  w <- 0.632 / (1 - 0.368 * R)
  c(apparent = app, oob = oob, c632p = (1 - w) * app + w * oob)
}

res <- rbind(prostate = boot632p(f_pc, pc, full, 500), pbc = boot632p(f_pb, pb, full, 500))
write.csv(round(res, 3), "results/realdata_cox_cindex.csv")
res_sel <- rbind(prostate = boot632p(f_pc, pc, aic(~ rx + rx:age + rx:bm + rx:ap), 200),
                 pbc = boot632p(f_pb, pb, aic(~ trt), 200))
write.csv(round(res_sel, 3), "results/realdata_cox_cindex_aicselection.csv")
