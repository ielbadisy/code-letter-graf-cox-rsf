# example-res.R -- real-data check of the Cox C index reported by
# Graf, Todd & Baksh (2026), Tables 2-3 (0.632+ bootstrap): PBC 0.776, prostate 0.521.
# Same full models (all main effects + treatment; prostate adds rx:age, rx:bm, rx:ap),
# Harrell C, apparent / 10-fold CV / 0.632+ (Efron & Tibshirani 1997, noinfo = 0.5).
# Run from the repo root: Rscript example-res.R

suppressMessages(library(survival))

pc <- read.csv("data-ref/adv_prostate_ca/adv_prostate_ca.csv")
pc$time <- pc$survtime; pc$status <- pc$cens
f_pc <- Surv(time, status) ~ rx + age + wt + sbp + dbp + sz + ap + hg + sg + pf + hx + bm + ekg +
  rx:age + rx:bm + rx:ap

pb <- subset(survival::pbc, !is.na(trt))
pb$status <- as.integer(pb$status > 0)
for (v in c("chol", "copper", "trig", "platelet")) pb[[v]][is.na(pb[[v]])] <- mean(pb[[v]], na.rm = TRUE)
f_pb <- Surv(time, status) ~ trt + age + sex + ascites + hepato + spiders + edema + bili + chol +
  albumin + copper + alk.phos + ast + trig + platelet + protime + stage

cidx <- function(fit, newdata) {
  lp <- predict(fit, newdata = newdata)
  concordance(Surv(newdata$time, newdata$status) ~ lp, reverse = TRUE)$concordance
}

boot632p <- function(formula, data, B = 500, seed = 1) {
  set.seed(seed)
  fit <- coxph(formula, data)
  app <- cidx(fit, data)
  oob <- replicate(B, {
    i <- sample.int(nrow(data), replace = TRUE); o <- setdiff(seq_len(nrow(data)), i)
    tryCatch(cidx(coxph(formula, data[i, ]), data[o, ]), error = function(e) NA)
  })
  oob <- mean(oob, na.rm = TRUE)
  # 0.632+ for a "higher is better" measure: overfitting = app - oob, noinfo = 0.5
  R <- if (app > oob) min(1, (app - oob) / (app - 0.5)) else 0
  w <- 0.632 / (1 - 0.368 * R)
  c(apparent = app, oob = oob, c632p = (1 - w) * app + w * oob)
}

cv10 <- function(formula, data, seed = 1) {
  set.seed(seed); k <- sample(rep(1:10, length.out = nrow(data))); lp <- numeric(nrow(data))
  for (j in 1:10) lp[k == j] <- predict(coxph(formula, data[k != j, ]), data[k == j, ])
  concordance(Surv(data$time, data$status) ~ lp, reverse = TRUE)$concordance
}

res <- rbind(
  prostate = c(boot632p(f_pc, pc), cv10 = cv10(f_pc, pc), graf_reported = 0.521),
  pbc      = c(boot632p(f_pb, pb), cv10 = cv10(f_pb, pb), graf_reported = 0.776)
)
print(round(res, 3))
dir.create("results", showWarnings = FALSE)
write.csv(round(res, 3), "results/realdata_cox_cindex.csv")

# sensitivity: backward AIC selection repeated inside every bootstrap sample
# (Graf et al. select variables by AIC/BIC before fitting the Cox model)
sel_fit <- function(formula, data) {
  full <- do.call(coxph, list(formula, data = data))  # data embedded so step() can refit
  # keep treatment and its interactions in scope (always retained)
  step(full, scope = list(lower = ~ rx + rx:age + rx:bm + rx:ap, upper = formula), trace = 0)
}
sel_fit_pb <- function(formula, data) step(do.call(coxph, list(formula, data = data)), scope = list(lower = ~ trt, upper = formula), trace = 0)
boot632p_sel <- function(formula, data, selector, B = 200, seed = 1) {
  set.seed(seed)
  app <- cidx(selector(formula, data), data)
  oob <- mean(replicate(B, {
    i <- sample.int(nrow(data), replace = TRUE); o <- setdiff(seq_len(nrow(data)), i)
    tryCatch(cidx(selector(formula, data[i, ]), data[o, ]), error = function(e) NA)
  }), na.rm = TRUE)
  R <- if (app > oob) min(1, (app - oob) / (app - 0.5)) else 0
  w <- 0.632 / (1 - 0.368 * R)
  c(apparent = app, oob = oob, c632p = (1 - w) * app + w * oob)
}
res_sel <- rbind(prostate = boot632p_sel(f_pc, pc, sel_fit), pbc = boot632p_sel(f_pb, pb, sel_fit_pb))
print(round(res_sel, 3))
write.csv(round(res_sel, 3), "results/realdata_cox_cindex_aicselection.csv")
