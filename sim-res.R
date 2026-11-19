# sim-res.R -- summarise sim-run.R output.
# Outputs: full summary (mean, MCSE, failure rate) per cell x model x metric,
# a compact letter table, and the counts quoted in the letter (results/letter_numbers.csv):
# how often the C index / IBS and the truth-based error disagree on the best model, etc.
# Run from the repo root: Rscript sim-res.R [reps]


args <- commandArgs(TRUE)
reps <- if (length(args)) as.integer(args[1]) else 200L
raw <- read.csv(sprintf("results/sim_raw_%dreps.csv.gz", reps), stringsAsFactors = FALSE)
keys <- c("dataset", "hazard", "n", "cens", "model", "metric")

cell_id <- do.call(paste, raw[keys])
summ <- do.call(rbind, lapply(split(raw, cell_id), function(d) {
  v <- d$value[is.finite(d$value)]
  data.frame(d[1, keys], mean = mean(v), mcse = sd(v) / sqrt(length(v)),
             fail = 1 - length(v) / nrow(d), reps = nrow(d))
}))
rownames(summ) <- NULL
summ <- summ[order(summ$dataset, summ$hazard, summ$n, summ$cens, summ$metric, summ$model), ]
write.csv(summ, "results/sim_summary.csv", row.names = FALSE)

# best model per cell under each metric (lower is better except the C indices)
higher_better <- c(harrell_c = TRUE, uno_c = TRUE, ibs = FALSE, mae_true = FALSE, mae_true_med = FALSE)
cells <- split(summ, do.call(paste, summ[c("dataset", "hazard", "n", "cens")]))
best <- do.call(rbind, lapply(cells, function(d) {
  b <- sapply(names(higher_better), function(m) {
    x <- d[d$metric == m, ]
    x$model[if (higher_better[[m]]) which.max(x$mean) else which.min(x$mean)]
  })
  data.frame(d[1, c("dataset", "hazard", "n", "cens")], t(b))
}))
rownames(best) <- NULL
write.csv(best, "results/sim_best_by_cell.csv", row.names = FALSE)
# counts quoted in the letter (Sections 3-4)
mt <- summ[summ$metric == "mae_true" & summ$hazard == "nph", ]
w <- reshape(mt[c("dataset", "n", "cens", "model", "mean")], idvar = c("dataset", "n", "cens"),
             timevar = "model", direction = "wide")
rsf <- grepl("^rsf", best$mae_true)
num <- c(n_cells = nrow(best), c_vs_mae_differ = sum(best$harrell_c != best$mae_true),
         ibs_vs_mae_differ = sum(best$ibs != best$mae_true), forest_best = sum(rsf),
         forest_best_max_n = max(best$n[rsf]), nph_cells = nrow(w),
         nph_spec_beats_forests = sum(w$mean.cox_spec < pmin(w$mean.rsf_logrank, w$mean.rsf_extratrees)))
write.csv(data.frame(quantity = names(num), value = num), "results/letter_numbers.csv", row.names = FALSE)
print(num)

# letter table: N = 200, 30% censoring, all hazards, both trials
lt <- summ[summ$n == 200 & summ$cens == 0.3 & summ$metric %in% c("harrell_c", "ibs", "mae_true"), ]
lt$cell <- sprintf("%.3f (%.3f)", lt$mean, lt$mcse)
wide <- reshape(lt[c("dataset", "hazard", "model", "metric", "cell")],
                idvar = c("dataset", "hazard", "model"), timevar = "metric", direction = "wide")
names(wide) <- sub("^cell\\.", "", names(wide))
wide$hazard <- factor(wide$hazard, c("dec", "const", "inc", "nph"))
wide <- wide[order(wide$dataset, wide$hazard, wide$model), ]
write.csv(wide, "results/letter_table.csv", row.names = FALSE)
print(wide, row.names = FALSE)
