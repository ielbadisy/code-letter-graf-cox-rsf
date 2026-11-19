# sim-res.R -- summarise sim-run.R output. Input: results/sim_raw_<reps>reps.csv.gz.
# Outputs in results/: sim_summary.csv (mean, MCSE per scenario x model x metric),
# sim_best_by_cell.csv (best model per scenario and metric), letter_table.csv
# (N = 200, 30% censoring), letter_numbers.csv (counts quoted in Sections 3-4).
# Run from the repo root: Rscript sim-res.R 200
reps <- as.integer(commandArgs(TRUE)[1])
raw <- read.csv(sprintf("results/sim_raw_%dreps.csv.gz", reps))
CELL <- c("dataset", "hazard", "n", "cens")
HIGHER_BETTER <- c(harrell_c = TRUE, ibs = FALSE, mae_true = FALSE)

summ <- aggregate(value ~ ., raw[c(CELL, "model", "metric", "value")],
                  function(v) c(mean = mean(v), mcse = sd(v) / sqrt(length(v))))
summ <- do.call(data.frame, summ); names(summ) <- sub("^value\\.", "", names(summ))
summ <- summ[do.call(order, summ[c(CELL, "metric", "model")]), ]
write.csv(summ, "results/sim_summary.csv", row.names = FALSE)

best <- do.call(rbind, lapply(split(summ, summ[CELL], drop = TRUE), function(d)
  data.frame(d[1, CELL], as.list(vapply(names(HIGHER_BETTER), function(m) {
    x <- d[d$metric == m, ]; x$model[which.max(if (HIGHER_BETTER[[m]]) x$mean else -x$mean)]
  }, ""))))); rownames(best) <- NULL
write.csv(best, "results/sim_best_by_cell.csv", row.names = FALSE)

nph <- with(subset(summ, metric == "mae_true" & hazard == "nph"), tapply(mean, list(paste(dataset, n, cens), model), c))
rsf <- grepl("^rsf", best$mae_true)
num <- c(n_cells = nrow(best), c_vs_mae_differ = sum(best$harrell_c != best$mae_true),
         ibs_vs_mae_differ = sum(best$ibs != best$mae_true), forest_best = sum(rsf),
         forest_best_max_n = max(best$n[rsf]), nph_cells = nrow(nph),
         nph_spec_beats_forests = sum(nph[, "cox_spec"] < pmin(nph[, "rsf_logrank"], nph[, "rsf_extratrees"])))
write.csv(data.frame(quantity = names(num), value = num), "results/letter_numbers.csv", row.names = FALSE)

lt <- subset(summ, n == 200 & cens == 0.3)
lt$cell <- sprintf("%.3f (%.3f)", lt$mean, lt$mcse)
wide <- reshape(lt[c("dataset", "hazard", "model", "metric", "cell")], direction = "wide",
                idvar = c("dataset", "hazard", "model"), timevar = "metric")
names(wide) <- sub("^cell\\.", "", names(wide))
wide <- wide[order(wide$dataset, match(wide$hazard, c("dec", "const", "inc", "nph")), wide$model), ]
write.csv(wide, "results/letter_table.csv", row.names = FALSE)
