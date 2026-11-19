# install-deps.R -- CRAN packages used by the letter scripts.
p <- c("survival", "randomForestSRC", "ranger", "mcstatsim")
install.packages(setdiff(p, rownames(installed.packages())))
for (x in p) cat(x, if (requireNamespace(x, quietly = TRUE)) "ok" else "MISSING", "\n")
