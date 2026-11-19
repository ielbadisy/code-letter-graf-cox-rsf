# code-letter-graf-cox-rsf

Code for our Letter to the Editor (Biometrical Journal) on Graf, Todd and Baksh (2026),
"Comparison of the Cox proportional hazards model and random survival forest algorithm for
predicting patient-specific survival probabilities in clinical trial data", 68:e70171.

The two reference trials are public: PBC (`survival::pbc`) and the VACURG advanced prostate
cancer trial (`data-ref/adv_prostate_ca/adv_prostate_ca.csv`, 475 patients).

## Scripts and the letter sections they support

| Script | Letter | Output in `results/` |
|---|---|---|
| `example-res.R` | 1. Cox C index on the real trials (apparent, out-of-bag, 0.632+, 10-fold CV, AIC selection inside the bootstrap) | `realdata_cox_cindex.csv`, `realdata_cox_cindex_aicselection.csv` |
| `dgp-check.R` | 2. Published prostate DGP: uncentred refit, treatment shift at the covariate means, simulated vs trial median | `dgp_check.csv` |
| `validate-dgp.R` | Check of `R/dgp-rct.R` on 20 000 patients before the simulation (KM vs true S, censoring rate) | console |
| `sim-run.R` | 3-4. Simulation: 48 scenarios x 200 replicates, 4 models | `sim_raw_200reps.csv.gz` |
| `sim-res.R` | 3-4. Means and MCSE, best model per scenario, counts quoted in the letter | `sim_summary.csv`, `sim_best_by_cell.csv`, `letter_table.csv`, `letter_numbers.csv` |

`R/dgp-rct.R` holds the data-generating models: Weibull event times on resampled trial
covariates, re-randomised 1:1 treatment, uniform censoring calibrated to 30% or 60%, and the
true survival function of every simulated patient. Covariates are centred at the trial means
(`centre = FALSE` gives the published uncentred version).

Models: Cox with main effects (as in Graf et al.), Cox stratified by arm with the true
treatment-covariate interactions, RSF with log-rank splitting (`randomForestSRC`) and RSF with
extremely randomised trees (`ranger`). Metrics: Harrell C on predicted risk at the median
event time, IPCW integrated Brier score, and mean absolute error against the true S(t | x)
over nine time points.

The letter (`letter.qmd`) reads every number it quotes from these `results/` files.

## Reproduce

Run everything from the repo root.

```bash
Rscript install-deps.R
Rscript example-res.R         # a few minutes
Rscript dgp-check.R
Rscript validate-dgp.R
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NCORES=14 Rscript sim-run.R 200
Rscript sim-res.R 200         # seconds; uses the committed raw results if sim-run.R is skipped
```

The simulation takes about 1 hour on 14 cores. Its raw output is committed, so `sim-res.R`
runs without it. Seeds are fixed (`runsim(seed = 2026)`).
