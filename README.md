# matchfamily: Matching and Weighting Family for jamovi

A jamovi module for treatment effects of a binary treatment when all confounders are observed, in the menu group **Match**. The estimation
engine is written in base R (no dependency on `MatchIt`, `WeightIt`, `Matching` or `cobalt`) and is checked against Stata (`teffects`,
`ebalance`) and R packages (`WeightIt`, `cobalt`, `Matching`, `EValue`, `rbounds`, `PSweight` for overlap weights, `CBPS`; `tests/reference_m5.R`).

| Analysis | What it does |
|---|---|
| **Propensity Score and Overlap** | Logit or probit propensity score with robust standard errors, c statistic, score distribution by group, common support, trimming (common range, Crump 0.1 to 0.9, own bounds) |
| **Covariate Balance** | Standardised mean difference, variance ratio and Kolmogorov-Smirnov statistic for every covariate, before and after inverse probability weights, entropy balancing or weights from a variable; Love plot; count above a threshold; effective sample size |
| **Weighting and Doubly Robust Estimators** | Regression adjustment, IPW, AIPW and IPWRA for the ATE, ATT or ATC; M-estimation standard errors that include the estimated propensity score and regressions; sampling weights, clusters, trimming, weight diagnostics |
| **Matching Estimators** | k nearest neighbours with replacement and ties, Mahalanobis, inverse-variance or Euclidean distance, exact matching, caliper, bias adjustment, Abadie-Imbens standard errors; matches per unit and balance after matching |
| **Entropy Balancing** | Weights that match means (optionally variances and skewness) exactly; ATE, ATT, ATC; M-estimation standard errors |
| **Compare Estimators** | Unadjusted difference, RA, IPW, AIPW, IPWRA, matching and entropy balancing side by side with a forest plot and the balance each achieves |
| **Sensitivity to Unobserved Confounding** | E-value for an estimate from the data or typed in (standardised mean difference, risk, odds and hazard ratios); Rosenbaum bounds for matched pairs (sign and signed-rank tests, critical Gamma, plot) |

Every analysis has *Method notes* and *Data and setup guide* in its Output section. Four example data sets are in the Data Library.

## Example data (`data/`)

| File | Source | Used for |
|---|---|---|
| `lalonde.csv`, `lalonde.omv` | `MatchIt::lalonde` (Dehejia and Wahba sample: NSW treated, PSID controls) | Propensity score, balance, all estimators; the textbook case of poor overlap |
| `sim_confounded.csv` | simulated | Weighting, matching, entropy balancing, comparison (good overlap, heterogeneous effect) |
| `sim_overlap.csv` | simulated | Overlap and trimming |
| `sim_hidden.csv` | simulated, with a confounder left out of the data | Sensitivity analysis |

`dev/make_data.R` rebuilds them. `data/` documents the expected results in each analysis's *Data and setup guide*.

## Validation

| Part | Reference | Agreement |
|---|---|---|
| RA, IPW, AIPW, IPWRA (ATE, ATT, ATC; logit, probit) | Stata `teffects ra / ipw / aipw / ipwra` | Estimates and standard errors to 1e-8, except IPWRA standard errors on LaLonde (1 to 2 per cent), where Stata reports `convergence not achieved` |
| Entropy balancing weights (1 to 3 moments) | `WeightIt`, Stata `ebalance` | Weights to 1e-6 (`ebalance` stops at its tolerance, 0.015 by default) |
| Balance statistics | `cobalt` | Exact for continuous variables; binary variables differ by the factor sqrt((n-1)/n) because cobalt uses sqrt(p(1-p)) |
| Nearest-neighbour matching (3 distances, k = 1 and 3, ATE, ATT, ATC, bias adjustment, exact matching, ties) | Stata `teffects nnmatch`; point estimates also `Matching` | Estimates and standard errors to 1e-8 (28 cases) |
| E-value | `EValue` | 1e-8 (the package rounds 1.96 x 0.91 to 1.78 for standardised differences) |
| Rosenbaum signed-rank bounds | `rbounds::psens` | Equal to the package's printed precision |

`../stata_validation/` holds the Stata scripts (`14_teffects.do`, `15_ebalance.do`, `16_nnmatch.do`) and the comparison scripts
(`compare_teffects.R`, `compare_ebal.R`, `compare_nnmatch.R`).

## Design notes and limits

- All estimators assume that the covariates contain every confounder and that treated and control units overlap. The sensitivity analysis
  says how fragile that assumption is; it does not test it.
- The ATC is obtained by exchanging the groups.
- Matching is with replacement. Matching on the estimated propensity score (`teffects psmatch`) is **not** offered: the point estimates agree
  with Stata but the adjusted standard error of Abadie and Imbens (2016) could not be reproduced (see `PLAN.md`; `mf_psmatch` in
  `R/engine_match.R` is experimental). Use *Propensity Score and Overlap* with *Weighting and Doubly Robust Estimators* instead.
- With a caliper or exact-match variables, units without a valid match are left out and counted (Stata stops with an error).
- The bootstrap is not used for matching estimators with a fixed number of neighbours (it is not valid for them).
- Not included: stratification, overlap weights (ATO), CBPS, coarsened exact matching, optimal and full matching, matching without
  replacement, Cinelli-Hazlett robustness values and Oster's delta.

## Build and test

```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")
```

`Rscript tests/smoke.R` runs every analysis headless (about 85 scenarios, including bad input). `tests/reference_m2.R` (balance and entropy
balancing against `cobalt` and `WeightIt`), `tests/reference_m3.R` (matching against `Matching`) and `tests/reference_m4.R` (E-value and
Rosenbaum bounds against `EValue` and `rbounds`) need those packages. Run the smoke test in jamovi's own R with the bundled libraries too
(see `PLAN.md`). `tests/ui_checklist.md` lists what to look at by hand in jamovi.

## Ribbon icon

`icon/analysis-match.svg` (74 x 44, same size and colours as the other families): treated units (orange) joined to the comparison units
(blue) they are matched to, with circles of different size for the weight a unit carries. jmvtools 28.17 cannot attach an icon to a module, so
`dev/icon_preview.ps1` adds one CSS rule for the menu group Matching to jamovi's client stylesheet (run it in an administrator PowerShell, see the
comments at the top; restart jamovi afterwards).
