# matchfamily — matching and weighting in jamovi

A jamovi module for the effect of a binary treatment when all confounders are observed: propensity scores and overlap, covariate balance, regression
adjustment, weighting, doubly robust estimators, matching, entropy balancing and sensitivity to unobserved confounding. The estimation engine is
written in base R (no dependency on `MatchIt`, `WeightIt`, `Matching` or `cobalt`) and is checked against Stata and R packages. Chinese summary: [README.zh.md](README.zh.md). Japanese summary: [README.ja.md](README.ja.md).

## Analyses (jamovi menu: Matching)

| Analysis | What it does |
|---|---|
| **Propensity Score and Overlap** | Logit, probit or covariate-balancing (CBPS) propensity score with robust standard errors, c statistic, score distribution by group, common support, trimming (common range, Crump 0.1 to 0.9, own bounds) |
| **Covariate Balance** | Standardised mean difference, variance ratio and Kolmogorov-Smirnov statistic for every covariate, before and after inverse probability weights, overlap weights, entropy balancing or weights from a variable; Love plot; count above a threshold; effective sample size |
| **Weighting and Doubly Robust Estimators** | Regression adjustment, IPW, AIPW and IPWRA for the ATE, ATT, ATC or the overlap-weighted ATO; M-estimation standard errors that include the estimated propensity score and regressions; sampling weights, clusters, trimming, weight diagnostics |
| **Matching Estimators** | k nearest neighbours with replacement and ties, Mahalanobis, inverse-variance or Euclidean distance, exact matching, caliper, bias adjustment, Abadie-Imbens standard errors; matches per unit and balance after matching |
| **Entropy Balancing** | Weights that match means (optionally variances and skewness) exactly; ATE, ATT, ATC; M-estimation standard errors |
| **Compare Estimators** | Unadjusted difference, RA, IPW, AIPW, IPWRA, matching and entropy balancing side by side with a forest plot and the balance each achieves |
| **Sensitivity to Unobserved Confounding** | E-value for an estimate from the data or typed in (standardised mean difference, risk, odds and hazard ratios); Rosenbaum bounds for matched pairs (sign and signed-rank tests, critical Gamma, plot) |

Every analysis has *Method notes* and *Data and setup guide* in its Output section.

## Install
- **jamovi library**: in jamovi open the modules button (top right), *jamovi library*, find *Matching and Weighting Family* and install (once the module is listed).
- **Sideload**: download the `.jmo` file and use *Install from file* in the same menu. A `.jmo` works only on the operating system, processor type and jamovi version it was built for.
- **From source** (needs R and the jamovi compiler `jmvtools`; adjust the path to your jamovi folder):
```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")
```

## Use from R (Rj editor, scripts)
Every analysis is also an R function, for example
```r
matchfamily::psmodel(data = d, treat = "treat", covs = c("x1", "x2", "x3"), link = "logit")
```

## Data format
One row per unit. A 0/1 treatment, a numeric outcome (not needed for the propensity score and balance analyses) and numeric or categorical covariates.
Optional sampling weights and cluster variable. Rows with missing values are dropped and reported.

## Example data
Shown in jamovi under **Open > Data Library > Matching and Weighting Family**; each entry says which analysis to use and which variable goes where.
Files in `data/` (`dev/make_data.R` rebuilds the simulated ones):

| File | Source | Used for |
|---|---|---|
| `lalonde.csv`, `lalonde.omv` | `MatchIt::lalonde` (Dehejia and Wahba sample: NSW treated, PSID controls) | Propensity score, balance, all estimators; the textbook case of poor overlap |
| `sim_confounded.csv` | simulated | Weighting, matching, entropy balancing, comparison (good overlap, heterogeneous effect) |
| `sim_overlap.csv` | simulated | Overlap and trimming |
| `sim_hidden.csv` | simulated, with a confounder left out of the data | Sensitivity analysis |

## References
Papers and R packages behind each analysis appear in the References section of the jamovi output (defined in `jamovi/00refs.yaml`).

## Validation
| Part | Reference | Agreement |
|---|---|---|
| RA, IPW, AIPW, IPWRA (ATE, ATT, ATC; logit, probit) | Stata `teffects ra / ipw / aipw / ipwra` | Estimates and standard errors to 1e-8, except IPWRA standard errors on LaLonde (1 to 2 per cent), where Stata reports `convergence not achieved` |
| Entropy balancing weights (1 to 3 moments) | `WeightIt`, Stata `ebalance` | Weights to 1e-6 (`ebalance` stops at its tolerance, 0.015 by default) |
| Balance statistics | `cobalt` | Exact for continuous variables; binary variables differ by the factor sqrt((n-1)/n) because cobalt uses sqrt(p(1-p)) |
| Nearest-neighbour matching (3 distances, k = 1 and 3, ATE, ATT, ATC, bias adjustment, exact matching, ties) | Stata `teffects nnmatch`; point estimates also `Matching` | Estimates and standard errors to 1e-8 (28 cases) |
| E-value | `EValue` | 1e-8 (the package rounds 1.96 x 0.91 to 1.78 for standardised differences) |
| Rosenbaum signed-rank bounds | `rbounds::psens` | Equal to the package's printed precision |
| Overlap weights (ATO) | `PSweight` | Estimates and standard errors agree |
| CBPS (exactly identified) | `CBPS` | Balance error 3e-15 (the package itself reaches about 6e-4) |

`../stata_validation/` holds the Stata scripts (`14_teffects.do`, `15_ebalance.do`, `16_nnmatch.do`) and the comparison scripts. `tests/smoke.R` runs every analysis headless (dozens of scenarios, including bad input).

## Known limitations
- All estimators assume that the covariates contain every confounder and that treated and control units overlap. The sensitivity analysis says how fragile that assumption is; it does not test it.
- The ATC is obtained by exchanging the groups.
- Matching is with replacement. Matching on the estimated propensity score (`teffects psmatch`) is **not** offered: the point estimates agree with Stata but the adjusted standard error of Abadie and Imbens (2016) could not be reproduced. Use *Propensity Score and Overlap* with *Weighting and Doubly Robust Estimators* instead.
- With a caliper or exact-match variables, units without a valid match are left out and counted (Stata stops with an error).
- The bootstrap is not used for matching estimators with a fixed number of neighbours (it is not valid for them).
- Not included: stratification, coarsened exact matching, optimal and full matching, matching without replacement, Cinelli-Hazlett robustness values and Oster's delta (these are in `sensfamily`).

## Notes for developers
- Layout: `R/engine_*.R` is plain R with no jamovi dependency; `R/*.b.R` and `jamovi/*.yaml` are the jamovi front end.
- Tests: `Rscript tests/smoke.R`, `tests/reference_m2.R` to `reference_m5.R` (need `cobalt`, `WeightIt`, `Matching`, `EValue`, `rbounds`, `PSweight`, `CBPS`, in a library named by `MF_REFLIB`); `tests/ui_checklist.md` lists what to check by hand in jamovi. Run the smoke test in jamovi's own R with the bundled libraries as well.
- Ribbon icon: `icon/analysis-match.svg` (74 x 44): treated units joined to the comparison units they are matched to. jmvtools 28.17 cannot attach an icon to a module, so `dev/icon_preview.ps1` adds one CSS rule for the menu group to jamovi's client stylesheet (administrator PowerShell, restart jamovi afterwards).

## How to cite
Sun, M. *matchfamily: Matching and Weighting Family for jamovi* (jamovi module). https://github.com/mingchao-sun/matchfamily. Please also cite the papers behind the analysis you use (listed in the output).

## License
GPL-3, as stated in `DESCRIPTION`; third-party code and data are listed in `NOTICE.md`. Bug reports and suggestions: https://github.com/mingchao-sun/matchfamily/issues
