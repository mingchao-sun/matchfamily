# matchfamily 1.0.0

* First public release (jamovi library submission). Contains everything listed below; the earlier 0.x entries are development history.

# matchfamily 1.1.0

* **Overlap weights (ATO)** as an estimand in *Weighting and Doubly Robust Estimators* (RA, IPW and AIPW) and in *Covariate Balance* (Li, Morgan and Zaslavsky, 2018). The estimates and the M-estimation standard errors equal the R package `PSweight`.
* **Covariate balancing propensity score (CBPS, exactly identified)** as a propensity score model in the same two analyses (Imai and Ratkovic, 2014), separately for the ATE, the ATT and the ATC; the weighted covariate means balance exactly (the R package `CBPS` stops at about 1e-3), and the standard errors use the balancing conditions.

# matchfamily 1.0.0

First release. Seven analyses in the menu group Match:

* **Propensity Score and Overlap**: logit or probit score, c statistic, overlap by group, common support, trimming rules.
* **Covariate Balance**: standardised differences, variance ratios and KS statistics before and after IPW, entropy balancing or your own weights; Love plot.
* **Weighting and Doubly Robust Estimators**: RA, IPW, AIPW, IPWRA for the ATE, ATT and ATC with M-estimation standard errors, sampling weights, clusters and trimming.
* **Matching Estimators**: k nearest neighbours with replacement and ties, three distances, exact matching, caliper, bias adjustment, Abadie-Imbens standard errors.
* **Entropy Balancing**: exact balance of means, variances or skewness; effects with M-estimation standard errors.
* **Compare Estimators**: all of the above side by side, with a forest plot and the balance each method achieves.
* **Sensitivity to Unobserved Confounding**: E-value and Rosenbaum bounds.

Every analysis offers *Method notes* and *Data and setup guide* (Output section). Estimates agree with Stata `teffects` (RA, IPW, AIPW, IPWRA, `nnmatch`) to 1e-8,
with `ebalance` and R `WeightIt`, with `cobalt` for balance, with `EValue` and `rbounds` for the sensitivity analysis (see `README.md`).

Known limits: matching on the estimated propensity score is not offered (standard error not reproduced); IPWRA standard errors on the LaLonde data differ from
Stata's by 1 to 2 per cent, where Stata does not converge.
