# matchfamily: checklist for testing in jamovi

Menu: Matching (analysis group). Install: `jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`.
Open sample data from File > Open > Data Library > Matching and Weighting Family (`lalonde` is an .omv file with `age` and `educ` Continuous and `race` Nominal; if you import the CSV yourself, jamovi guesses `age` and `educ` as Nominal, which makes `nodegree` collinear: set them to Continuous).

## 1. Propensity Score and Overlap
Dataset: LaLonde. Treatment `treat`; Covariates `age`, `educ`, `race`, `married`, `nodegree`, `re74`, `re75`.
Expect: c statistic 0.874; median score 0.65 (treated) and 0.08 (control); the histogram shows controls piled near zero. Try trimming "Outside 0.1 to 0.9".
Try a perfect predictor (a copy of `treat`) to see the plain-language message.

## 2. Covariate Balance
Dataset: LaLonde, same variables. Adjustment: Inverse probability weights; Estimand: ATT.
Expect: before, largest absolute standardised difference 1.877 and 8 of 9 variables above 0.1; after, 0.119 (age) and 1 above; effective sample size of the controls about 100.
Then: Entropy balancing (everything about 0), Weights from a variable (select any positive variable), turn off group means, variance ratio, KS; the Love plot has two symbols and a dashed line.

## 3. Weighting and Doubly Robust Estimators
Dataset: Simulated confounded treatment. Outcome `y`; Treatment `treat`; Covariates `x1`, `x2`, `x3`; Estimand: All three (the default shows ATE, ATT and ATC rows for each estimator).
Expect (RA, IPW, AIPW, IPWRA): ATE 1.012, 1.037, 1.010, 1.010; ATT 1.173, 1.191, 1.172, 1.171; ATC 0.859, 0.886, 0.857, 0.858 (standard errors about 0.04 to 0.06).
Then: Probit, tick Potential-outcome means, add a cluster variable, a sampling weight; on LaLonde the estimators disagree widely and the weight diagnostics warn.

## 4. Matching Estimators
Dataset: Simulated confounded. Outcome `y`; Treatment `treat`; Matching covariates `x1`, `x2`, `x3`; one neighbour; Mahalanobis. The default Estimand shows ATE, ATT and ATC (tables carry an Estimand column); choose ATT for the numbers below.
Expect: 1.208 (s.e. 0.060); k = 3: 1.189 (0.052); with Bias adjustment and k = 1: 1.196 (0.059).
On LaLonde with ties: ATT k = 1 gives 262 (s.e. 1207), k = 3 gives 1112 (955); Exact-match variables `married`, `nodegree`; a caliper small enough to drop units shows the note.

## 5. Entropy Balancing
Dataset: Simulated confounded, Means; the default shows all three estimands, choose ATT for: 1.171 (s.e. 0.046) and a balance table with identical means. On LaLonde ATC with "Means and variances": the message that the target cannot be reached.

## 6. Compare Estimators
Dataset: Simulated confounded; the default shows all three estimands (the forest plot has one panel each), choose ATT for: unadjusted 1.704; RA 1.173; IPW 1.191; AIPW 1.172; IPWRA 1.171; matching 1.196; entropy balancing 1.171; the forest plot lists them. On LaLonde ATC with two moments, entropy balancing is reported in the notes and the others still show.

## 7. Sensitivity to Unobserved Confounding
Dataset: Simulated hidden confounder. Outcome `y`; Treatment `treat`; Covariates `x1`, `x2`; ATT.
Expect: AIPW effect 1.477; E-value 3.96 and 3.64 for the confidence limit; Rosenbaum: 908 pairs, critical Gamma 6.97 (signed rank) and 3.59 (sign test); the plot shows the upper bounds rising with a dashed line at 0.05.
Then: "Enter an estimate" with Risk ratio 1.5 (1.2, 1.9); an interval that includes 1 gives E-value 1 with the note.

## Overlap weights and CBPS (Weighting and Doubly Robust Estimators; Covariate Balance)
Dataset: Simulated confounded treatment. Outcome `y`; Treatment `treat`; Covariates `x1`, `x2`, `x3`. Estimand "Overlap-weighted effect (ATO)", all four estimators ticked.
Expect: three rows: RA 1.018 (s.e. 0.040), IPW 1.018 (0.040), AIPW 1.016 (0.040), and a note that IPWRA was skipped; with "Weight diagnostics" the weights lie between 0.04 and 0.94 and the effective sample sizes are about 1230 and 1290.
Then: Estimand "All three" with the model "Covariate balancing (CBPS, exact)": IPW 1.014 (0.042) ATE, 1.171 (0.046) ATT, 0.858 (0.047) ATC (the logit values were 1.037, 1.191, 0.886); the note explains CBPS. In Covariate Balance choose the same score with adjustment IPW: the weighted standardised differences are exactly 0 for the means.
Messages: ATO with only IPWRA ticked -> message to tick another estimator.

## Things to look at in the interface
- Collapsed groups open and close; options that do not apply (trimming bounds, manual input) are disabled or hidden.
- Messages for wrong input are in plain words (non-binary treatment, a variable in two roles, no estimator ticked).
- Tables fill as variables are dropped into the boxes and clear when a variable is removed.
- Number formats, p-values (< .001), long variable names, nominal covariates that expand into dummies.
- Plots (histogram, Love plot, forest, p-value bounds): legend, axis labels, theme and dark mode.
- The ribbon icon (after running `dev/icon_preview.ps1` as administrator).
