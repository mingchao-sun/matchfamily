# Third-party software and data

matchfamily is distributed under the GNU General Public License, version 3 (see `LICENSE`).

## Code

The estimators and tests are implemented from the published formulas in base R. No code of other packages is copied. Where a statistic
is defined by an algorithm rather than a closed form, the implementation follows the published description and was checked against the
reference implementation:

| Method | Reference | Checked against |
|---|---|---|
| RA, IPW, AIPW, IPWRA and their M-estimation standard errors | Wooldridge (2007), Robins, Rotnitzky and Zhao (1994), Imbens and Rubin (2015) | Stata `teffects` |
| Entropy balancing | Hainmueller (2012) | `WeightIt` (R), `ebalance` (Stata) |
| Nearest-neighbour matching, bias adjustment, standard errors | Abadie and Imbens (2006, 2011); the conditional variance and tie handling follow the description in the Stata manual for `teffects nnmatch` | Stata `teffects nnmatch`, `Matching` (R) |
| Balance statistics | Austin (2009) | `cobalt` |
| E-value | VanderWeele and Ding (2017) | `EValue` |
| Rosenbaum bounds | Rosenbaum (2002) | `rbounds` |

jmvcore, R6 and ggplot2 are used for the jamovi interface and the plots. The reference packages and Stata are used only for testing.

## Example data

`lalonde.csv` is derived from the data set `lalonde` of the CRAN package `MatchIt` (Ho, Imai, King and Stuart, GPL-2 | GPL-3), reproduced here under that
licence. Original source: LaLonde (1986) and Dehejia and Wahba (1999) for the sample. `sim_confounded.csv`, `sim_overlap.csv` and `sim_hidden.csv` are simulated for this module.
