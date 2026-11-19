# Checks of the E-value and Rosenbaum bounds against the EValue and rbounds packages. Run from the module root: Rscript tests/reference_m4.R
suppressMessages({library(EValue); library(rbounds)})
source("R/engine_sens.R")
fails <- 0
ok <- function(label, a, b, tol = 1e-8) { pass <- is.finite(a) && is.finite(b) && abs(a - b) < tol; if (!pass) fails <<- fails + 1; cat(sprintf("%-40s %.6f vs %.6f  %s\n", label, a, b, if (pass) "ok" else "FAIL")) }
suppressMessages({
    e <- mf_evalue(1.8, 1.2, 2.7, "rr"); r <- evalues.RR(1.8, 1.2, 2.7); ok("RR: point", e$e_point, r["E-values", "point"]); ok("RR: confidence limit", e$e_ci, r["E-values", "lower"])
    e <- mf_evalue(0.6, 0.4, 0.9, "rr"); r <- evalues.RR(0.6, 0.4, 0.9); ok("RR < 1: point", e$e_point, r["E-values", "point"]); ok("RR < 1: confidence limit", e$e_ci, r["E-values", "upper"])
    e <- mf_evalue(0.9, 0.7, 1.2, "rr"); ok("RR, interval includes the null", e$e_ci, 1)
    e <- mf_evalue(2.2, 1.4, 3.1, "or_common"); r <- evalues.OR(2.2, 1.4, 3.1, rare = FALSE); ok("OR common: point", e$e_point, r["E-values", "point"]); ok("OR common: limit", e$e_ci, r["E-values", "lower"])
    e <- mf_evalue(2.2, 1.4, 3.1, "or_rare"); r <- evalues.OR(2.2, 1.4, 3.1, rare = TRUE); ok("OR rare: point", e$e_point, r["E-values", "point"])
    e <- mf_evalue(1.7, 1.2, 2.4, "hr_common"); r <- evalues.HR(1.7, 1.2, 2.4, rare = FALSE); ok("HR common: point", e$e_point, r["E-values", "point"]); ok("HR common: limit", e$e_ci, r["E-values", "lower"])
    e <- mf_evalue(0.35, 0.12, 0.58, "smd"); r <- evalues.MD(0.35, se = (0.58 - 0.12) / (2 * 1.96)); ok("SMD: point", e$e_point, r["E-values", "point"])
    ok("SMD: limit (EValue rounds 1.96 * 0.91 to 1.78)", e$e_ci, r["E-values", "lower"], 2e-3)
})
set.seed(4); n <- 60; x <- rnorm(n, 1, 1.5); y <- rnorm(n, 0, 1.5)
rb <- psens(x, y, Gamma = 3, GammaInc = 0.5)$bounds
tb <- mf_rosenbaum(x - y, seq(1, 3, 0.5))$table
ok("Rosenbaum signed rank, upper bounds", max(abs(tb$wilcoxon_hi - rb[["Upper bound"]])), 0, 5e-5)
ok("Rosenbaum signed rank, lower bounds", max(abs(tb$wilcoxon_lo - rb[["Lower bound"]])), 0, 5e-5)
# the sign-test bound is the exact binomial tail P(S >= s); rbounds::binarysens uses P(S > s) and is therefore not used as a reference
s <- sum(x - y > 0); ok("sign test at Gamma = 1 equals the binomial test", tb$sign_hi[1], stats::binom.test(s, n, 0.5, alternative = "greater")$p.value)
cat(if (fails == 0) "\nAll sensitivity checks passed.\n" else sprintf("\n%d check(s) failed.\n", fails))
if (fails > 0) quit(status = 1)
