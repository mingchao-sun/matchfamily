# Checks of the nearest-neighbour matching engine against the Matching package. Run from the module root: Rscript tests/reference_m3.R
# Point estimates must agree where there are no near-ties. Standard errors (and bias-adjusted estimates) differ from Matching, which follows
# other conventions; the module agrees with Stata's teffects nnmatch to 1e-8 on all of them (stata_validation/16_nnmatch.do, compare_nnmatch.R).
suppressMessages(library(Matching))
for (f in c("engine_ps.R", "engine_match.R")) source(file.path("R", f))
set.seed(5); n <- 400
x1 <- rnorm(n); x2 <- sample(1:4, n, TRUE); x3 <- rnorm(n)
d <- rbinom(n, 1, plogis(0.6 * x1 + 0.3 * (x2 - 2.5))); y <- 1 + d + x1 + 0.5 * x2 + x3 + rnorm(n); X <- cbind(x1, x2, x3)
fails <- 0
for (est in c("ATE", "ATT", "ATC")) for (M in 1:4) {
    m <- Match(Y = y, Tr = d, X = X, estimand = est, M = M, Weight = 2, ties = TRUE, Var.calc = 2, distance.tolerance = 1e-12)
    r <- mf_nnmatch(y, d, X, tolower(est), k = M, "maha", J = 2)
    diff <- abs(r$effect - m$est); pass <- diff < 1e-8; if (!pass) fails <- fails + 1
    cat(sprintf("%s M=%d  effect %.6f vs %.6f  diff %.1e  %s   (se %.4f vs %.4f)\n", est, M, r$effect, m$est, diff, if (pass) "ok" else "FAIL", r$se, m$se))
}
cat(if (fails == 0) "\nPoint estimates agree with Matching.\n" else sprintf("\n%d check(s) failed.\n", fails))
if (fails > 0) quit(status = 1)
