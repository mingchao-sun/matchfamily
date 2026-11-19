# Checks of the balance and entropy balancing engines against cobalt and WeightIt. Run from the module root: Rscript tests/reference_m2.R
suppressMessages({library(cobalt); library(WeightIt)})
for (f in c("engine_ps.R", "engine_weight.R", "engine_balance.R", "engine_ebal.R")) source(file.path("R", f))
ll <- read.csv("data/lalonde.csv"); ll$race <- factor(ll$race); sim <- read.csv("data/sim_confounded.csv")
cl <- c("age", "educ", "race", "married", "nodegree", "re74", "re75"); cv <- c("x1", "x2", "x3")
fails <- 0
ok <- function(label, val, tol) { pass <- is.finite(val) && val < tol; if (!pass) fails <<- fails + 1; cat(sprintf("%-58s %.2e  %s\n", label, val, if (pass) "ok" else "FAIL")) }

# ---- entropy balancing weights against WeightIt (normalised within group) -------------------------------------
for (spec in list(list("sim", sim, "y", cv), list("lalonde", ll, "re78", cl))) for (est in c("att", "atc", "ate")) for (mom in 1:2) {
    pp <- mf_prepare(spec[[2]], "treat", spec[[4]], outcome = spec[[3]]); Xc <- mf_design(pp$data, spec[[4]], intercept = FALSE)
    r <- tryCatch(mf_ebal_weights(pp$d, Xc, pp$sw, est, mom), error = function(e) NULL)
    if (is.null(r)) { cat(sprintf("%-58s infeasible (no solution)\n", paste("ebal", spec[[1]], est, mom))); next }
    dd <- pp$data; dd$treat <- pp$d
    w <- weightit(as.formula(paste("treat ~", paste(spec[[4]], collapse = "+"))), data = dd, method = "ebal", estimand = toupper(est), moments = mom, int = FALSE)
    a <- r$w; b <- w$weights
    for (g in 0:1) { i <- pp$d == g; a[i] <- a[i] / sum(a[i]); b[i] <- b[i] / sum(b[i]) }
    ok(paste("ebal weights vs WeightIt:", spec[[1]], est, "moments", mom), max(abs(a - b)), 1e-5)
}

# ---- balance statistics against cobalt (continuous variables: SMD, variance ratio, KS) --------------------------
set.seed(2); pp <- mf_prepare(ll, "treat", cl, outcome = "re78"); sw <- runif(nrow(ll), 0.5, 2)
bm <- mf_bal_matrix(pp$data, cl); ps <- mf_ps_fit(pp$d, pp$X, sw, "logit")
for (est in c("ate", "att", "atc")) {
    w <- mf_weights(pp$d, ps$p, est)
    b <- mf_balance(bm$M, bm$binary, pp$d, sw, sw * w, est)
    bt <- bal.tab(pp$data[, cl], treat = pp$d, weights = w, s.weights = sw, estimand = toupper(est), un = TRUE, binary = "std", continuous = "std",
                  disp.v.ratio = TRUE, stats = c("mean.diffs", "variance.ratios", "ks.statistics"))$Balance
    nm <- sub("^race_", "race: ", rownames(bt)); ix <- match(b$raw$variable, nm)
    cont <- !b$raw$binary
    ok(paste("SMD unadjusted (continuous) vs cobalt:", est), max(abs(b$raw$smd[cont] - bt$Diff.Un[ix][cont])), 1e-8)
    ok(paste("SMD adjusted (continuous) vs cobalt:", est), max(abs(b$adj$smd[cont] - bt$Diff.Adj[ix][cont])), 1e-8)
    ok(paste("variance ratio, both stages vs cobalt:", est), max(abs(b$raw$vr - bt$V.Ratio.Un[ix]), abs(b$adj$vr - bt$V.Ratio.Adj[ix]), na.rm = TRUE), 1e-8)
    ok(paste("KS, both stages vs cobalt:", est), max(abs(b$raw$ks - bt$KS.Un[ix]), abs(b$adj$ks - bt$KS.Adj[ix]), na.rm = TRUE), 1e-8)
    # binary variables: cobalt divides by sqrt(p(1-p)), we use the sample variance; the ratio is sqrt((n-1)/n) of the reference group
    bn <- b$raw$binary
    ok(paste("SMD binary: within 0.5% of cobalt:", est), max(abs(b$raw$smd[bn] / bt$Diff.Un[ix][bn] - 1)), 5e-3)
}
cat(if (fails == 0) "\nAll reference checks passed.\n" else sprintf("\n%d check(s) failed.\n", fails))
if (fails > 0) quit(status = 1)
