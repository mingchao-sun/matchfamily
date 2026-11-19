# Overlap weights (ATO) against PSweight, and the covariate balancing propensity score (CBPS, exact) against the CBPS package; run from the module root:
#   Rscript tests/reference_m5.R   (needs PSweight and CBPS, e.g. in a library named by MF_REFLIB)
lib <- Sys.getenv("MF_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
for (f in c("engine_ps.R", "engine_weight.R")) source(file.path("R", f))
if (!all(vapply(c("PSweight", "CBPS"), requireNamespace, TRUE, quietly = TRUE))) { cat("PSweight or CBPS is not installed: skipped\n"); quit(status = 0) }
suppressMessages({library(PSweight); library(CBPS)})
nfail <- 0
cmp <- function(label, a, b, tol) { ok <- all(is.finite(a)) && max(abs(a - b)) < tol; cat(sprintf("[%s] %-62s max diff %.1e\n", if (ok) "PASS" else "FAIL", label, max(abs(a - b)))); if (!ok) nfail <<- nfail + 1 }
chk <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) nfail <<- nfail + 1 }

d <- read.csv("data/sim_confounded.csv"); pp <- mf_prepare(d, "treat", c("x1", "x2", "x3"), outcome = "y")

# ---- overlap weights (ATO) against PSweight
for (cfg in list(list("ipw", FALSE), list("aipw", TRUE))) {
    r <- mf_estimate(pp$y, pp$d, pp$X, pp$X, pp$sw, NULL, cfg[[1]], "ato", "logit")
    ps <- if (cfg[[2]]) PSweight(ps.formula = treat ~ x1 + x2 + x3, yname = "y", data = d, weight = "overlap", augmentation = TRUE, out.formula = y ~ x1 + x2 + x3)
          else PSweight(ps.formula = treat ~ x1 + x2 + x3, yname = "y", data = d, weight = "overlap")
    s <- summary(ps)$estimates
    cmp(sprintf("ATO %s: estimate and standard error equal PSweight", toupper(cfg[[1]])), c(r$effect, r$se), c(s[1, "Estimate"], s[1, "Std.Error"]), 5e-6)
}
# regression with overlap weights: the weighted mean of the outcome-model difference
r <- mf_estimate(pp$y, pp$d, pp$X, pp$X, pp$sw, NULL, "ra", "ato", "logit")
f1 <- lm(y ~ x1 + x2 + x3, data = d[d$treat == 1, ]); f0 <- lm(y ~ x1 + x2 + x3, data = d[d$treat == 0, ])
pr <- mf_ps_fit(pp$d, pp$X, pp$sw, "logit")$p; h <- pr * (1 - pr)
cmp("ATO regression adjustment equals the h-weighted outcome-model contrast", r$effect, sum(h * (predict(f1, d) - predict(f0, d))) / sum(h), 1e-8)
# exact balance of the weighted covariate means with overlap weights
w <- mf_weights(pp$d, pr, "ato"); X <- as.matrix(d[, c("x1", "x2", "x3")])
bal <- sapply(1:3, function(j) weighted.mean(X[d$treat == 1, j], w[d$treat == 1]) - weighted.mean(X[d$treat == 0, j], w[d$treat == 0]))
cmp("overlap weights balance the covariate means exactly (logit)", bal, rep(0, 3), 1e-9)
chk(inherits(try(mf_estimate(pp$y, pp$d, pp$X, pp$X, pp$sw, NULL, "ipwra", "ato", "logit"), silent = TRUE), "try-error"), "IPWRA with overlap weights is refused")

# ---- CBPS (exact) against the CBPS package: coefficients and weights
for (att in c(0, 1)) {
    est <- if (att == 1) "att" else "ate"
    mine <- mf_ps_fit(pp$d, pp$X, pp$sw, "cbps", est)
    ref <- suppressMessages(CBPS(treat ~ x1 + x2 + x3, data = d, ATT = att, method = "exact"))
    cmp(sprintf("CBPS %s: propensity scores agree with the package (its optimiser stops at 1e-3)", toupper(est)), mine$p, ref$fitted.values, 3e-3)
    wref <- mf_weights(pp$d, ref$fitted.values, est)
    imb <- function(w) if (att == 1) max(abs(sapply(1:3, function(j) mean(X[d$treat == 1, j]) - weighted.mean(X[d$treat == 0, j], w[d$treat == 0])))) else max(abs(sapply(1:3, function(j) weighted.mean(X[d$treat == 1, j], w[d$treat == 1]) - weighted.mean(X[d$treat == 0, j], w[d$treat == 0]))))
    chk(imb(mf_weights(pp$d, mine$p, est)) < imb(wref) / 100, sprintf("CBPS %s: the balance error here (%.0e) is far below the package's (%.0e)", toupper(est), imb(mf_weights(pp$d, mine$p, est)), imb(wref)))
    w <- mf_weights(pp$d, mine$p, est)
    bal <- sapply(1:3, function(j) weighted.mean(X[d$treat == 1, j], w[d$treat == 1]) - weighted.mean(X[d$treat == 0, j], w[d$treat == 0]))
    if (att == 1) bal <- sapply(1:3, function(j) mean(X[d$treat == 1, j]) - weighted.mean(X[d$treat == 0, j], w[d$treat == 0]))
    cmp(sprintf("CBPS %s: covariate means are balanced exactly", toupper(est)), bal, rep(0, 3), 1e-8)
}
# ATC: the swapped ATT; the weights make the treated mean equal the control mean
mc <- mf_ps_fit_for(pp$d, pp$X, pp$sw, "cbps", "atc"); wc <- mf_weights(pp$d, mc$p, "atc")
bal <- sapply(1:3, function(j) mean(X[d$treat == 0, j]) - weighted.mean(X[d$treat == 1, j], wc[d$treat == 1]))
cmp("CBPS ATC: the treated are balanced to the control means", bal, rep(0, 3), 1e-8)
# estimator: with the cbps link the IPW point estimate equals the weighted difference
for (est in c("ate", "att")) {
    r <- mf_estimate(pp$y, pp$d, pp$X, pp$X, pp$sw, NULL, "ipw", est, "cbps")
    p <- mf_ps_fit(pp$d, pp$X, pp$sw, "cbps", est)$p; w <- mf_weights(pp$d, p, est)
    tau <- weighted.mean(d$y[d$treat == 1], w[d$treat == 1]) - weighted.mean(d$y[d$treat == 0], w[d$treat == 0])
    cmp(sprintf("CBPS IPW %s: equals the difference of the weighted means", toupper(est)), r$effect, tau, 1e-8)
}
# standard errors: M-estimation (with the balancing moments as the propensity block) against the nonparametric bootstrap
set.seed(2); B <- 300
boot <- replicate(B, { i <- sample(nrow(d), replace = TRUE); q <- mf_prepare(d[i, ], "treat", c("x1", "x2", "x3"), outcome = "y")
                       tryCatch(mf_estimate(q$y, q$d, q$X, q$X, q$sw, NULL, "ipw", "ate", "cbps")$effect, error = function(e) NA) })
rr <- mf_estimate(pp$y, pp$d, pp$X, pp$X, pp$sw, NULL, "ipw", "ate", "cbps")
chk(abs(rr$se / sd(boot, na.rm = TRUE) - 1) < 0.15, sprintf("CBPS ATE standard error %.4f against bootstrap %.4f", rr$se, sd(boot, na.rm = TRUE)))
cat(if (nfail == 0) "\nAll overlap-weight and CBPS checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
