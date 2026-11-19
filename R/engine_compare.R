# Side-by-side comparison of the estimators of the matching and weighting family on the same data and estimand. Base R only.

mf_compare_methods <- function() c(naive = "Unadjusted difference in means", ra = "Regression adjustment (RA)", ipw = "Inverse probability weighting (IPW)",
                                   aipw = "Augmented IPW (AIPW)", ipwra = "IPW regression adjustment (IPWRA)",
                                   nn = "Nearest-neighbour matching", ebal = "Entropy balancing")

# pp: prepared sample (mf_prepare with an outcome), covs: covariate names, methods: subset of names(mf_compare_methods()).
# Each estimator that fails is reported with its message instead of stopping the comparison.
mf_compare <- function(pp, covs, methods, estimand, link = "logit", k = 1, biasadj = TRUE, moments = 1, maha = "maha") {
    y <- pp$y; d <- pp$d; sw <- pp$sw
    Xc <- mf_design(pp$data, covs, intercept = FALSE)
    bm <- mf_bal_matrix(pp$data, covs)
    out <- list()
    one <- function(m) {
        switch(m,
               naive = { r <- mf_naive(y, d, sw); list(effect = r$effect, se = r$se, w = NULL) },
               ra = , aipw = , ipwra = { r <- mf_estimate(y, d, pp$X, pp$X, sw, NULL, m, estimand, link); list(effect = r$effect, se = r$se, w = NULL) },
               ipw = {
                   r <- mf_estimate(y, d, pp$X, pp$X, sw, NULL, "ipw", estimand, link)
                   ps <- mf_ps_fit(d, pp$X, sw, link)
                   list(effect = r$effect, se = r$se, w = sw * mf_weights(d, ps$p, estimand))
               },
               nn = {
                   r <- mf_nnmatch(y, d, Xc, estimand, k, maha, NULL, NULL, if (biasadj) pp$X else NULL, 2)
                   list(effect = r$effect, se = r$se, w = r$weights, dropped = r$n_dropped)
               },
               ebal = { r <- mf_ebal(y, d, Xc, sw, NULL, estimand, moments); list(effect = r$effect, se = r$se, w = r$weights) })
    }
    for (m in methods) {
        res <- tryCatch(one(m), error = function(e) list(error = conditionMessage(e)))
        if (!is.null(res$error)) { out[[m]] <- list(method = m, error = res$error); next }
        smd_max <- NA_real_; above <- NA_integer_
        w <- if (m == "naive") sw else res$w
        if (!is.null(w)) {
            b <- mf_balance(bm$M, bm$binary, d, sw, w, estimand, "auto", TRUE)
            tb <- if (m == "naive") b$raw else b$adj
            smd_max <- max(abs(tb$smd), na.rm = TRUE); above <- sum(abs(tb$smd) > 0.1, na.rm = TRUE)
        }
        out[[m]] <- list(method = m, effect = res$effect, se = res$se, smd_max = smd_max, above = above, dropped = res$dropped)
    }
    out
}
