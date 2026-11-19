# jamovi analysis class: entropy balancing.

ebalClass <- R6::R6Class(
    "ebalClass",
    inherit = ebalBase,
    private = list(

        .run = function() {
            opts <- self$options
            mf_show_text(self$results, opts, "ebal")
            if (is.null(opts$outcome) || is.null(opts$treat) || length(opts$covs) == 0) {
                self$results$notes$setContent("<p>Select an outcome, a binary treatment and at least one covariate.</p>")
                return()
            }
            pp <- mf_guard(mf_prepare(self$data, opts$treat, opts$covs, outcome = opts$outcome, weights = opts$weights, cluster = opts$cluster))
            Xc <- mf_design(pp$data, opts$covs, intercept = FALSE)
            mom <- as.numeric(substring(opts$moments, 2))
            crit <- stats::qnorm(1 - (1 - opts$level / 100) / 2)
            ests <- if (opts$estimand == "all") c("ate", "att", "atc") else opts$estimand
            elab <- c(ate = "ATE", att = "ATT", atc = "ATC")
            d <- pp$d
            bm <- mf_bal_matrix(pp$data, opts$covs)
            tb <- self$results$effects; pm <- self$results$pom; wt <- self$results$wdiag; bt <- self$results$bal
            notes <- character(0); nfail <- 0; kw <- 0; kb <- 0; nfeat <- NA
            for (i in seq_along(ests)) {
                e <- ests[i]
                r <- tryCatch(mf_ebal(pp$y, d, Xc, pp$sw, pp$cl, e, mom), error = function(err) conditionMessage(err))
                if (is.character(r)) {
                    if (length(ests) == 1) jmvcore::reject(r)
                    nfail <- nfail + 1
                    tb$addRow(rowKey = i, values = list(method = elab[[e]]))
                    notes <- c(notes, sprintf("%s could not be computed: %s", elab[[e]], r))
                    next
                }
                nfeat <- r$n_features
                z <- r$effect / r$se
                tb$addRow(rowKey = i, values = list(method = elab[[e]], estimate = r$effect, se = r$se, z = z, p = 2 * stats::pnorm(-abs(z)),
                                                    lo = r$effect - crit * r$se, hi = r$effect + crit * r$se))
                pm$addRow(rowKey = 2 * i - 1, values = list(estimand = elab[[e]], group = "Y(1)", estimate = r$pom1, se = r$se_pom1))
                pm$addRow(rowKey = 2 * i, values = list(estimand = elab[[e]], group = "Y(0)", estimate = r$pom0, se = r$se_pom0))
                if (isTRUE(opts$showWeights)) {
                    ws <- mf_weight_summary(r$weights, d)
                    for (j in seq_len(nrow(ws))) {
                        kw <- kw + 1
                        wt$addRow(rowKey = kw, values = list(estimand = elab[[e]], group = ws$group[j], n = ws$n[j], min = ws$min[j], mean = ws$mean[j], max = ws$max[j],
                                                             ess = ws$ess[j], share = ws$share_max[j]))
                    }
                }
                if (isTRUE(opts$showBalance)) {
                    b <- mf_balance(bm$M, bm$binary, d, pp$sw, r$weights, e, "auto", TRUE)
                    for (j in seq_len(nrow(b$raw))) {
                        kb <- kb + 1
                        bt$addRow(rowKey = kb, values = list(estimand = elab[[e]], variable = b$raw$variable[j], m1 = b$raw$mean1[j], m0 = b$raw$mean0[j], smd = b$raw$smd[j],
                                                             am1 = b$adj$mean1[j], am0 = b$adj$mean0[j], asmd = b$adj$smd[j], avr = b$adj$vr[j]))
                    }
                }
                ess <- c(mf_ess(r$weights[d == 1]), mf_ess(r$weights[d == 0])); n_g <- c(sum(d == 1), sum(d == 0))
                if (any(ess < 0.5 * n_g)) notes <- c(notes, sprintf("%s: the effective sample size of a group is below half its size, a few units carry much of the weight.", elab[[e]]))
                if (r$imbalance > 1e-6) notes <- c(notes, sprintf("%s: the balance conditions are not met to full precision (%s).", elab[[e]], format(r$imbalance, digits = 2)))
                if (r$moment_check > 1e-6 * max(1, stats::sd(pp$y))) notes <- c(notes, sprintf("%s: the estimating equations are not solved to full precision (%s).", elab[[e]], format(r$moment_check, digits = 2)))
            }
            if (nfail == length(ests)) jmvcore::reject(paste(notes, collapse = " "))
            tb$setNote("e", sprintf("Weighted difference in means after balancing %d moment(s) of the covariates (%s balancing conditions). Standard errors are M-estimation (sandwich) standard errors that include the estimation of the balancing weights%s. Normal-based intervals.",
                                    mom, if (is.na(nfeat)) "-" else as.character(nfeat), if (is.null(pp$cl)) "" else ", clustered"))
            if (isTRUE(opts$showWeights))
                wt$setNote("w", "Weights for each estimand; the group that is not reweighted keeps its sampling weights (1 without them). Weights are scaled to sum to the group size.")
            if (isTRUE(opts$showBalance))
                bt$setNote("b", "Standardised difference = difference in means divided by the standard deviation of the target group (pooled for the ATE) in the unadjusted sample.")
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
