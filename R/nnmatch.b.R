# jamovi analysis class: nearest-neighbour matching estimators.

nnmatchClass <- R6::R6Class(
    "nnmatchClass",
    inherit = nnmatchBase,
    private = list(

        .run = function() {
            opts <- self$options
            mf_show_text(self$results, opts, "nnmatch")
            if (is.null(opts$outcome) || is.null(opts$treat) || length(opts$covs) == 0) {
                self$results$notes$setContent("<p>Select an outcome, a binary treatment and at least one matching covariate.</p>")
                return()
            }
            adj <- if (isTRUE(opts$biasadj)) (if (length(opts$adjvars)) opts$adjvars else opts$covs) else NULL
            pp <- mf_guard(mf_prepare(self$data, opts$treat, opts$covs, outcome = opts$outcome, extra = opts$adjvars, exact = opts$exact))
            Xm <- mf_design(pp$data, opts$covs, intercept = FALSE)
            Zb <- if (!is.null(adj)) mf_guard(mf_design(pp$data, adj, intercept = TRUE)) else NULL
            key <- if (length(opts$exact)) do.call(paste, c(lapply(opts$exact, function(v) as.character(pp$data[[v]])), sep = "|")) else NULL
            cal <- if (opts$caliper > 0) opts$caliper else NULL
            crit <- stats::qnorm(1 - (1 - opts$level / 100) / 2)
            ests <- if (opts$estimand == "all") c("ate", "att", "atc") else opts$estimand
            elab <- c(ate = "ATE", att = "ATT", atc = "ATC")
            d <- pp$d
            bm <- mf_bal_matrix(pp$data, opts$covs)
            tb <- self$results$effects; sm <- self$results$summary; ut <- self$results$uses; bt <- self$results$bal
            notes <- character(0); nfail <- 0; ks <- 0; ku <- 0; kb <- 0
            for (i in seq_along(ests)) {
                e <- ests[i]
                r <- tryCatch(mf_nnmatch(pp$y, d, Xm, e, opts$k, opts$metric, cal, key, Zb, opts$varNn), error = function(err) conditionMessage(err))
                if (is.character(r)) {
                    if (length(ests) == 1) jmvcore::reject(r)
                    nfail <- nfail + 1
                    tb$addRow(rowKey = i, values = list(method = elab[[e]]))
                    notes <- c(notes, sprintf("%s could not be computed: %s", elab[[e]], r))
                    next
                }
                z <- r$effect / r$se
                tb$addRow(rowKey = i, values = list(method = elab[[e]], estimate = r$effect, se = r$se, z = z, p = 2 * stats::pnorm(-abs(z)),
                                                    lo = r$effect - crit * r$se, hi = r$effect + crit * r$se))
                nm <- r$n_matches
                items <- list(c("Observations", format(length(d), big.mark = ",")),
                              c("Treated / control", sprintf("%s / %s", sum(d == 1), sum(d == 0))),
                              c("Units with a counterfactual from matching", format(r$n_used, big.mark = ",")),
                              c("Units without a valid match (left out)", format(r$n_dropped, big.mark = ",")),
                              c("Matches per unit (mean / maximum, ties included)", sprintf("%s / %d", mf_fmt(mean(nm), 2), max(nm))),
                              c("Units with more matches than requested because of ties", format(sum(nm > opts$k), big.mark = ",")))
                for (j in seq_along(items)) { ks <- ks + 1; sm$addRow(rowKey = ks, values = list(estimand = elab[[e]], item = items[[j]][1], value = items[[j]][2])) }
                if (isTRUE(opts$showMatches)) {
                    for (g in c(1, 0)) {
                        kk <- r$K[d == g]; ku <- ku + 1
                        ut$addRow(rowKey = ku, values = list(estimand = elab[[e]], group = if (g == 1) "Treated" else "Control", n0 = sum(kk < 1e-9), n1 = sum(kk > 1e-9 & kk < 1.5),
                                                             n2 = sum(kk >= 1.5 & kk < 2.5), n3 = sum(kk >= 2.5 & kk < 4.5), n5 = sum(kk >= 4.5), max = max(kk)))
                    }
                }
                if (isTRUE(opts$showBalance)) {
                    b <- mf_balance(bm$M, bm$binary, d, rep(1, length(d)), r$weights, e, "auto", TRUE)
                    for (j in seq_len(nrow(b$raw))) {
                        kb <- kb + 1
                        bt$addRow(rowKey = kb, values = list(estimand = elab[[e]], variable = b$raw$variable[j], m1 = b$raw$mean1[j], m0 = b$raw$mean0[j], smd = b$raw$smd[j],
                                                             am1 = b$adj$mean1[j], am0 = b$adj$mean0[j], asmd = b$adj$smd[j], avr = b$adj$vr[j]))
                    }
                }
                if (r$n_dropped > 0) notes <- c(notes, sprintf("%s: %d unit(s) have no valid match (caliper or exact-match conditions) and were left out of the estimate.", elab[[e]], r$n_dropped))
                if (max(r$K) > 0.2 * length(d)) notes <- c(notes, sprintf("%s: one unit serves as the match for more than 20%% of the sample; the estimate leans heavily on it.", elab[[e]]))
            }
            if (nfail == length(ests)) jmvcore::reject(paste(notes, collapse = " "))
            tb$setNote("e", sprintf("Matching with replacement on %s distance, %d nearest neighbour(s) with ties included%s. Abadie-Imbens heteroskedasticity-consistent standard error using %d neighbours for the conditional variances. Normal-based intervals. The ATC is obtained by exchanging the groups.",
                                    c(maha = "Mahalanobis", ivar = "inverse-variance", eucl = "Euclidean")[[opts$metric]], opts$k,
                                    if (isTRUE(opts$biasadj)) ", with bias adjustment" else "", opts$varNn))
            if (isTRUE(opts$showMatches))
                ut$setNote("u", "Number of times a unit is used as a match (tied matches share the use). Heavy reuse of a few units makes the estimate depend on them.")
            if (isTRUE(opts$showBalance))
                bt$setNote("b", "Matched columns use the matching weights: a unit counts once if it has a counterfactual, plus the number of times it serves as a match.")
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (opts$metric == "maha" && ncol(Xm) > 10)
                notes <- c(notes, "With many covariates the nearest neighbours are far away and matching is biased; consider the bias adjustment.")
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
