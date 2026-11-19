# jamovi analysis class: regression adjustment, IPW, AIPW and IPWRA.

wtestClass <- R6::R6Class(
    "wtestClass",
    inherit = wtestBase,
    private = list(

        .run = function() {
            opts <- self$options
            mf_show_text(self$results, opts, "wtest")
            if (is.null(opts$outcome) || is.null(opts$treat) || length(opts$covs) == 0) {
                self$results$notes$setContent("<p>Select an outcome, a binary treatment and at least one covariate.</p>")
                return()
            }
            methods <- c(ra = isTRUE(opts$doRa), ipw = isTRUE(opts$doIpw), aipw = isTRUE(opts$doAipw), ipwra = isTRUE(opts$doIpwra))
            if (!any(methods)) jmvcore::reject("Tick at least one estimator.")
            if (identical(opts$estimand, "ato") && !any(methods[c("ra", "ipw", "aipw")])) jmvcore::reject("Overlap weights are not available for IPWRA: tick regression adjustment, IPW or AIPW.")
            pp <- mf_guard(mf_prepare(self$data, opts$treat, opts$covs, outcome = opts$outcome, extra = opts$extra,
                                      weights = opts$weights, cluster = opts$cluster))
            tr <- mf_guard(mf_trimmed(pp, opts$link, opts$trim, opts$trimLo, opts$trimHi))
            q <- tr$pp
            crit <- stats::qnorm(1 - (1 - opts$level / 100) / 2)
            labs <- c(ra = "Regression adjustment (RA)", ipw = "Inverse probability weighting (IPW)", aipw = "Augmented IPW (AIPW)",
                      ipwra = "IPW regression adjustment (IPWRA)")
            ests <- if (opts$estimand == "all") c("ate", "att", "atc") else opts$estimand
            elab <- c(ate = "ATE", att = "ATT", atc = "ATC", ato = "ATO (overlap)")
            tb <- self$results$effects
            pm <- self$results$pom
            notes <- character(0)
            k <- 0; kw <- 0
            for (e in ests) for (m in names(methods)[methods]) {
                if (e == "ato" && m == "ipwra") { notes <- c(notes, "IPWRA is not available for overlap weights and was skipped."); next }
                private$.checkpoint()
                r <- mf_guard(mf_estimate(q$y, q$d, q$X, q$W, q$sw, q$cl, m, e, opts$link))
                k <- k + 1
                z <- r$effect / r$se
                tb$addRow(rowKey = k, values = list(estimand = elab[[e]], method = labs[[m]], estimate = r$effect, se = r$se, z = z, p = 2 * stats::pnorm(-abs(z)),
                                                    lo = r$effect - crit * r$se, hi = r$effect + crit * r$se))
                pm$addRow(rowKey = 2 * k - 1, values = list(estimand = elab[[e]], method = labs[[m]], group = "Y(1)", estimate = r$pom1, se = r$se_pom1))
                pm$addRow(rowKey = 2 * k, values = list(estimand = elab[[e]], method = labs[[m]], group = "Y(0)", estimate = r$pom0, se = r$se_pom0))
                if (r$moment_check > 1e-6 * max(1, stats::sd(q$y)))
                    notes <- c(notes, sprintf("%s, %s: the estimating equations are not solved to full precision (%s).", elab[[e]], labs[[m]], format(r$moment_check, digits = 2)))
            }
            tb$setNote("e", sprintf("Standard errors are M-estimation (sandwich) standard errors that include the estimation of the propensity score and the outcome regressions%s. Normal-based intervals. ATE: effect on the whole sample; ATT: on the treated; ATC: on the controls.",
                                    if (is.null(q$cl)) "" else ", clustered"))
            if (isTRUE(opts$showWeights)) {
                wt <- self$results$wdiag
                for (e in ests) {
                    ps <- mf_guard(mf_ps_fit_for(q$d, q$X, q$sw, opts$link, e))
                    w <- mf_weights(q$d, ps$p, e)
                    ws <- mf_weight_summary(w, q$d)
                    for (i in seq_len(nrow(ws))) {
                        kw <- kw + 1
                        wt$addRow(rowKey = kw, values = list(estimand = elab[[e]], group = ws$group[i], n = ws$n[i], min = ws$min[i], mean = ws$mean[i], max = ws$max[i],
                                                             ess = ws$ess[i], share = ws$share_max[i]))
                    }
                    if (any(ws$max / ws$mean > 20)) notes <- c(notes, sprintf("%s: some weights are more than 20 times the average: overlap is poor, IPW is unreliable here.", elab[[e]]))
                    if (any(ws$ess < 0.5 * ws$n)) notes <- c(notes, sprintf("%s: the effective sample size of a group is below half its size.", elab[[e]]))
                }
                wt$setNote("w", "Weights for each estimand; they are 1 for the group that is not reweighted. With a few very large weights the IPW estimate is unstable: consider trimming.")
            }
            if ("ato" %in% ests) notes <- c(notes, "ATO: the effect for the population that the overlap weights favour (units whose treatment is most uncertain); it needs no trimming, the weights are bounded, and with the logit model the weighted covariate means are balanced exactly.")
            if (identical(opts$link, "cbps")) notes <- c(notes, "CBPS: the propensity score solves the balancing conditions for each estimand (exactly identified), so the weighted covariate means of the groups match exactly; the standard errors use those conditions in place of the likelihood score.")
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (opts$trim != "none") notes <- c(notes, sprintf("Trimmed by the propensity score: %d treated and %d control observation(s) removed.", tr$n_trim[["treated"]], tr$n_trim[["control"]]))
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        }
    )
)
