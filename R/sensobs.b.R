# jamovi analysis class: sensitivity to unobserved confounding (E-value, Rosenbaum bounds).

sensobsClass <- R6::R6Class(
    "sensobsClass",
    inherit = sensobsBase,
    private = list(

        .run = function() {
            opts <- self$options
            mf_show_text(self$results, opts, "sensobs")
            crit <- stats::qnorm(1 - (1 - opts$level / 100) / 2)
            notes <- character(0)
            ev <- self$results$evalue
            pp <- NULL
            if (opts$mode == "manual") {
                if (!all(is.finite(c(opts$mEst, opts$mLo, opts$mHi)))) jmvcore::reject("Enter the estimate and both confidence limits.")
                r <- mf_guard(mf_evalue(opts$mEst, opts$mLo, opts$mHi, opts$mType))
                est <- opts$mEst; lo <- min(opts$mLo, opts$mHi); hi <- max(opts$mLo, opts$mHi); type <- opts$mType
            } else {
                if (is.null(opts$outcome) || is.null(opts$treat) || length(opts$covs) == 0) {
                    self$results$notes$setContent("<p>Select an outcome, a binary treatment and at least one covariate, or switch to manual input.</p>")
                    return()
                }
                pp <- mf_guard(mf_prepare(self$data, opts$treat, opts$covs, outcome = opts$outcome))
                fit <- mf_guard(mf_estimate(pp$y, pp$d, pp$X, pp$X, pp$sw, NULL, "aipw", opts$estimand, opts$link))
                s1 <- stats::sd(pp$y[pp$d == 1]); s0 <- stats::sd(pp$y[pp$d == 0]); sdp <- sqrt((s1^2 + s0^2) / 2)
                if (!is.finite(sdp) || sdp <= 0) jmvcore::reject("The outcome has no variation.")
                est <- fit$effect / sdp; lo <- (fit$effect - crit * fit$se) / sdp; hi <- (fit$effect + crit * fit$se) / sdp; type <- "smd"
                r <- mf_guard(mf_evalue(est, lo, hi, "smd"))
                self$results$effect$addRow(rowKey = 1, values = list(method = c(ate = "ATE", att = "ATT", atc = "ATC")[[opts$estimand]], estimate = fit$effect, se = fit$se,
                                                                     lo = fit$effect - crit * fit$se, hi = fit$effect + crit * fit$se, d = est, dlo = lo, dhi = hi))
                self$results$effect$setNote("a", "Doubly robust (AIPW) estimate; the standardised difference divides by the pooled within-group standard deviation of the outcome.")
            }
            null <- if (type == "smd") 0 else 1
            near <- if (r$ci_includes_null) NA_real_ else if (est > null) lo else hi
            ev$addRow(rowKey = 1, values = list(item = "Estimate", value = est, rr = r$rr, e = r$e_point))
            ev$addRow(rowKey = 2, values = list(item = "Confidence limit closest to the null", value = if (is.na(near)) NaN else near,
                                                rr = if (is.na(near)) NaN else mf_to_rr(near, type), e = r$e_ci))
            conv <- if (type == "smd") "A standardised mean difference d is converted to a risk ratio by exp(0.91 d), an approximation that holds best for common outcomes. "
                    else if (type %in% c("or_common", "hr_common")) "The ratio is converted to an approximate risk ratio for a common outcome. " else ""
            ev$setNote("e", paste0("The E-value is the smallest association, on the risk-ratio scale, that an unmeasured confounder would need with both the treatment and the outcome (beyond the measured covariates) to explain the estimate away. ",
                                   conv, if (r$ci_includes_null) "The confidence interval includes the null, so the E-value of the limit is 1: the result is already compatible with no effect." else ""))
            if (opts$mode == "data" && isTRUE(opts$doRosenbaum)) {
                Xc <- mf_design(pp$data, opts$covs, intercept = FALSE)
                m <- mf_guard(mf_nnmatch(pp$y, pp$d, Xc, "att", opts$k, "maha", NULL, NULL, NULL, 2))
                tr <- m$used[pp$d[m$used] == 1]
                dif <- pp$y[tr] - m$ycf[tr]
                sgn <- if (sum(dif < 0) > sum(dif > 0)) -1 else 1
                dif <- sgn * dif
                gam <- seq(1, opts$gammaMax, by = opts$gammaStep)
                rb <- mf_guard(mf_rosenbaum(dif, gam))
                rt <- self$results$rosen
                for (i in seq_len(nrow(rb$table))) {
                    q <- rb$table[i, ]
                    rt$addRow(rowKey = i, values = list(gamma = q$gamma, wlo = q$wilcoxon_lo, whi = q$wilcoxon_hi, slo = q$sign_lo, shi = q$sign_hi))
                }
                rt$setNote("r", sprintf("One-sided tests of a %s effect among the treated. Gamma is the factor by which unobserved confounding could change the odds of treatment between two matched units. The upper bound is the worst case; at Gamma = 1 both bounds equal the usual p-value.",
                                        if (sgn > 0) "positive" else "negative"))
                cg_w <- mf_critical_gamma(dif, "wilcoxon", 0.05); cg_s <- mf_critical_gamma(dif, "sign", 0.05)
                fg <- function(x) if (is.infinite(x)) "not reached" else mf_fmt(x, 2)
                st <- self$results$rsum
                items <- list(c("Treated units matched", format(length(tr), big.mark = ",")), c("Pairs with a non-zero difference", format(rb$n, big.mark = ",")),
                              c("Pairs in the direction of the effect", sprintf("%d of %d", if (sgn > 0) rb$n_pos else rb$n - rb$n_pos, rb$n)),
                              c("Gamma at which the signed-rank upper bound passes 0.05", fg(cg_w)), c("Gamma at which the sign-test upper bound passes 0.05", fg(cg_s)))
                for (i in seq_along(items)) st$addRow(rowKey = i, values = list(item = items[[i]][1], value = items[[i]][2]))
                if (isTRUE(opts$plot)) self$results$rosPlot$setState(list(tb = rb$table))
                if (m$n_dropped > 0) notes <- c(notes, sprintf("%d treated unit(s) had no valid match.", m$n_dropped))
                notes <- c(notes, "The matched pairs use the ATT with replacement; a control used by several treated units appears in several pairs, which the bounds do not account for.")
            }
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
        },

        .plotRos = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- rbind(data.frame(gamma = st$tb$gamma, p = st$tb$wilcoxon_hi, test = "Signed rank"), data.frame(gamma = st$tb$gamma, p = st$tb$sign_hi, test = "Sign test"))
            g <- ggplot2::ggplot(df, ggplot2::aes(x = gamma, y = p, colour = test, linetype = test)) +
                ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed") + ggplot2::geom_line() + ggplot2::geom_point() +
                ggplot2::labs(x = "Gamma (hidden bias)", y = "Upper bound of the p-value", colour = NULL, linetype = NULL) + ggtheme
            print(g)
            TRUE
        }
    )
)
