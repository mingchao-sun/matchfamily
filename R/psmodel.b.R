# jamovi analysis class: propensity score model, overlap and trimming.

psmodelClass <- R6::R6Class(
    "psmodelClass",
    inherit = psmodelBase,
    private = list(

        .run = function() {
            opts <- self$options
            mf_show_text(self$results, opts, "psmodel")
            if (is.null(opts$treat) || length(opts$covs) == 0) {
                self$results$notes$setContent("<p>Select a binary treatment and at least one covariate.</p>")
                return()
            }
            pp <- mf_guard(mf_prepare(self$data, opts$treat, opts$covs, weights = opts$weights))
            fit <- mf_guard(mf_ps_fit(pp$d, pp$X, pp$sw, opts$link))
            V <- mf_guard(mf_ps_vcov(pp$d, pp$X, pp$sw, fit$coef, opts$link))
            se <- sqrt(diag(V)); z <- fit$coef / se
            crit <- stats::qnorm(1 - (1 - opts$level / 100) / 2)
            tb <- self$results$coef
            for (i in seq_along(fit$coef))
                tb$addRow(rowKey = i, values = list(term = names(fit$coef)[i], estimate = fit$coef[[i]], se = se[i], z = z[i],
                                                    p = 2 * stats::pnorm(-abs(z[i])), lo = fit$coef[[i]] - crit * se[i], hi = fit$coef[[i]] + crit * se[i]))
            tb$setNote("l", sprintf("%s model of the probability of treatment (%s is the treated group). Robust standard errors.",
                                    if (opts$link == "logit") "Logit" else "Probit", if (is.null(pp$labels)) "1" else pp$labels[2]))

            p <- fit$p; d <- pp$d
            ot <- self$results$overlap
            for (g in c(1, 0)) {
                x <- p[d == g]; q <- stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
                ot$addRow(rowKey = g + 1, values = list(group = if (g == 1) "Treated" else "Control", n = length(x), min = min(x), q05 = q[1], median = q[2], q95 = q[3], max = max(x)))
            }
            a <- max(min(p[d == 1]), min(p[d == 0])); b <- min(max(p[d == 1]), max(p[d == 0]))
            out_t <- sum(d == 1 & (p < a | p > b)); out_c <- sum(d == 0 & (p < a | p > b))
            keep <- mf_trim(p, d, opts$trim, opts$trimLo, opts$trimHi)
            items <- list(
                c("Observations", format(length(d), big.mark = ",")),
                c("Treated / control", sprintf("%s / %s", format(sum(d == 1), big.mark = ","), format(sum(d == 0), big.mark = ","))),
                c("c statistic (area under the ROC curve)", mf_fmt(mf_auc(p, d), 3)),
                c("Common support", sprintf("[%s, %s]", mf_fmt(a, 3), mf_fmt(b, 3))),
                c("Outside the common support (treated / control)", sprintf("%d / %d", out_t, out_c)))
            if (opts$trim != "none")
                items[[length(items) + 1]] <- c("Removed by the trimming rule (treated / control)", sprintf("%d / %d", sum(!keep & d == 1), sum(!keep & d == 0)))
            sm <- self$results$summary
            for (i in seq_along(items)) sm$addRow(rowKey = i, values = list(item = items[[i]][1], value = items[[i]][2]))

            notes <- character(0)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (a >= b) notes <- c(notes, "The two groups have no common support: there is no overlap in the estimated propensity score.")
            else if ((out_t + out_c) / length(d) > 0.2) notes <- c(notes, "More than 20% of the observations lie outside the common support; overlap is poor and estimates will rely on extrapolation.")
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
            if (isTRUE(opts$plot))
                self$results$psPlot$setState(list(p = p, d = d, a = a, b = b, trim = opts$trim, lo = opts$trimLo, hi = opts$trimHi))
        },

        .plotPs = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- data.frame(p = st$p, group = factor(ifelse(st$d == 1, "Treated", "Control"), levels = c("Control", "Treated")))
            g <- ggplot2::ggplot(df, ggplot2::aes(x = p, fill = group, colour = group)) +
                ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30, alpha = 0.35, position = "identity") +
                ggplot2::labs(x = "Estimated propensity score", y = "Density", fill = NULL, colour = NULL) +
                ggplot2::coord_cartesian(xlim = c(0, 1)) + ggtheme
            if (st$a < st$b) g <- g + ggplot2::geom_vline(xintercept = c(st$a, st$b), linetype = "dashed")
            print(g)
            TRUE
        }
    )
)
