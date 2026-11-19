# jamovi analysis class: comparison of the estimators of the matching and weighting family.

matchcmpClass <- R6::R6Class(
    "matchcmpClass",
    inherit = matchcmpBase,
    private = list(

        .run = function() {
            opts <- self$options
            mf_show_text(self$results, opts, "matchcmp")
            if (is.null(opts$outcome) || is.null(opts$treat) || length(opts$covs) == 0) {
                self$results$notes$setContent("<p>Select an outcome, a binary treatment and at least one covariate.</p>")
                return()
            }
            flags <- c(naive = opts$doNaive, ra = opts$doRa, ipw = opts$doIpw, aipw = opts$doAipw, ipwra = opts$doIpwra, nn = opts$doNn, ebal = opts$doEbal)
            methods <- names(flags)[vapply(flags, isTRUE, logical(1))]
            if (!length(methods)) jmvcore::reject("Tick at least one estimator.")
            pp <- mf_guard(mf_prepare(self$data, opts$treat, opts$covs, outcome = opts$outcome))
            labs <- mf_compare_methods()
            elab <- c(ate = "ATE", att = "ATT", atc = "ATC")
            ests <- if (opts$estimand == "all") c("ate", "att", "atc") else opts$estimand
            crit <- stats::qnorm(1 - (1 - opts$level / 100) / 2)
            tb <- self$results$est
            notes <- character(0); plotdf <- NULL; k <- 0
            for (e in ests) {
                res <- mf_guard(mf_compare(pp, opts$covs, methods, e, opts$link, opts$k, isTRUE(opts$biasadj), as.numeric(substring(opts$moments, 2))))
                for (i in seq_along(res)) {
                    r <- res[[i]]; k <- k + 1
                    if (!is.null(r$error)) {
                        tb$addRow(rowKey = k, values = list(estimand = elab[[e]], method = labs[[r$method]]))
                        notes <- c(notes, sprintf("%s, %s could not be computed: %s", elab[[e]], labs[[r$method]], r$error))
                        next
                    }
                    z <- r$effect / r$se
                    tb$addRow(rowKey = k, values = list(estimand = elab[[e]], method = labs[[r$method]], estimate = r$effect, se = r$se, p = 2 * stats::pnorm(-abs(z)),
                                                        lo = r$effect - crit * r$se, hi = r$effect + crit * r$se, smd = r$smd_max, above = r$above))
                    plotdf <- rbind(plotdf, data.frame(estimand = elab[[e]], method = labs[[r$method]], est = r$effect, lo = r$effect - crit * r$se, hi = r$effect + crit * r$se,
                                                       stringsAsFactors = FALSE))
                    if (!is.null(r$dropped) && r$dropped > 0) notes <- c(notes, sprintf("%s, matching: %d unit(s) without a valid match were left out.", elab[[e]], r$dropped))
                }
            }
            tb$setNote("e", "The unadjusted difference does not depend on the estimand. Balance columns are measured after the weighting or matching of each method; regression-based estimators (RA, AIPW, IPWRA) do not reweight the sample and have none. Normal-based intervals.")
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
            if (isTRUE(opts$plot) && !is.null(plotdf)) self$results$forest$setState(plotdf)
        },

        .plotForest = function(image, ggtheme, theme, ...) {
            df <- image$state
            if (is.null(df)) return(FALSE)
            df$method <- factor(df$method, levels = rev(unique(df$method)))
            g <- ggplot2::ggplot(df, ggplot2::aes(x = est, y = method)) +
                ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
                ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.2) +
                ggplot2::geom_point(size = 3) +
                ggplot2::labs(x = "Effect estimate", y = NULL) + ggtheme
            if (length(unique(df$estimand)) > 1) g <- g + ggplot2::facet_wrap(~ estimand, ncol = 1)
            print(g)
            TRUE
        }
    )
)
