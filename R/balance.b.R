# jamovi analysis class: covariate balance before and after weighting.

balanceClass <- R6::R6Class(
    "balanceClass",
    inherit = balanceBase,
    private = list(

        .run = function() {
            opts <- self$options
            mf_show_text(self$results, opts, "balance")
            if (is.null(opts$treat) || length(opts$covs) == 0) {
                self$results$notes$setContent("<p>Select a binary treatment and at least one covariate.</p>")
                return()
            }
            adjust <- opts$adjust
            if (adjust == "var" && is.null(opts$adjW)) {
                self$results$notes$setContent("<p>Select the weights variable.</p>")
                return()
            }
            pp <- mf_guard(mf_prepare(self$data, opts$treat, opts$covs, weights = opts$weights, aux = if (adjust == "var") opts$adjW else NULL))
            d <- pp$d; sw <- pp$sw
            bm <- mf_bal_matrix(pp$data, opts$covs)
            notes <- character(0)
            w_adj <- NULL
            if (adjust == "ipw") {
                ps <- mf_guard(mf_ps_fit_for(d, pp$X, sw, opts$link, opts$estimand))
                w_adj <- sw * mf_weights(d, ps$p, opts$estimand)
            } else if (adjust == "ebal") {
                Xc <- mf_design(pp$data, opts$covs, intercept = FALSE)
                w_adj <- mf_guard(mf_ebal_weights(d, Xc, sw, opts$estimand, as.numeric(substring(opts$moments, 2))))$w
            } else if (adjust == "var") {
                a <- pp$aux
                if (any(a < 0)) jmvcore::reject("The weights must not be negative.")
                if (sum(a * (d == 1)) <= 0 || sum(a * (d == 0)) <= 0) jmvcore::reject("The weights are zero for every unit of a treatment group.")
                w_adj <- sw * a
            }
            b <- mf_guard(mf_balance(bm$M, bm$binary, d, sw, w_adj, opts$estimand, opts$denom, isTRUE(opts$stdBinary)))
            thr <- opts$threshold
            tb <- self$results$bal
            for (i in seq_len(nrow(b$raw))) {
                r <- b$raw[i, ]
                v <- list(variable = r$variable, m1 = r$mean1, m0 = r$mean0, smd = r$smd, vr = r$vr, ks = r$ks)
                if (!is.null(b$adj)) {
                    a <- b$adj[i, ]
                    v <- c(v, list(am1 = a$mean1, am0 = a$mean0, asmd = a$smd, avr = a$vr, aks = a$ks))
                }
                tb$addRow(rowKey = i, values = v)
            }
            dn <- c(pooled = "the pooled standard deviation of the two groups", treated = "the standard deviation of the treated group",
                    control = "the standard deviation of the control group")[[b$denom]]
            tb$setNote("s", sprintf("Standardised difference = difference in means divided by %s in the unadjusted sample%s. Variance ratio = treated / control. KS = largest distance between the two distribution functions.",
                                    dn, if (isTRUE(opts$stdBinary)) "" else " (binary variables: difference in proportions)"))

            sb <- mf_balance_summary(b$raw, thr); sa <- if (!is.null(b$adj)) mf_balance_summary(b$adj, thr) else NULL
            ess <- function(w) c(mf_ess(w[d == 1]), mf_ess(w[d == 0]))
            e0 <- ess(sw); e1 <- if (!is.null(w_adj)) ess(w_adj) else NULL
            f2 <- function(x) mf_fmt(x, 3)
            items <- list(
                list("Observations", format(length(d), big.mark = ","), format(length(d), big.mark = ",")),
                list("Treated / control", sprintf("%s / %s", sum(d == 1), sum(d == 0)), sprintf("%s / %s", sum(d == 1), sum(d == 0))),
                list("Largest absolute standardised difference", f2(sb[["max"]]), if (!is.null(sa)) f2(sa[["max"]]) else ""),
                list("Mean absolute standardised difference", f2(sb[["mean"]]), if (!is.null(sa)) f2(sa[["mean"]]) else ""),
                list(sprintf("Variables above %s", f2(thr)), sprintf("%d of %d", sb[["above"]], sb[["n"]]), if (!is.null(sa)) sprintf("%d of %d", sa[["above"]], sa[["n"]]) else ""),
                list("Effective sample size (treated / control)", sprintf("%s / %s", f2(e0[1]), f2(e0[2])), if (!is.null(e1)) sprintf("%s / %s", f2(e1[1]), f2(e1[2])) else ""))
            sm <- self$results$summary
            for (i in seq_along(items)) sm$addRow(rowKey = i, values = list(item = items[[i]][[1]], before = items[[i]][[2]], after = items[[i]][[3]]))
            if (!is.null(sa) && sa[["above"]] > 0)
                notes <- c(notes, sprintf("%d variable(s) remain above the threshold after adjustment: the adjustment has not balanced them.", sa[["above"]]))
            if (!is.null(e1) && any(e1 < 0.5 * c(sum(d == 1), sum(d == 0))))
                notes <- c(notes, "The effective sample size of a group fell below half its size: a few units carry much of the weight.")
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (any(is.na(b$raw$smd))) notes <- c(notes, "A variable has no variation in the reference group, its standardised difference is not defined.")
            self$results$notes$setContent(if (length(notes)) paste0("<p>", paste(notes, collapse = " "), "</p>") else "")
            if (isTRUE(opts$plot))
                self$results$lovePlot$setState(list(raw = b$raw, adj = b$adj, thr = thr))
        },

        .plotLove = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            lv <- st$raw$variable[order(abs(st$raw$smd), na.last = FALSE)]
            df <- data.frame(variable = factor(st$raw$variable, levels = lv), smd = abs(st$raw$smd), stage = "Unadjusted")
            if (!is.null(st$adj)) df <- rbind(df, data.frame(variable = factor(st$adj$variable, levels = lv), smd = abs(st$adj$smd), stage = "Adjusted"))
            df$stage <- factor(df$stage, levels = c("Unadjusted", "Adjusted"))
            g <- ggplot2::ggplot(df, ggplot2::aes(x = smd, y = variable, shape = stage, colour = stage)) +
                ggplot2::geom_vline(xintercept = st$thr, linetype = "dashed") +
                ggplot2::geom_point(size = 3) +
                ggplot2::labs(x = "Absolute standardised difference", y = NULL, shape = NULL, colour = NULL) + ggtheme
            print(g)
            TRUE
        }
    )
)
