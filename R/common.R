# Helpers shared by the matching and weighting analyses.

mf_fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, drop0trailing = TRUE)

# Turn engine errors into readable jamovi messages.
mf_guard <- function(expr) tryCatch(expr, error = function(e) jmvcore::reject(conditionMessage(e)))

mf_show_text <- function(results, opts, key) {
    if (isTRUE(opts$showGuide)) results$guide$setContent(mf_guide_html(key))
    if (isTRUE(opts$showHelp)) results$helpText$setContent(mf_notes_html(key))
}

# Estimation sample after the optional trimming on the propensity score. Returns the prepared data restricted to
# the retained observations, the propensity score fitted on all observations and the number trimmed in each group.
mf_trimmed <- function(pp, link, rule, lo, hi) {
    ps0 <- mf_ps_fit(pp$d, pp$X, pp$sw, link)
    keep <- mf_trim(ps0$p, pp$d, rule, lo, hi)
    if (sum(keep & pp$d == 1) < 5 || sum(keep & pp$d == 0) < 5)
        stop("Trimming leaves fewer than 5 observations in a treatment group; relax the trimming rule.", call. = FALSE)
    out <- pp
    out$y <- if (!is.null(pp$y)) pp$y[keep] else NULL
    out$d <- pp$d[keep]; out$X <- pp$X[keep, , drop = FALSE]; out$W <- pp$W[keep, , drop = FALSE]; out$sw <- pp$sw[keep]
    out$cl <- if (!is.null(pp$cl)) pp$cl[keep] else NULL
    list(pp = out, ps_all = ps0, keep = keep, n_trim = c(treated = sum(!keep & pp$d == 1), control = sum(!keep & pp$d == 0)))
}
