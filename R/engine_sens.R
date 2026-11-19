# Sensitivity to unobserved confounding: E-value (VanderWeele and Ding 2017) and Rosenbaum bounds for matched pairs. Base R only.

# Risk-ratio equivalent of an effect measure. `type`: smd (standardised mean difference), rr, or_rare, or_common, hr_rare, hr_common.
mf_to_rr <- function(x, type) {
    if (!is.finite(x)) return(NA_real_)
    switch(type,
           smd = exp(0.91 * x),
           rr = x, or_rare = x, hr_rare = x,
           or_common = sqrt(x),
           hr_common = (1 - 0.5^sqrt(x)) / (1 - 0.5^sqrt(1 / x)),
           stop("Unknown effect measure.", call. = FALSE))
}

# E-value for a risk ratio: the minimum strength of association (on the risk-ratio scale) that an unmeasured confounder would need
# with both treatment and outcome to explain the estimate away.
mf_evalue_rr <- function(rr) {
    if (!is.finite(rr) || rr <= 0) return(NA_real_)
    if (rr < 1) rr <- 1 / rr
    rr + sqrt(rr * (rr - 1))
}

# Estimate with its confidence limits. Returns the risk-ratio equivalents and the E-values of the estimate and of the confidence limit
# closest to the null (1 when the interval includes the null).
mf_evalue <- function(est, lo, hi, type = c("smd", "rr", "or_rare", "or_common", "hr_rare", "hr_common")) {
    type <- match.arg(type)
    null <- if (type == "smd") 0 else 1
    if (type != "smd" && (est <= 0 || lo <= 0 || hi <= 0)) stop("Ratio measures must be positive.", call. = FALSE)
    if (lo > hi) { t <- lo; lo <- hi; hi <- t }
    rr <- mf_to_rr(est, type); rl <- mf_to_rr(lo, type); rh <- mf_to_rr(hi, type)
    cross <- lo <= null && hi >= null
    e_ci <- if (cross) 1 else mf_evalue_rr(if (est > null) rl else rh)
    list(rr = rr, rr_lo = rl, rr_hi = rh, e_point = mf_evalue_rr(rr), e_ci = e_ci, ci_includes_null = cross)
}

# Rosenbaum bounds for paired differences d (treated minus matched control). Under hidden bias of size Gamma the probability that the
# treated unit of a pair has the higher outcome lies between 1/(1+Gamma) and Gamma/(1+Gamma). Returns upper and lower bounds of the
# one-sided p-value of the sign test (exact binomial) and of the Wilcoxon signed-rank test (normal approximation, mid-ranks).
mf_rosenbaum <- function(d, gammas) {
    d <- d[is.finite(d) & d != 0]
    n <- length(d)
    if (n < 5) stop("Fewer than 5 pairs with a non-zero difference; Rosenbaum bounds cannot be computed.", call. = FALSE)
    s <- sum(d > 0)
    r <- rank(abs(d)); tstat <- sum(r[d > 0]); sr <- sum(r); sr2 <- sum(r^2)
    rows <- lapply(gammas, function(g) {
        pu <- g / (1 + g); pl <- 1 / (1 + g)
        sign_u <- stats::pbinom(s - 1, n, pu, lower.tail = FALSE); sign_l <- stats::pbinom(s - 1, n, pl, lower.tail = FALSE)
        wil <- function(p) stats::pnorm((tstat - p * sr) / sqrt(p * (1 - p) * sr2), lower.tail = FALSE)
        data.frame(gamma = g, sign_lo = sign_l, sign_hi = sign_u, wilcoxon_lo = wil(pl), wilcoxon_hi = wil(pu))
    })
    list(table = do.call(rbind, rows), n = n, n_pos = s, statistic = tstat)
}

# Smallest Gamma at which the upper bound of the p-value exceeds `alpha`, by bisection; Inf when it never does up to `gmax`.
mf_critical_gamma <- function(d, test = c("wilcoxon", "sign"), alpha = 0.05, gmax = 50) {
    test <- match.arg(test)
    f <- function(g) { tb <- mf_rosenbaum(d, g)$table; (if (test == "sign") tb$sign_hi else tb$wilcoxon_hi) - alpha }
    if (f(1) > 0) return(1)
    if (f(gmax) < 0) return(Inf)
    stats::uniroot(f, c(1, gmax), tol = 1e-8)$root
}
