# Covariate balance: standardised mean difference, variance ratio and Kolmogorov-Smirnov statistic, before and after weighting.
# Base R only.

# Covariate columns for the balance table: numeric variables as they are, a dummy for every level of a factor or text variable.
# `binary` marks columns that only take the values 0 and 1.
mf_bal_matrix <- function(data, covs) {
    cols <- list(); names_ <- character(0); bin <- logical(0)
    for (v in covs) {
        x <- data[[v]]
        if (is.character(x) || is.logical(x)) x <- factor(x)
        if (is.factor(x)) {
            x <- droplevels(x)
            for (lv in levels(x)) { cols[[length(cols) + 1]] <- as.numeric(x == lv); names_ <- c(names_, paste0(v, ": ", lv)); bin <- c(bin, TRUE) }
        } else {
            x <- as.numeric(x)
            cols[[length(cols) + 1]] <- x; names_ <- c(names_, v); bin <- c(bin, all(x %in% c(0, 1)))
        }
    }
    M <- do.call(cbind, cols); colnames(M) <- names_
    list(M = M, binary = bin)
}

mf_wmean <- function(x, w) sum(w * x) / sum(w)

# Weighted variance with reliability-weight (frequency-like) correction; equals var() for equal weights.
mf_wvar <- function(x, w) {
    m <- mf_wmean(x, w); den <- sum(w) - sum(w^2) / sum(w)
    if (den <= 0) return(NA_real_)
    sum(w * (x - m)^2) / den
}

# Weighted Kolmogorov-Smirnov distance between the two groups.
mf_wks <- function(x, d, w) {
    o <- order(x); x <- x[o]; d <- d[o]; w <- w[o]
    w1 <- w * (d == 1); w0 <- w * (d == 0)
    last <- c(x[-1] != x[-length(x)], TRUE)                  # evaluate the distribution functions after the last tie
    f1 <- cumsum(w1) / sum(w1); f0 <- cumsum(w0) / sum(w0)
    max(abs(f1[last] - f0[last]))
}

# Balance table. `w_adj` are the adjusting weights (already multiplied by the sampling weights); `sw` are the sampling weights
# used for the unadjusted comparison. The standardising denominator always comes from the unadjusted sample.
mf_balance <- function(M, binary, d, sw, w_adj = NULL, estimand = c("ate", "att", "atc", "ato"), denom = c("auto", "pooled", "treated", "control"),
                       std_binary = TRUE) {
    estimand <- match.arg(estimand); denom <- match.arg(denom)
    if (denom == "auto") denom <- c(ate = "pooled", att = "treated", atc = "control", ato = "pooled")[[estimand]]
    one <- function(w) {
        t(vapply(seq_len(ncol(M)), function(j) {
            x <- M[, j]; i1 <- d == 1; i0 <- d == 0
            m1 <- mf_wmean(x[i1], w[i1]); m0 <- mf_wmean(x[i0], w[i0])
            v1 <- mf_wvar(x[i1], w[i1]); v0 <- mf_wvar(x[i0], w[i0])
            c(m1 = m1, m0 = m0, v1 = v1, v0 = v0, ks = mf_wks(x, d, w))
        }, numeric(5)))
    }
    raw <- one(sw)
    sd_den <- apply(raw, 1, function(r) switch(denom, pooled = sqrt((r["v1"] + r["v0"]) / 2), treated = sqrt(r["v1"]), control = sqrt(r["v0"])))
    sd_den[!is.finite(sd_den) | sd_den <= 0] <- NA
    if (!std_binary) sd_den[binary] <- 1
    fin <- function(tb) data.frame(variable = colnames(M), mean1 = tb[, "m1"], mean0 = tb[, "m0"], smd = (tb[, "m1"] - tb[, "m0"]) / sd_den,
                                   vr = tb[, "v1"] / tb[, "v0"], ks = tb[, "ks"], binary = binary, stringsAsFactors = FALSE)
    out <- list(raw = fin(raw), denom = denom)
    if (!is.null(w_adj)) out$adj <- fin(one(w_adj))
    out
}

# Counts used in the summary table.
mf_balance_summary <- function(tb, threshold) {
    a <- abs(tb$smd)
    c(max = max(a, na.rm = TRUE), mean = mean(a, na.rm = TRUE), above = sum(a > threshold, na.rm = TRUE), n = nrow(tb))
}
