# Nearest-neighbour matching with replacement (Abadie and Imbens 2006, 2011), as in Stata's teffects nnmatch:
# k nearest neighbours with all ties included, optional exact matching, caliper and bias adjustment, and the Abadie-Imbens
# heteroskedasticity-consistent standard error. Base R only.

# Coordinates in which the Euclidean distance equals the chosen metric.
mf_match_coords <- function(Xm, metric = c("maha", "ivar", "eucl")) {
    metric <- match.arg(metric)
    Xm <- as.matrix(Xm)
    switch(metric,
           eucl = Xm,
           ivar = { s <- apply(Xm, 2, stats::var); if (any(!is.finite(s) | s <= 0)) stop("A matching covariate has no variation.", call. = FALSE); sweep(Xm, 2, sqrt(s), "/") },
           maha = {
               S <- stats::cov(Xm)
               ch <- tryCatch(chol(S), error = function(e) stop("The covariates are perfectly collinear or constant, the Mahalanobis distance cannot be computed; remove a redundant covariate.", call. = FALSE))
               Xm %*% backsolve(ch, diag(ncol(Xm)))
           })
}

# For each unit in `from` (indices), the nearest `k` units of `pool` (indices) with ties included.
# Returns a list of integer vectors (empty when no candidate satisfies the exact and caliper conditions).
mf_nn <- function(Z, from, pool, k, key = NULL, caliper = NULL, self = FALSE) {
    tol <- 1e-9
    lapply(from, function(i) {
        cand <- if (self) pool else pool[pool != i]
        if (!is.null(key)) cand <- cand[key[cand] == key[i]]
        if (!length(cand)) return(integer(0))
        d2 <- colSums((t(Z[cand, , drop = FALSE]) - Z[i, ])^2)
        if (!is.null(caliper)) { keep <- d2 <= caliper^2 * (1 + tol); cand <- cand[keep]; d2 <- d2[keep] }
        if (!length(cand)) return(integer(0))
        thr <- sort(d2, partial = min(k, length(d2)))[min(k, length(d2))]
        cand[d2 <= thr * (1 + tol) + 1e-14]
    })
}

# estimand: "ate", "att" or "atc" (the last by exchanging the groups). y outcome, d treatment 0/1, Xm matching covariates,
# Zb design (with intercept) for the bias adjustment or NULL, exact: vector of group keys or NULL, caliper in distance units or NULL,
# k neighbours, J neighbours for the variance.
mf_nnmatch <- function(y, d, Xm, estimand = c("ate", "att", "atc"), k = 1, metric = c("maha", "ivar", "eucl"), caliper = NULL,
                       exact = NULL, Zb = NULL, J = 2) {
    estimand <- match.arg(estimand); metric <- match.arg(metric)
    flip <- estimand == "atc"
    if (flip) { d <- 1 - d; estimand <- "att" }
    n <- length(y)
    Z <- mf_match_coords(Xm, metric)
    key <- if (!is.null(exact)) as.integer(factor(exact)) else NULL
    i1 <- which(d == 1); i0 <- which(d == 0)
    # units that need a counterfactual: all (ATE) or the treated (ATT)
    from <- if (estimand == "ate") seq_len(n) else i1
    m_opp <- vector("list", n)
    for (g in c(1, 0)) {
        fg <- from[d[from] == g]
        if (length(fg)) m_opp[fg] <- mf_nn(Z, fg, if (g == 1) i0 else i1, k, key, caliper)
    }
    has <- vapply(from, function(i) length(m_opp[[i]]) > 0, logical(1))
    used <- from[has]; dropped <- from[!has]
    if (!length(used)) stop("No unit has a valid match; relax the caliper or the exact-matching variables.", call. = FALSE)
    if (estimand == "att" && !any(d[used] == 1)) stop("No treated unit has a valid match.", call. = FALSE)
    # number of times each unit serves as a match: sum of 1/m (m = number of matches of the unit that uses it) and sum of 1/m^2
    Kc <- numeric(n); Kp <- numeric(n)
    for (i in used) { j <- m_opp[[i]]; Kc[j] <- Kc[j] + 1 / length(j); Kp[j] <- Kp[j] + 1 / length(j)^2 }
    # outcome models for the bias adjustment, one per group, fitted on the units used as matches, weighted by Kc
    mu <- NULL
    if (!is.null(Zb)) {
        mu <- lapply(0:1, function(g) {
            if (estimand == "att" && g == 1) return(NULL)
            ig <- which(d == g & Kc > 0)
            Zg <- Zb[ig, , drop = FALSE]
            if (length(ig) < 2) stop("Too few units are used as matches to fit the bias-adjustment regression.", call. = FALSE)
            b <- stats::lm.wfit(Zg, y[ig], Kc[ig])$coefficients
            b[is.na(b)] <- 0                                      # collinear covariates among the matched units are dropped, as regress does
            drop(Zb %*% b)
        })
    }
    ycf <- rep(NA_real_, n)                                   # counterfactual (opposite arm) for the used units
    for (i in used) {
        j <- m_opp[[i]]
        ycf[i] <- mean(y[j]) + if (!is.null(mu)) { g <- 1 - d[i]; mu[[g + 1]][i] - mean(mu[[g + 1]][j]) } else 0
    }
    y1 <- ifelse(d == 1, y, ycf); y0 <- ifelse(d == 0, y, ycf)
    tau_i <- (y1 - y0)[used]
    nu <- length(used)
    tau <- mean(tau_i)
    # conditional variance of each unit: sample variance of the unit and its J nearest neighbours in the same group (ties included)
    sig2 <- rep(0, n)
    need_sig <- if (estimand == "ate") used else which(d == 0 & Kc > 0)
    for (g in c(1, 0)) {
        ig <- need_sig[d[need_sig] == g]
        if (!length(ig)) next
        nb <- mf_nn(Z, ig, which(d == g), J + 1, key, NULL, self = TRUE)
        for (a in seq_along(ig)) {
            jj <- nb[[a]]
            sig2[ig[a]] <- if (length(jj) > 1) sum((y[jj] - mean(y[jj]))^2) / (length(jj) - 1) else 0
        }
    }
    if (estimand == "ate") {
        v <- (sum((tau_i - tau)^2) + sum((Kc^2 + 2 * Kc - Kp)[used] * sig2[used])) / nu^2
    } else {
        v <- (sum((tau_i - tau)^2) + sum(((Kc^2 - Kp) * sig2)[d == 0])) / nu^2
    }
    w <- rep(0, n)
    if (estimand == "ate") { w[used] <- 1; w <- w + Kc } else { w[used] <- 1; w[d == 0] <- Kc[d == 0] }
    out <- list(effect = if (flip) -tau else tau, se = sqrt(v), n_used = nu, n_dropped = length(dropped), dropped = dropped,
                weights = w, K = Kc, ycf = ycf, matches = m_opp, n_matches = vapply(m_opp[used], length, integer(1)), flipped = flip, used = used)
    out
}

# EXPERIMENTAL, not used by any analysis. Propensity-score matching (teffects psmatch): point estimates agree with Stata to 1e-8 but the
# adjusted standard error does NOT (see PLAN.md); only the unadjusted standard error (se_unadjusted) is verified, against nnmatch on the score.
# Original note:  nearest neighbours on the estimated propensity score with the Abadie-Imbens (2016)
# standard error, which is adjusted because the score is estimated. X is the propensity-score design with an intercept, `Zc` the
# covariates without it, used to match on z for the derivative term of the ATET variance.
mf_psmatch <- function(y, d, X, estimand = c("ate", "att", "atc"), k = 1, link = c("logit", "probit"), caliper = NULL, J = 2) {
    estimand <- match.arg(estimand); link <- match.arg(link)
    n <- length(y); sw <- rep(1, n)
    fit <- mf_ps_fit(d, X, sw, link)
    p <- fit$p
    Vg <- mf_ps_vcov(d, X, sw, fit$coef, link)
    flip <- estimand == "atc"
    if (flip) { d <- 1 - d; p <- 1 - p; estimand <- "att" }
    eta <- drop(X %*% fit$coef)
    f <- if (link == "logit") stats::plogis(eta) * (1 - stats::plogis(eta)) else stats::dnorm(eta)
    base <- mf_nnmatch(y, if (flip) 1 - d else d, matrix(p, ncol = 1), if (flip) "atc" else estimand, k, "eucl", caliper, NULL, NULL, J)
    if (base$n_dropped > 0) stop("Some units have no match within the caliper; the propensity-score variance adjustment needs every unit to be matched. Widen the caliper.", call. = FALSE)
    tau <- if (flip) -base$effect else base$effect
    v0 <- base$se^2
    Z <- matrix(p, ncol = 1)
    same <- lapply(seq_len(n), function(i) mf_nn(Z, i, which(d == d[i]), J + 1, NULL, NULL, self = TRUE)[[1]])
    opp <- lapply(seq_len(n), function(i) mf_nn(Z, i, which(d != d[i]), J, NULL, NULL)[[1]])
    Xc <- X[, -1, drop = FALSE]
    covz <- function(S) {
        if (length(S) < 2) return(numeric(ncol(Xc)))
        zc <- sweep(Xc[S, , drop = FALSE], 2, colMeans(Xc[S, , drop = FALSE])); yc <- y[S] - mean(y[S])
        drop(crossprod(zc, yc)) / (length(S) - 1)
    }
    cov1 <- t(vapply(seq_len(n), function(i) covz(if (d[i] == 1) same[[i]] else opp[[i]]), numeric(ncol(Xc))))   # cov(z_i, y_i(1))
    cov0 <- t(vapply(seq_len(n), function(i) covz(if (d[i] == 0) same[[i]] else opp[[i]]), numeric(ncol(Xc))))   # cov(z_i, y_i(0))
    pad <- function(M) cbind(0, M)                                                                                  # no covariance for the intercept
    if (estimand == "ate") {
        ct <- colMeans(f * (pad(cov1) / p + pad(cov0) / (1 - p)))
        v <- v0 + drop(t(ct) %*% Vg %*% ct)
        tauhat <- tau
    } else {
        n1 <- sum(d == 1)
        ysame <- vapply(seq_len(n), function(i) { s <- setdiff(mf_nn(Z, i, which(d == d[i]), J + 1, NULL, NULL, self = TRUE)[[1]], i); mean(y[s]) }, numeric(1))
        yopp <- vapply(opp, function(s) mean(y[s]), numeric(1))
        y1t <- ifelse(d == 1, ysame, yopp); y0t <- ifelse(d == 0, ysame, yopp)
        delta <- if (flip) -tau else tau
        c1 <- colSums(X * (f * (y1t - y0t - delta))) / n1
        c2 <- colSums(f * (pad(cov1) + pad(cov0) * (p / (1 - p)))) / n1
        cd <- c1 + c2
        # derivative term: opposite-group neighbours found on the covariates z
        Zm <- mf_match_coords(Xc, "maha")
        yz <- vapply(seq_len(n), function(i) mean(y[mf_nn(Zm, i, which(d != d[i]), k, NULL, NULL)[[1]]]), numeric(1))
        dd <- colSums(X * (f * ((2 * d - 1) * (y - yz) - delta))) / n1
        v <- v0 - drop(t(cd) %*% Vg %*% cd) + drop(t(dd) %*% Vg %*% dd)
    }
    c(base[c("n_used", "n_dropped", "weights", "K", "matches", "n_matches", "used")], list(effect = base$effect, se = sqrt(max(v, 0)), se_unadjusted = base$se, ps = if (flip) 1 - p else p, link = link))
}
