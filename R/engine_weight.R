# Regression adjustment, inverse probability weighting, augmented IPW and IPW regression adjustment.
# Every estimator is a set of moment conditions in the parameters theta = (propensity score, outcome regressions, potential-outcome
# means); the covariance is the M-estimation sandwich, so the uncertainty of the estimated propensity score and of the
# outcome regressions is taken into account (as Stata's teffects does).

# Centre and scale the non-constant columns of a design matrix. The fitted values, effects and standard errors do not change
# (the models are linear re-parametrisations), but the numerical derivatives and the matrix inversions become well conditioned.
mf_scale <- function(M) {
    for (j in seq_len(ncol(M))) {
        s <- stats::sd(M[, j])
        if (is.finite(s) && s > 0) M[, j] <- (M[, j] - mean(M[, j])) / s
    }
    M
}

mf_wls <- function(W, y, w) {
    ok <- w > 0
    if (sum(ok) < ncol(W)) stop("Too few observations to fit an outcome regression in one treatment group.", call. = FALSE)
    fit <- stats::lm.wfit(W[ok, , drop = FALSE], y[ok], w[ok])
    if (anyNA(fit$coefficients)) stop("An outcome regression has a collinear covariate in one treatment group; remove it or merge sparse categories.", call. = FALSE)
    fit$coefficients
}

# Estimate one estimator for one estimand. d is 0/1, X the propensity score design, W the outcome-model design (both with
# an intercept), sw sampling weights, cl optional cluster ids.
mf_estimate <- function(y, d, X, W, sw = rep(1, length(y)), cl = NULL, method = c("ra", "ipw", "aipw", "ipwra"),
                        estimand = c("ate", "att", "atc", "ato"), link = c("logit", "probit", "cbps")) {
    method <- match.arg(method); estimand <- match.arg(estimand); link <- match.arg(link)
    if (estimand == "ato" && method == "ipwra") stop("Overlap weights are available for regression adjustment, IPW and AIPW, not for IPWRA.", call. = FALSE)
    if (estimand == "ato" && link == "cbps") link <- "logit"          # the logit score already balances the covariates exactly under overlap weights
    flip <- estimand == "atc"
    if (flip) { d <- 1 - d; estimand <- "att" }
    X <- mf_scale(X); W <- mf_scale(W)
    ysc <- stats::sd(y); if (!is.finite(ysc) || ysc <= 0) ysc <- 1
    ymu <- mean(y); y <- (y - ymu) / ysc                    # outcome on a unit scale; the effect is multiplied back below
    n <- length(y); kx <- ncol(X); kw <- ncol(W)
    ps <- method != "ra" || estimand == "ato"
    # ---- point estimates ----------------------------------------------------------------------------------------------
    gamma <- NULL; p <- NULL
    if (ps) { fit <- mf_ps_fit(d, X, sw, link, estimand); gamma <- fit$coef; p <- fit$p }
    ipw1 <- if (ps) d / p else NULL
    ipw0 <- if (ps) (if (estimand == "ate") (1 - d) / (1 - p) else (1 - d) * p / (1 - p)) else NULL
    b1 <- b0 <- NULL
    if (method != "ipw") {
        w0 <- switch(method, ra = 1, aipw = 1, ipwra = if (estimand == "ate") 1 / (1 - p) else p / (1 - p))
        b0 <- mf_wls(W, y, sw * (1 - d) * w0)
        if (estimand %in% c("ate", "ato")) b1 <- mf_wls(W, y, sw * d * (if (method == "ipwra") 1 / p else 1))
    }

    # ---- layout of theta and the moment function -----------------------------------------------------------------------------
    pos <- 0; blk <- function(k) { r <- pos + seq_len(k); pos <<- pos + k; r }
    ig <- if (ps) blk(kx) else integer(0)
    ib1 <- if (estimand %in% c("ate", "ato") && method != "ipw") blk(kw) else integer(0)
    ib0 <- if (method != "ipw") blk(kw) else integer(0)
    i1 <- blk(1); i0 <- blk(1)
    theta <- numeric(pos)
    if (ps) theta[ig] <- gamma
    if (length(ib1)) theta[ib1] <- b1
    if (length(ib0)) theta[ib0] <- b0

    psi_fun <- function(th) {
        g <- if (ps) th[ig] else NULL
        pp <- if (ps) mf_ps_prob(X, g, link) else NULL
        mu1 <- if (length(ib1)) drop(W %*% th[ib1]) else NULL
        mu0 <- if (length(ib0)) drop(W %*% th[ib0]) else NULL
        m1 <- th[i1]; m0 <- th[i0]
        cols <- list()
        if (ps) cols$ps <- mf_ps_score(d, X, g, link, estimand)
        w1 <- if (ps) d / pp else NULL
        w0 <- if (ps) (if (estimand == "ate") (1 - d) / (1 - pp) else (1 - d) * pp / (1 - pp)) else NULL
        sc1 <- d * (if (method == "ipwra") w1 else 1)
        sc0 <- if (estimand %in% c("ate", "ato")) (1 - d) * (if (method == "ipwra") 1 / (1 - pp) else 1)
               else (1 - d) * (if (method == "ipwra") pp / (1 - pp) else 1)
        if (method %in% c("ra", "aipw", "ipwra")) {
            if (length(ib1)) cols$b1 <- W * (sc1 * (y - mu1))
            cols$b0 <- W * (sc0 * (y - mu0))
        }
        if (estimand == "ato") {                                 # overlap weights: h = p (1 - p); treated weight 1 - p, control weight p
            h <- pp * (1 - pp)
            e1 <- switch(method, ra = h * (mu1 - m1), ipw = d * (1 - pp) * (y - m1), aipw = h * mu1 + d * (1 - pp) * (y - mu1) - h * m1)
            e0 <- switch(method, ra = h * (mu0 - m0), ipw = (1 - d) * pp * (y - m0), aipw = h * mu0 + (1 - d) * pp * (y - mu0) - h * m0)
        } else if (estimand == "ate") {
            e1 <- switch(method, ra = mu1 - m1, ipw = w1 * (y - m1), aipw = mu1 + w1 * (y - mu1) - m1, ipwra = mu1 - m1)
            e0 <- switch(method, ra = mu0 - m0, ipw = w0 * (y - m0), aipw = mu0 + w0 * (y - mu0) - m0, ipwra = mu0 - m0)
        } else {
            e1 <- d * (y - m1)
            e0 <- switch(method, ra = d * (mu0 - m0), ipw = w0 * (y - m0), aipw = d * mu0 + w0 * (y - mu0) - d * m0, ipwra = d * (mu0 - m0))
        }
        cols$m1 <- e1; cols$m0 <- e0
        do.call(cbind, cols) * sw
    }
    # the potential-outcome means are linear in the moments: solve them from the moments at 0 and 1
    for (i in c(i1, i0)) {
        t0 <- theta; t1 <- theta; t0[i] <- 0; t1[i] <- 1
        s0 <- colSums(psi_fun(t0)); s1 <- colSums(psi_fun(t1))
        col <- ncol(psi_fun(theta)) - (if (i == i1) 1 else 0)
        theta[i] <- -s0[col] / (s1[col] - s0[col])
    }
    psi <- psi_fun(theta)
    chk <- max(abs(colMeans(psi)))
    # ---- M-estimation covariance ---------------------------------------------------------------------------------------------------
    q <- length(theta)
    G <- matrix(0, q, q)
    for (j in seq_len(q)) {
        step <- 1e-6 * (1 + abs(theta[j])); t1 <- theta; t2 <- theta; t1[j] <- t1[j] + step; t2[j] <- t2[j] - step
        G[, j] <- (colMeans(psi_fun(t1)) - colMeans(psi_fun(t2))) / (2 * step)
    }
    Gi <- tryCatch(solve(G), error = function(e) stop("The variance matrix could not be computed (singular Jacobian); check collinearity and overlap.", call. = FALSE))
    Bm <- if (is.null(cl)) crossprod(psi) / n else crossprod(rowsum(psi, cl)) / n
    V <- Gi %*% Bm %*% t(Gi) / n
    cvec <- numeric(q); cvec[i1] <- 1; cvec[i0] <- -1
    eff <- theta[i1] - theta[i0]
    se <- sqrt(drop(t(cvec) %*% V %*% cvec))
    pom <- c(theta[i1], theta[i0]); pse <- sqrt(c(V[i1, i1], V[i0, i0]))
    eff <- eff * ysc; se <- se * ysc; pom <- pom * ysc + ymu; pse <- pse * ysc; chk <- chk * ysc
    if (flip) { eff <- -eff; pom <- rev(pom); pse <- rev(pse) }
    list(effect = eff, se = se, pom1 = pom[1], pom0 = pom[2], se_pom1 = pse[1], se_pom0 = pse[2], moment_check = chk,
         ps = p, V = V, theta = theta, method = method, estimand = if (flip) "atc" else estimand,
         n = n, flipped = flip)
}

# Inverse-probability weights for diagnostics: weights of the treated and of the controls for the chosen estimand.
mf_weights <- function(d, p, estimand = c("ate", "att", "atc", "ato")) {
    estimand <- match.arg(estimand)
    switch(estimand,
           ato = ifelse(d == 1, 1 - p, p),
           ate = ifelse(d == 1, 1 / p, 1 / (1 - p)),
           att = ifelse(d == 1, 1, p / (1 - p)),
           atc = ifelse(d == 1, (1 - p) / p, 1))
}

mf_weight_summary <- function(w, d) {
    do.call(rbind, lapply(c(1, 0), function(g) {
        x <- w[d == g]
        data.frame(group = if (g == 1) "Treated" else "Control", n = length(x), min = min(x), mean = mean(x), max = max(x),
                   ess = mf_ess(x), share_max = max(x) / sum(x), stringsAsFactors = FALSE)
    }))
}

# Unadjusted difference in means with a robust standard error (the naive comparison).
mf_naive <- function(y, d, sw = rep(1, length(y))) {
    m1 <- sum(sw * d * y) / sum(sw * d); m0 <- sum(sw * (1 - d) * y) / sum(sw * (1 - d))
    e1 <- y - m1; e0 <- y - m0
    v1 <- sum((sw * d * e1)^2) / sum(sw * d)^2; v0 <- sum((sw * (1 - d) * e0)^2) / sum(sw * (1 - d))^2
    list(effect = m1 - m0, se = sqrt(v1 + v0), pom1 = m1, pom0 = m0)
}
