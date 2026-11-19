# Propensity score engine: sample preparation, logit / probit propensity score, overlap, trimming.
# Base R only; nothing here depends on jamovi, so it can be tested from a plain R session.

mf_num <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)

# Design matrix with an intercept: numeric variables as they are, factors and text as treatment-coded dummies.
mf_design <- function(data, vars, intercept = TRUE) {
    if (!length(vars)) return(if (intercept) matrix(1, nrow(data), 1, dimnames = list(NULL, "(Intercept)")) else NULL)
    df <- data[, vars, drop = FALSE]
    for (v in vars) {
        if (is.character(df[[v]]) || is.logical(df[[v]])) df[[v]] <- factor(df[[v]])
        else if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]])
    }
    mm <- stats::model.matrix(~ ., data = df)
    if (!intercept) mm <- mm[, -1, drop = FALSE]
    mm
}

# Build the estimation sample. A two-level factor or text treatment becomes 0/1 with the second level as treated.
# Rows with a missing value in any used variable are dropped.
mf_prepare <- function(data, treat, covs, outcome = NULL, extra = NULL, weights = NULL, cluster = NULL, exact = NULL, aux = NULL) {
    if (is.null(treat)) stop("Select a treatment variable.", call. = FALSE)
    if (!length(covs)) stop("Select at least one covariate for the propensity score.", call. = FALSE)
    if (any(duplicated(c(outcome, treat, covs, extra, exact))))
        stop("A variable appears in more than one role (outcome, treatment, covariates).", call. = FALSE)
    used <- unique(c(outcome, treat, covs, extra, exact, weights, cluster, aux))
    keep <- stats::complete.cases(data[, used, drop = FALSE])
    tv <- data[[treat]]
    labels <- NULL
    if (is.factor(tv) || is.character(tv) || is.logical(tv)) {
        lv <- levels(factor(tv[keep]))
        if (length(lv) != 2) stop(sprintf("The treatment %s must have exactly two values.", treat), call. = FALSE)
        d_all <- as.numeric(factor(tv, levels = lv)) - 1
        labels <- lv
    } else d_all <- as.numeric(tv)
    if (!is.null(outcome)) keep <- keep & is.finite(mf_num(data[[outcome]]))
    keep <- keep & is.finite(d_all)
    n_drop <- sum(!keep)
    dd <- data[keep, , drop = FALSE]
    d <- d_all[keep]
    if (!all(d %in% c(0, 1)) || length(unique(d)) < 2)
        stop(sprintf("The treatment %s must be coded 0/1 with both groups present.", treat), call. = FALSE)
    if (min(table(d)) < 5) stop("Each treatment group needs at least 5 observations.", call. = FALSE)
    X <- mf_design(dd, covs)
    W <- if (length(extra)) cbind(X, mf_design(dd, extra, intercept = FALSE)) else X
    sw <- if (!is.null(weights)) mf_num(dd[[weights]]) else rep(1, nrow(dd))
    if (any(sw <= 0)) stop("Weights must be positive.", call. = FALSE)
    list(y = if (!is.null(outcome)) mf_num(dd[[outcome]]) else NULL, d = d, X = X, W = W, sw = sw,
         cl = if (!is.null(cluster)) as.character(dd[[cluster]]) else NULL,
         aux = if (!is.null(aux)) mf_num(dd[[aux]]) else NULL,
         ex = if (length(exact)) dd[, exact, drop = FALSE] else NULL,
         n_drop = n_drop, labels = labels, covs = covs, data = dd, keep = keep)
}

# Propensity score by maximum likelihood (logit or probit), with an informative error when it cannot be estimated.
# link "cbps": covariate balancing propensity score (Imai and Ratkovic, 2014), exactly identified: the logit model whose coefficients solve the
# balancing conditions  sum sw x (d / p - (1 - d) / (1 - p)) = 0  for the ATE and  sum sw x (d - (1 - d) p / (1 - p)) = 0  for the ATT (the ATC is the ATT with
# the groups swapped), so that the weighted covariate means of the two groups equal the target means exactly.
mf_ps_fit <- function(d, X, sw = rep(1, length(d)), link = c("logit", "probit", "cbps"), estimand = "ate") {
    link <- match.arg(link)
    q <- qr(X)
    if (q$rank < ncol(X)) {
        red <- colnames(X)[q$pivot[(q$rank + 1):ncol(X)]]
        stop(sprintf("The covariates are perfectly collinear: %s can be written from the other covariates. Remove a redundant covariate; a numeric variable set to Nominal (for example years of education next to an indicator derived from it) is a common cause.",
                     paste(sprintf("'%s'", red), collapse = ", ")), call. = FALSE)
    }
    fit <- suppressWarnings(stats::glm.fit(X, d, weights = sw, family = stats::binomial(if (link == "cbps") "logit" else link), control = list(maxit = 100, epsilon = 1e-12)))
    if (link == "cbps" && !anyNA(fit$coefficients)) {
        g <- fit$coefficients
        for (it in 1:100) {
            p <- stats::plogis(drop(X %*% g)); p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
            if (estimand == "att") { r <- d - (1 - d) * p / (1 - p); dr <- -(1 - d) * p / (1 - p) }
            else { r <- d / p - (1 - d) / (1 - p); dr <- -d * (1 - p) / p - (1 - d) * p / (1 - p) }
            gv <- colSums(sw * X * r)
            Jm <- crossprod(X, sw * dr * X)
            step <- tryCatch(solve(Jm, gv), error = function(e) rep(NA_real_, length(g)))
            if (anyNA(step)) break
            lam <- 1
            repeat {                                           # halve the step until the balance error does not grow
                gn <- g - lam * step; pn <- pmin(pmax(stats::plogis(drop(X %*% gn)), 1e-10), 1 - 1e-10)
                rn <- if (estimand == "att") d - (1 - d) * pn / (1 - pn) else d / pn - (1 - d) / (1 - pn)
                if (sqrt(sum(colSums(sw * X * rn)^2)) <= sqrt(sum(gv^2)) * (1 - 1e-4 * lam) || lam < 1e-6) break
                lam <- lam / 2
            }
            g <- gn
            if (max(abs(colMeans(sw * X * rn))) < 1e-10) break
        }
        fit$coefficients <- g; fit$fitted.values <- stats::plogis(drop(X %*% g)); fit$converged <- it < 100
    }
    p <- as.numeric(fit$fitted.values)
    sep <- !fit$converged || any(abs(fit$coefficients) > 25, na.rm = TRUE) || any(p < 1e-8 | p > 1 - 1e-8)
    if (anyNA(fit$coefficients)) stop("The propensity score model has a perfectly collinear covariate.", call. = FALSE)
    if (sep) stop("The propensity score model does not converge: a covariate (or a combination) predicts treatment almost perfectly. ",
                  "Remove it, merge sparse categories, or use fewer covariates.", call. = FALSE)
    list(coef = fit$coefficients, p = p, link = link, converged = fit$converged, X = X)
}

# P(D = 1 | X) for any estimand: for the ATC with the cbps link the balancing conditions are those of the swapped groups.
mf_ps_fit_for <- function(d, X, sw, link, estimand) {
    if (identical(link, "cbps") && estimand == "atc") { f <- mf_ps_fit(1 - d, X, sw, link, "att"); f$p <- 1 - f$p; return(f) }
    mf_ps_fit(d, X, sw, link, if (identical(link, "cbps") && estimand == "ato") "ate" else estimand)
}

# Score of the log-likelihood of one observation for the propensity score parameters (rows).
mf_ps_score <- function(d, X, gamma, link, estimand = "ate") {
    eta <- drop(X %*% gamma)
    if (link == "cbps") { p <- stats::plogis(eta); return(X * (if (estimand == "att") d - (1 - d) * p / (1 - p) else d / p - (1 - d) / (1 - p))) }
    if (link == "logit") { p <- stats::plogis(eta); X * (d - p) }
    else { p <- stats::pnorm(eta); X * (stats::dnorm(eta) * (d - p) / (p * (1 - p))) }
}
mf_ps_prob <- function(X, gamma, link) {
    eta <- drop(X %*% gamma)
    if (link %in% c("logit", "cbps")) stats::plogis(eta) else stats::pnorm(eta)
}

# Standard errors of the propensity score coefficients: the inverse of the observed information (robust sandwich).
mf_ps_vcov <- function(d, X, sw, gamma, link, estimand = "ate") {
    S <- mf_ps_score(d, X, gamma, link, estimand) * sw
    k <- ncol(X); n <- nrow(X)
    h <- 1e-6
    J <- matrix(0, k, k)
    for (j in seq_len(k)) {
        g1 <- gamma; g2 <- gamma; step <- h * (1 + abs(gamma[j])); g1[j] <- g1[j] + step; g2[j] <- g2[j] - step
        J[, j] <- (colMeans(mf_ps_score(d, X, g1, link, estimand) * sw) - colMeans(mf_ps_score(d, X, g2, link, estimand) * sw)) / (2 * step)
    }
    Ji <- solve(J)
    Ji %*% (crossprod(S) / n) %*% t(Ji) / n
}

# Area under the ROC curve of the propensity score (Mann-Whitney statistic).
mf_auc <- function(p, d) {
    r <- rank(p); n1 <- sum(d == 1); n0 <- sum(d == 0)
    (sum(r[d == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Common support and trimming. rule: "none", "minmax" (observations outside the overlap of the two groups' ranges),
# "crump" (propensity score outside [0.1, 0.9]) or "custom" with lo and hi. Returns a logical vector of retained observations.
mf_trim <- function(p, d, rule = c("none", "minmax", "crump", "custom"), lo = 0.05, hi = 0.95) {
    rule <- match.arg(rule)
    switch(rule,
           none = rep(TRUE, length(p)),
           minmax = { a <- max(min(p[d == 1]), min(p[d == 0])); b <- min(max(p[d == 1]), max(p[d == 0])); p >= a & p <= b },
           crump = p >= 0.1 & p <= 0.9,
           custom = p >= lo & p <= hi)
}

mf_ess <- function(w) sum(w)^2 / sum(w^2)
