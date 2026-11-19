# Entropy balancing (Hainmueller 2012): weights that exactly match the moments of the covariates, closest to the base weights in
# Kullback-Leibler distance. Solved through the convex dual; the effect and its standard error come from stacked moment
# conditions (the multipliers and the target means are estimated), as in an M-estimator. Base R only.

# Balancing features: the covariate columns, plus squares (second moment) and cubes (third) of the non-binary ones.
# Features are centred and scaled on the whole sample, and linearly dependent columns are dropped.
mf_ebal_features <- function(Xc, moments = 1) {
    nb <- apply(Xc, 2, function(z) !all(z %in% c(0, 1)))
    F <- Xc
    if (moments >= 2 && any(nb)) { A <- Xc[, nb, drop = FALSE]^2; colnames(A) <- paste0(colnames(A), "^2"); F <- cbind(F, A) }
    if (moments >= 3 && any(nb)) { A <- Xc[, nb, drop = FALSE]^3; colnames(A) <- paste0(colnames(A), "^3"); F <- cbind(F, A) }
    s <- apply(F, 2, stats::sd)
    if (any(!is.finite(s) | s <= 0)) F <- F[, is.finite(s) & s > 0, drop = FALSE]
    F <- scale(F)
    attr(F, "scaled:center") <- NULL; attr(F, "scaled:scale") <- NULL
    qrF <- qr(F, tol = 1e-9)
    if (qrF$rank < ncol(F)) F <- F[, sort(qrF$pivot[seq_len(qrF$rank)]), drop = FALSE]
    F
}

# Dual problem: minimise log sum(q exp(C lambda)). Returns the multipliers; the weights are q exp(C lambda).
mf_ebal_dual <- function(C, q, maxit = 200, tol = 1e-11) {
    k <- ncol(C); lam <- numeric(k)
    obj <- function(l) { z <- drop(C %*% l); mz <- max(z); mz + log(sum(q * exp(z - mz))) }
    f <- obj(lam)
    for (it in seq_len(maxit)) {
        z <- drop(C %*% lam); a <- q * exp(z - max(z)); a <- a / sum(a)
        g <- drop(crossprod(C, a))
        if (max(abs(g)) < tol) return(list(lambda = lam, converged = TRUE, iter = it))
        H <- crossprod(C, a * C) - tcrossprod(g)
        step <- tryCatch(-solve(H, g), error = function(e) -g)
        t <- 1
        repeat {
            fn <- obj(lam + t * step)
            if (is.finite(fn) && fn <= f + 1e-4 * t * sum(g * step)) break
            t <- t / 2
            if (t < 1e-12) break
        }
        lam <- lam + t * step; f <- obj(lam)
        if (max(abs(lam)) > 200) break
    }
    list(lambda = lam, converged = FALSE, iter = maxit)
}

mf_ebal_fail <- function() stop("Entropy balancing cannot match the requested moments: the target values lie outside the range of the group being reweighted ",
                                "(no positive weights can reach them). Use fewer covariates or moments, or trim the sample first.", call. = FALSE)

# Weights and target means (no outcome needed). `estimand`: the group(s) being reweighted follow it.
mf_ebal_weights <- function(d, Xc, sw = rep(1, length(d)), estimand = c("ate", "att", "atc"), moments = 1) {
    estimand <- match.arg(estimand)
    F <- mf_ebal_features(Xc, moments)
    tg <- switch(estimand, ate = rep(TRUE, length(d)), att = d == 1, atc = d == 0)
    m <- colSums(F[tg, , drop = FALSE] * sw[tg]) / sum(sw[tg])
    C <- sweep(F, 2, m)
    w <- sw; lam <- list()
    for (g in c(1, 0)) {
        if ((estimand == "att" && g == 1) || (estimand == "atc" && g == 0)) next
        i <- d == g
        r <- mf_ebal_dual(C[i, , drop = FALSE], sw[i])
        if (!r$converged) mf_ebal_fail()
        e <- sw[i] * exp(drop(C[i, , drop = FALSE] %*% r$lambda) - max(drop(C[i, , drop = FALSE] %*% r$lambda)))
        w[i] <- e * sum(sw[i]) / sum(e)
        lam[[as.character(g)]] <- r$lambda
    }
    imb <- max(vapply(c(1, 0), function(g) {
        i <- d == g; if (all(w[i] == sw[i])) return(0)
        max(abs(colSums(C[i, , drop = FALSE] * w[i]) / sum(w[i])))
    }, numeric(1)))
    list(w = w, F = F, target = m, lambda = lam, imbalance = imb)
}

# Effect after entropy balancing with M-estimation standard errors.
mf_ebal <- function(y, d, Xc, sw = rep(1, length(y)), cl = NULL, estimand = c("ate", "att", "atc"), moments = 1) {
    estimand <- match.arg(estimand)
    ew <- mf_ebal_weights(d, Xc, sw, estimand, moments)
    F <- ew$F; k <- ncol(F); n <- length(y)
    ysc <- stats::sd(y); if (!is.finite(ysc) || ysc <= 0) ysc <- 1
    ymu <- mean(y); ys <- (y - ymu) / ysc
    rw1 <- estimand != "att"; rw0 <- estimand != "atc"
    tg <- switch(estimand, ate = rep(1, n), att = d, atc = 1 - d)
    pos <- 0; blk <- function(m_) { r <- pos + seq_len(m_); pos <<- pos + m_; r }
    im <- blk(k); il1 <- if (rw1) blk(k) else integer(0); il0 <- if (rw0) blk(k) else integer(0); i1 <- blk(1); i0 <- blk(1)
    psi_fun <- function(th) {
        m <- th[im]; C <- sweep(F, 2, m)
        e1 <- if (rw1) exp(drop(C %*% th[il1])) else 1
        e0 <- if (rw0) exp(drop(C %*% th[il0])) else 1
        cols <- list(tg * (F - matrix(m, n, k, byrow = TRUE)))
        if (rw1) cols[[length(cols) + 1]] <- d * e1 * C
        if (rw0) cols[[length(cols) + 1]] <- (1 - d) * e0 * C
        cols[[length(cols) + 1]] <- d * e1 * (ys - th[i1])
        cols[[length(cols) + 1]] <- (1 - d) * e0 * (ys - th[i0])
        do.call(cbind, cols) * sw
    }
    theta <- numeric(pos)
    theta[im] <- ew$target
    if (rw1) theta[il1] <- ew$lambda[["1"]]
    if (rw0) theta[il0] <- ew$lambda[["0"]]
    # the weighted outcome means are the solutions of the last two (linear) moments
    C <- sweep(F, 2, ew$target)
    e1 <- if (rw1) exp(drop(C %*% theta[il1])) else 1
    e0 <- if (rw0) exp(drop(C %*% theta[il0])) else 1
    theta[i1] <- sum(sw * d * e1 * ys) / sum(sw * d * e1)
    theta[i0] <- sum(sw * (1 - d) * e0 * ys) / sum(sw * (1 - d) * e0)
    psi <- psi_fun(theta)
    chk <- max(abs(colMeans(psi)))
    q <- length(theta); G <- matrix(0, q, q)
    for (j in seq_len(q)) {
        step <- 1e-6 * (1 + abs(theta[j])); t1 <- theta; t2 <- theta; t1[j] <- t1[j] + step; t2[j] <- t2[j] - step
        G[, j] <- (colMeans(psi_fun(t1)) - colMeans(psi_fun(t2))) / (2 * step)
    }
    Gi <- tryCatch(solve(G), error = function(e) stop("The variance matrix could not be computed (singular Jacobian); check collinearity and overlap.", call. = FALSE))
    Bm <- if (is.null(cl)) crossprod(psi) / n else crossprod(rowsum(psi, cl)) / n
    V <- Gi %*% Bm %*% t(Gi) / n
    cvec <- numeric(q); cvec[i1] <- 1; cvec[i0] <- -1
    eff <- theta[i1] - theta[i0]; se <- sqrt(drop(t(cvec) %*% V %*% cvec))
    list(effect = eff * ysc, se = se * ysc, pom1 = theta[i1] * ysc + ymu, pom0 = theta[i0] * ysc + ymu,
         se_pom1 = sqrt(V[i1, i1]) * ysc, se_pom0 = sqrt(V[i0, i0]) * ysc, moment_check = chk * ysc, weights = ew$w,
         imbalance = ew$imbalance, n_features = k, n = n)
}
