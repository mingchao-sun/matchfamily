# Builds the example data. Run from the module root: Rscript dev/make_data.R
# lalonde.csv comes from the CRAN package MatchIt (GPL-2); the simulations are generated here.
suppressMessages(library(MatchIt))
data("lalonde", package = "MatchIt")
ll <- lalonde
ll$race <- as.character(ll$race)
write.csv(ll, "data/lalonde.csv", row.names = FALSE)
cat("lalonde.csv:", nrow(ll), "rows\n")

# Confounded treatment with a known effect: the effect rises with x1, so ATE, ATT and ATC differ. Good overlap.
set.seed(2026)
n <- 3000
x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rbinom(n, 1, 0.4)
p <- plogis(-0.3 + 0.8 * x1 - 0.5 * x2 + 0.6 * x3)
d <- rbinom(n, 1, p)
tau <- 1 + 0.5 * x1
y0 <- 2 + 1.2 * x1 + 0.8 * x2 + 0.5 * x3 + rnorm(n)
y <- y0 + d * tau
sim <- data.frame(y = round(y, 4), treat = d, x1 = round(x1, 4), x2 = round(x2, 4), x3 = x3)
write.csv(sim, "data/sim_confounded.csv", row.names = FALSE)
cat(sprintf("sim_confounded.csv: ATE %.4f, ATT %.4f, ATC %.4f (sample values), treated %d of %d\n",
            mean(tau), mean(tau[d == 1]), mean(tau[d == 0]), sum(d), n))

# Poor overlap: strong selection on x1, so some treated units have almost no comparable controls.
set.seed(2027)
n <- 3000
x1 <- rnorm(n); x2 <- rnorm(n)
p <- plogis(-0.5 + 2.4 * x1 + 0.5 * x2)
d <- rbinom(n, 1, p)
y <- 1 + 1.0 * x1 + 0.6 * x2 + 1.5 * d + rnorm(n)
write.csv(data.frame(y = round(y, 4), treat = d, x1 = round(x1, 4), x2 = round(x2, 4)), "data/sim_overlap.csv", row.names = FALSE)
cat(sprintf("sim_overlap.csv: treated %d of %d, true effect 1.5\n", sum(d), n))

# Hidden confounder u (not in the data): the naive and adjusted estimates are biased, to try the sensitivity analysis.
set.seed(2028)
n <- 2000
u <- rnorm(n); x1 <- rnorm(n) + 0.3 * u; x2 <- rbinom(n, 1, 0.5)
d <- rbinom(n, 1, plogis(0.6 * x1 + 0.7 * u - 0.2))
y <- 1 + 0.8 * x1 + 0.4 * x2 + 0.9 * u + 1.0 * d + rnorm(n)
write.csv(data.frame(y = round(y, 4), treat = d, x1 = round(x1, 4), x2 = x2), "data/sim_hidden.csv", row.names = FALSE)
cat(sprintf("sim_hidden.csv: treated %d of %d, true effect 1.0\n", sum(d), n))
# jamovi guesses integer CSV columns (age, educ) as Nominal; the .omv file fixes the measure types for the Data Library
if (requireNamespace("jmvReadWrite", quietly = TRUE)) {
    lo <- read.csv("data/lalonde.csv", stringsAsFactors = FALSE); lo$race <- factor(lo$race)
    jmvReadWrite::write_omv(lo, "data/lalonde.omv")
}
