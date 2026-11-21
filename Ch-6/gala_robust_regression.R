# Robust regression (Huber, LAD, least trimmed squares) and the bootstrap, using the gala data
# Packages: MASS (comes with R) and quantreg (install once with install.packages("quantreg"))

library(faraway)
library(MASS)       # rlm() for Huber M-estimation, ltsreg() for least trimmed squares
data(gala)

# Least squares, for reference
gl <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala)
print(summary(gl))

# ---- Huber M-estimation ----
# rlm(): robust linear model; Huber's method is the default.
# It down-weights cases with large residuals, refitting by iteratively reweighted least squares
gr <- rlm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala)
print(summary(gr))

# ---- Least absolute deviations (LAD / L1) ----
# rq(): quantile regression; the default tau = 0.5 (the median) is LAD regression
if (requireNamespace("quantreg", quietly = TRUE)) {
  library(quantreg)
  gq <- rq(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, data = gala)
  print(summary(gq))
} else {
  cat("Install the quantreg package for the LAD fit: install.packages('quantreg')\n")
}

# ---- Least trimmed squares (LTS) ----
# ltsreg(): minimises the sum of the q smallest squared residuals, ignoring the rest.
# The default search is random, so two runs can give different answers
set.seed(1)
g <- ltsreg(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala)
print(coef(g))
g <- ltsreg(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala)
print(coef(g))
# nsamp = "exact": search every possible subset - slower, but always the same answer
g <- ltsreg(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala, nsamp = "exact")
print(coef(g))


# ---- Bootstrap for the LTS coefficients (no standard errors are available otherwise) ----
# sample(10, rep = TRUE): 10 draws from 1..10 with replacement
print(sample(10, rep = TRUE))
# A random resample of the LTS residuals
print(head(residuals(g)[sample(30, rep = TRUE)]))

# The fixed X matrix (without the intercept, which ltsreg() adds itself)
x <- model.matrix(~ Area + Elevation + Nearest + Scruz + Adjacent, gala)[, -1]
nboot <- 1000
# matrix(0, 1000, 6): a 1000 x 6 matrix of zeros to hold the bootstrap coefficients
bcoef <- matrix(0, nboot, 6)
set.seed(123)
for (i in 1:nboot) {
  # New response = LTS fitted values + resampled LTS residuals
  newy <- predict(g) + residuals(g)[sample(30, rep = TRUE)]
  # nsamp = "best" is much faster than "exact" (needed 1000 times here)
  brg <- ltsreg(x, newy, nsamp = "best")
  bcoef[i, ] <- coef(brg)
}
# quantile(): the 2.5% and 97.5% points give a 95% bootstrap CI for the Area coefficient
print(quantile(bcoef[, 2], c(0.025, 0.975)))
plot(density(bcoef[, 2]), xlab = "Coefficient of Area", main = "Bootstrap distribution")
abline(v = quantile(bcoef[, 2], c(0.025, 0.975)))


# ---- Least squares without the very influential island Isabela ----
print(round(sort(cooks.distance(gl), decreasing = TRUE)[1:3], 2))
gli <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala,
          subset = (row.names(gala) != "Isabela"))
print(summary(gli))

# What it shows (numbers from the book; the random parts will differ):
# - Huber (rlm): the coefficients move a little (Elevation 0.32 -> 0.25,
#   Adjacent -0.075 -> -0.055) but the same variables matter. Only approximate
#   inference is available (t-values, no p-values or R^2).
# - LAD (rq): again a modest change; the rank-based CIs now suggest Adjacent
#   may not be significant.
# - LTS: very different. With the exact search the book got Area 1.54,
#   Elevation 0.024, Adjacent -0.198: Area is now large and Elevation small.
# - The bootstrap 95% CI for the Area coefficient was (1.49, 1.62) in the
#   book: zero is far outside, so under LTS there is a clear Area effect,
#   unlike the least squares fit. The bootstrap distribution is peaked with
#   longish tails.
# - Isabela has by far the largest Cook's distance. Without it, least squares
#   agrees with LTS: Area (p = 8e-05) and Adjacent (p = 2e-05) are both highly
#   significant, R^2 = 0.87. Diagnostics + least squares and robust methods
#   reach the same place; the Huber M-estimate did not spot Isabela.
