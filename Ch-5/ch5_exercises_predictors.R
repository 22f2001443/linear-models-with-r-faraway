# =============================================================================
# Linear Models with R (Faraway), Chapter 5 exercises 1-5
# Measurement error (SIMEX) on faithful, and collinearity diagnostics
# (condition numbers, correlations, VIFs) on divusa, longley and prostate.
# =============================================================================

library(faraway)

# -----------------------------------------------------------------------------
# collinearity(g): condition numbers, predictor correlations and VIFs for model g
# -----------------------------------------------------------------------------
collinearity <- function(g) {
  x <- model.matrix(g)[, -1]            # X without the intercept column
  e <- eigen(t(x) %*% x)
  cat("Eigenvalues of X'X:\n");   print(signif(e$values, 4))
  cat("Condition numbers sqrt(l1 / lj):\n"); print(round(sqrt(e$values[1] / e$values), 1))
  # The same on standardised predictors (cor matrix), so units do not drive the answer
  es <- eigen(cor(x))
  cat("Condition numbers, standardised predictors:\n")
  print(round(sqrt(es$values[1] / es$values), 1))
  cat("Correlations between predictors:\n"); print(round(cor(x), 3))
  cat("VIFs:\n"); print(round(vif(x), 2))
  invisible(x)
}


# =============================================================================
# Exercise 1: faithful - Old Faithful geyser, 272 eruptions (comes with R)
#   eruptions = eruption duration (min), waiting = waiting time to this eruption (min)
#   Assume waiting has measurement error with SD 30 seconds = 0.5 min, variance 0.25
# =============================================================================
data(faithful)
g <- lm(eruptions ~ waiting, data = faithful)
print(summary(g))

set.seed(123)
base_var <- 0.5^2                        # assumed error variance already in waiting
extra <- c(0.05, 0.10, 0.15, 0.20, 0.25) # extra error variance we add on top
nrep <- 1000
n <- nrow(faithful)
mean_slopes <- numeric(length(extra))
for (j in seq_along(extra)) {
  s <- numeric(nrep)
  for (i in 1:nrep) {
    s[i] <- coef(lm(eruptions ~ I(waiting + sqrt(extra[j]) * rnorm(n)), data = faithful))[2]
  }
  mean_slopes[j] <- mean(s)
}
betas <- c(coef(g)[2], mean_slopes)
variances <- c(0, extra) + base_var
print(cbind(variances, betas))

gv <- lm(betas ~ variances)
print(coef(gv))
plot(variances, betas, xlim = c(0, max(variances)),
     ylim = range(betas, coef(gv)[1]), xlab = "error variance", ylab = "mean slope",
     main = "SIMEX for faithful")
abline(gv, lty = 2)
points(0, coef(gv)[1], pch = 3, cex = 1.5)
# SIMEX estimate vs least squares estimate
print(c(least_squares = unname(coef(g)[2]), simex = unname(coef(gv)[1])))
# Why the correction is tiny: the attenuation factor is var(x) / (var(x) + error variance)
print(var(faithful$waiting) / (var(faithful$waiting) + base_var))

# What stands out (faithful):
# [numbers to be filled in from your run_all.R output]
# - The SIMEX slope is only very slightly steeper than least squares.
#   waiting has variance of about 185 min^2, so an error variance of 0.25
#   attenuates the slope by a factor of about 0.999: half a minute of
#   measurement error is negligible next to the spread of the waiting times.


# =============================================================================
# Exercise 2: SIMEX applied to error in the RESPONSE instead of the predictor
#   Demonstration on cars: add extra noise to dist and watch the mean slope
# =============================================================================
data(cars)
gc <- lm(dist ~ speed, cars)
extra_y <- c(10, 50, 100, 200, 400)        # extra variance added to dist
ms <- sapply(extra_y, function(v) {
  # sapply(): for each v, run 500 noisy fits and return the mean slope
  mean(replicate(500, coef(lm(I(dist + sqrt(v) * rnorm(50)) ~ speed, cars))[2]))
})
print(cbind(extra_variance = c(0, extra_y), mean_slope = c(coef(gc)[2], ms)))

# What stands out (response error):
# - The mean slope stays at about 3.93 however much noise is added. Error in
#   y is absorbed into the error term epsilon: least squares stays unbiased,
#   only its variance grows. The SIMEX line would be flat, so the
#   extrapolation would just return the least squares estimate - there is
#   no bias to correct.


# =============================================================================
# Exercise 3: divusa ("divorce data") - divorce on unemployed, femlab, marriage,
#   birth, military
# =============================================================================
data(divusa)
g3 <- lm(divorce ~ unemployed + femlab + marriage + birth + military, data = divusa)
print(summary(g3))

# (a) condition numbers, and (b) VIFs
collinearity(g3)

# (c) drop the predictors not significant at 5% and look again
coefs <- summary(g3)$coefficients
keep <- setdiff(rownames(coefs)[coefs[, 4] < 0.05], "(Intercept)")
print(keep)
g3b <- lm(reformulate(keep, response = "divorce"), data = divusa)
print(summary(g3b))
if (length(keep) >= 2) collinearity(g3b)

# What stands out (divusa collinearity):
# [numbers to be filled in from your run_all.R output]
# - (a) Raw condition numbers are driven partly by the different scales of the
#   predictors, which is why the standardised version is also shown.
# - (b) VIFs well below 10 would mean collinearity is NOT the reason some
#   predictors are non-significant; they would simply have little effect
#   once the others are in.
# - (c) Compare the VIFs and condition numbers before and after removal.


# =============================================================================
# Exercise 4: longley - Employed on the other 6 variables, 1947-1962 (comes with R)
# =============================================================================
data(longley)
g4 <- lm(Employed ~ ., data = longley)
print(summary(g4))
# (a) condition numbers, (b) correlations, (c) VIFs
collinearity(g4)

# What stands out (longley):
# [numbers to be filled in from your run_all.R output]
# - longley is the classic example of severe collinearity: GNP.deflator, GNP,
#   Population and Year all trend upward together over the 16 years, so
#   expect correlations close to 1, huge condition numbers and VIFs in the
#   hundreds or more.


# =============================================================================
# Exercise 5: prostate - lpsa on all other variables
# =============================================================================
data(prostate)
g5 <- lm(lpsa ~ ., data = prostate)
print(summary(g5))
# (a) condition numbers, (b) correlations, (c) VIFs
collinearity(g5)

# What stands out (prostate collinearity):
# [numbers to be filled in from your run_all.R output]
# - The raw condition numbers are large partly because age and pgg45 are on a
#   much bigger scale than the log variables; the standardised ones are the
#   fairer guide.
# - Expect moderate correlations (lcavol with lcp and svi; gleason with pgg45)
#   but VIFs that are all modest: collinearity is not a serious problem here.
