# Generalized least squares (GLS) with AR(1) errors, using Longley's data (US economy, 1947-1962)
# Response: Employed (millions); predictors: GNP and Population. longley comes with R itself.

data(longley)
g <- lm(Employed ~ GNP + Population, longley)
print(summary(g, correlation = TRUE))

# The predictors are strongly positively correlated; their coefficients strongly negatively
print(cor(longley$GNP, longley$Population))

# Successive errors in yearly data may be correlated. Estimate the lag-1 correlation
# from the residuals: residuals(g)[-1] are years 2..16, residuals(g)[-16] are years 1..15
rho <- cor(residuals(g)[-1], residuals(g)[-16])
print(rho)


# ---- GLS by hand, treating rho as known ----
# Under AR(1) errors the correlation between years i and j is rho^|i - j|
x <- model.matrix(g)
Sigma <- diag(16)                                   # diag(16): 16 x 16 identity matrix
# row(), col(): row and column number of each cell, so abs(row - col) = |i - j|
Sigma <- rho^abs(row(Sigma) - col(Sigma))
Sigi <- solve(Sigma)                                # Sigma^-1
xtxi <- solve(t(x) %*% Sigi %*% x)
# beta-hat_GLS = (X' Sigma^-1 X)^-1 X' Sigma^-1 y
beta <- solve(t(x) %*% Sigi %*% x, t(x) %*% Sigi %*% longley$Employed)
print(beta)
res <- longley$Employed - x %*% beta
sig <- drop(sqrt((t(res) %*% Sigi %*% res) / df.residual(g)))
print(sig)
print(sqrt(diag(xtxi)) * sig)                       # GLS standard errors

# Equivalent: transform to an ordinary regression. chol(Sigma) gives S' with S S' = Sigma;
# multiplying y and X by S^-1 makes the errors independent
sm <- chol(Sigma)
smi <- solve(t(sm))
sx <- smi %*% x
sy <- smi %*% longley$Employed
# sy ~ sx - 1: no extra intercept, since sx already contains the transformed intercept column
print(summary(lm(sy ~ sx - 1)))

# rho should now be re-estimated from the GLS residuals, and the process repeated
print(cor(res[-1], res[-16]))


# ---- The convenient way: gls() from the nlme package ----
# library(nlme): gls() estimates beta and the correlation together (by REML by default)
library(nlme)
# correlation = corAR1(form = ~ Year): AR(1) errors ordered by Year
g2 <- gls(Employed ~ GNP + Population, correlation = corAR1(form = ~ Year), data = longley)
print(summary(g2))
# intervals(): approximate 95% CIs for the coefficients, phi (the correlation) and sigma
print(intervals(g2))

# What it shows (numbers from the book):
# - OLS: GNP 0.0632, Population -0.410, residual SE 0.546, R^2 = 0.979.
#   cor(GNP, Population) = 0.991 while the correlation of their estimates is
#   -0.99: positively correlated predictors give negatively correlated estimates.
# - Lag-1 residual correlation = 0.31. With rho = 0.31, GLS gives
#   94.90 + 0.0674 GNP - 0.474 Population, with sigma 0.542. The transformed
#   regression gives the same numbers.
# - Re-estimating rho from the GLS residuals gives 0.356, so the by-hand
#   method would have to be iterated. gls() does it properly.
# - gls(): phi = 0.644, coefficients 101.86, 0.0721, -0.549. But the 95% CI for
#   phi is (-0.44, 0.96), which includes 0: with only 16 years there is no
#   real evidence of serial correlation.
