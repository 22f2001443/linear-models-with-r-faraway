# Least squares estimation, using the gala data (30 Galapagos islands)
# Columns: Species (number of plant species), Endemics, Area (km^2), Elevation (m),
# Nearest (km to nearest island), Scruz (km to Santa Cruz), Adjacent (km^2, area of adjacent island)

# library(): loads the faraway package, which contains the gala dataset
library(faraway)
# data(): loads the gala dataset into the workspace
data(gala)
# head(): first 6 rows, to see what the data look like
print(head(gala))

# lm(): fits a linear model by least squares
# The formula Species ~ Area + ... means "Species explained by these five predictors";
# data = gala tells lm() where to find the columns
mdl <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, data = gala)
# summary(): coefficients, standard errors, t-tests, residual SE, R^2 and the overall F-test
print(summary(mdl))


# ---- Doing the same calculation by hand with matrices ----

# model.matrix(): builds the X matrix (a column of 1s for the intercept, then one column per predictor)
x <- model.matrix(~ Area + Elevation + Nearest + Scruz + Adjacent, gala)
# The response vector y
y <- gala$Species

# t(): transpose.  %*%: matrix multiplication.  solve(A): the inverse of A
xtxi <- solve(t(x) %*% x)              # (X'X)^-1
# beta-hat = (X'X)^-1 X'y : the least squares estimates, computed directly
print(xtxi %*% t(x) %*% y)

# A better (more accurate) way: solve(A, b) solves A beta = b without forming the inverse.
# crossprod(x, y) computes t(x) %*% y more efficiently
print(solve(crossprod(x, x), crossprod(x, y)))


# ---- Pulling quantities out of the fitted model ----

# names(): lists what is stored inside the lm object
print(names(mdl))
# Store the summary so we can take things out of it too
mdls <- summary(mdl)
print(names(mdls))

# deviance(): the residual sum of squares (RSS); df.residual(): n - p
# sigma-hat = sqrt(RSS / (n - p))
print(sqrt(deviance(mdl) / df.residual(mdl)))
# The same value is stored in the summary as $sigma
print(mdls$sigma)

# $cov.unscaled is (X'X)^-1. Its diagonal times sigma-hat^2 gives the variances of the
# coefficients, so sqrt(diag()) * sigma-hat gives their standard errors
xtxi <- mdls$cov.unscaled
# diag(): takes the diagonal of a matrix
print(sqrt(diag(xtxi)) * mdls$sigma)
# The same standard errors from the coefficient table (column 2)
print(mdls$coefficients[, 2])

# R^2 = 1 - RSS / TSS, where TSS = sum of (y - mean(y))^2
print(1 - deviance(mdl) / sum((y - mean(y))^2))
print(mdls$r.squared)

# Other extractors worth knowing: residuals(), fitted(), coef()
print(head(residuals(mdl)))   # observed - fitted
print(head(fitted(mdl)))      # fitted (predicted) values
print(coef(mdl))              # the estimates

# What it shows (numbers from the book):
# - Fitted model: Species = 7.07 - 0.024 Area + 0.319 Elevation + 0.009 Nearest
#   - 0.241 Scruz - 0.075 Adjacent. Only Elevation (p = 3.8e-06) and Adjacent
#   (p = 0.0003) are clearly significant.
# - Residual SE = 60.98 on 24 df, R^2 = 0.766, adjusted R^2 = 0.717,
#   overall F = 15.7 on 5 and 24 df (p = 6.8e-07).
# - The direct matrix formula and solve(crossprod()) give exactly the same
#   coefficients as lm() here because the data are well-behaved. The direct
#   inverse is still a bad habit: it is slow and inaccurate when predictors are
#   strongly correlated (see exercise 3).
# - sigma-hat, the standard errors and R^2 computed by hand all match summary().
