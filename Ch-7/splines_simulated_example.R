# Polynomials vs regression splines on simulated data with a known true curve
# True function: funky(x) = sin(2 * pi * x^3)^3, observed with N(0, 0.1^2) noise

library(splines)    # comes with R; splineDesign() builds B-spline basis functions

# function(x): defines our own function
funky <- function(x) sin(2 * pi * x^3)^3
x <- seq(0, 1, by = 0.01)          # 101 points from 0 to 1
set.seed(1)
y <- funky(x) + 0.1 * rnorm(101)
# matplot(..., type = "pl"): first column as points, second as a line
matplot(x, cbind(y, funky(x)), type = "pl", ylab = "y", pch = 18, lty = 1, col = 1,
        main = "Data and true curve")

# Orthogonal polynomials of degree 4 and 12
g4 <- lm(y ~ poly(x, 4))
g12 <- lm(y ~ poly(x, 12))
matplot(x, cbind(y, fitted(g4), fitted(g12)), type = "pll", ylab = "y", pch = 18,
        lty = c(1, 2, 1), col = 1, main = "Polynomials: degree 4 (dashed), 12 (solid)")
# (lty gives one line type per column: column 1 is the points, so 2 = dashed for degree 4)

# Cubic B-splines with 12 basis functions. The four repeated knots at each end are
# needed for the boundary; the interior knots are placed closer together where the
# curve bends most
knots <- c(0, 0, 0, 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1, 1, 1, 1)
# splineDesign(): one column per basis function, evaluated at each x (101 x 12)
bx <- splineDesign(knots, x)
# The B-splines add up to 1 at every x, so together with the intercept one column is
# redundant: lm() reports one NA coefficient, which is harmless here
gs <- lm(y ~ bx)
matplot(x, bx, type = "l", main = "The 12 B-spline basis functions")
matplot(x, cbind(y, fitted(gs)), type = "pl", ylab = "y", pch = 18, lty = 1, col = 1,
        main = "B-spline fit")
lines(x, funky(x), lty = 2, col = "red")   # the truth, for comparison

# How close is each fit to the truth? (root mean squared error against funky(x))
rmse_true <- function(f) sqrt(mean((f - funky(x))^2))
print(c(poly4 = rmse_true(fitted(g4)), poly12 = rmse_true(fitted(g12)),
        spline12 = rmse_true(fitted(gs))))

# What it shows:
# - Degree 4 clearly underfits. Degree 12 is better, but too wiggly on the
#   left and it misses the bend around x = 0.8: every point affects a
#   polynomial everywhere.
# - The B-spline fit with the same number of parameters (12) comes very close
#   to the truth: each basis function is non-zero only over a short stretch,
#   so the fit adapts locally while staying smooth. The RMSE against the
#   true curve confirms this.
