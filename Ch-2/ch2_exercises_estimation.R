# =============================================================================
# Linear Models with R (Faraway), Chapter 2 exercises 1-5
# Fitting linear models with lm(), reading the output, and computing the
# estimates by hand. Each section: fit, extract what the question asks for,
# and comment on it.
# =============================================================================

library(faraway)   # library(): loads the faraway package, which contains the datasets


# =============================================================================
# Exercise 1: teengamb - gambling expenditure of 47 British teenagers
#   sex 0 = male, 1 = female; status, income (pounds/week), verbal (score out of 12)
#   gamble = expenditure on gambling, pounds per year (the response)
# =============================================================================
data(teengamb)

# lm(): least squares fit of gamble on the four predictors
g1 <- lm(gamble ~ sex + status + income + verbal, data = teengamb)
# summary(): the full regression output ("present the output")
print(summary(g1))

# (a) R^2 is the proportion of variation in the response explained by the model
print(summary(g1)$r.squared * 100)   # as a percentage

# (b) Largest positive residual
# residuals(): observed - fitted.  which.max(): position of the largest value
print(which.max(residuals(g1)))
print(max(residuals(g1)))

# (c) Mean and median of the residuals
print(mean(residuals(g1)))     # always 0 (up to rounding) when the model has an intercept
print(median(residuals(g1)))

# (d) Correlation of the residuals with the fitted values
# cor(): correlation between two vectors; fitted(): the fitted values
print(cor(residuals(g1), fitted(g1)))

# (e) Correlation of the residuals with income
print(cor(residuals(g1), teengamb$income))

# (f) sex is 0 for males and 1 for females, so its coefficient is (female - male).
# coef(): the estimates; ["sex"] picks the one called sex. Male - female = -coef
print(-coef(g1)["sex"])

# What stands out (teengamb):
# [numbers to be filled in from your run_all.R output]
# - (a) R^2 x 100 is the percentage of the variation in gamble explained.
# - (b) The case number printed by which.max() is the teenager the model
#   under-predicts the most.
# - (c) The mean residual is 0 by construction: least squares with an intercept
#   forces the residuals to sum to zero. The median is negative because the
#   residuals are right-skewed (a few large positive ones), like gamble itself.
# - (d), (e) Both correlations are 0 (up to rounding error). This is not a
#   finding about the data: the residual vector is orthogonal to every column
#   of X and to the fitted values, so these correlations are always zero.
# - (f) Holding status, income and verbal fixed, a male is predicted to spend
#   -coef(sex) pounds a year more on gambling than a female.


# =============================================================================
# Exercise 2: uswages - weekly wages of 2000 US male workers (CPS 1988)
#   wage (dollars/week), educ and exper (years of education / experience)
# =============================================================================
data(uswages)

# Weekly wage on the original scale
g2 <- lm(wage ~ educ + exper, data = uswages)
print(summary(g2))

# Same model with log(wage) as the response
g2log <- lm(log(wage) ~ educ + exper, data = uswages)
print(summary(g2log))

# On the log scale a coefficient b multiplies the wage by exp(b) per unit change.
# exp(b) - 1 is the proportional change; x 100 gives a percentage
print((exp(coef(g2log)["educ"]) - 1) * 100)

# What stands out (uswages):
# [numbers to be filled in from your run_all.R output]
# - Original scale: one extra year of education is associated with an increase
#   of coef(educ) dollars in weekly wage, holding experience fixed. This says
#   the same dollar gain whether the worker starts from a low or a high wage.
# - Log scale: one extra year of education multiplies the expected wage by
#   exp(coef(educ)), i.e. raises it by the percentage printed above, holding
#   experience fixed.
# - The log interpretation is more natural: raises and pay differences are
#   usually thought of in percentages, and Chapter 1 showed wage is strongly
#   right-skewed while log(wage) is roughly symmetric.
# - Note: as found in Chapter 1, 33 workers have negative (impossible)
#   experience. They are left in here because the exercise uses the data as given.


# =============================================================================
# Exercise 3: computing the coefficients two ways, lm() vs the direct formula
#   x = 1:20, y = x + noise. Fit polynomials of increasing degree and find the
#   degree at which (X'X)^-1 X'y can no longer be computed.
# =============================================================================
set.seed(123)            # set.seed(): makes the random numbers reproducible
x <- 1:20
y <- x + rnorm(20)       # rnorm(20): 20 draws from N(0, 1)

# Data frame to collect the results, one row per degree
results <- data.frame(degree = 1:12, lm_ok = NA, direct_ok = NA, max_diff = NA)

for (d in 1:12) {
  # poly(x, d, raw = TRUE): the columns x, x^2, ..., x^d
  # (the same as writing x + I(x^2) + ... + I(x^d) in the formula)
  fit <- lm(y ~ poly(x, d, raw = TRUE))
  # anyNA(coef()): TRUE if lm() had to drop a column because it was too collinear
  results$lm_ok[d] <- !anyNA(coef(fit))

  X <- model.matrix(fit)   # the X matrix lm() used: 1, x, x^2, ..., x^d
  # tryCatch(): run the direct formula; if solve() throws an error, return NULL instead
  b_direct <- tryCatch(solve(t(X) %*% X) %*% t(X) %*% y,
                       error = function(e) NULL)
  results$direct_ok[d] <- !is.null(b_direct)
  if (!is.null(b_direct) && results$lm_ok[d]) {
    # Compare the fitted values from the two methods
    results$max_diff[d] <- max(abs(X %*% b_direct - fitted(fit)))
  }
}
print(results)

# The error solve() gives at the first failing degree:
d_fail <- which(!results$direct_ok)[1]
if (!is.na(d_fail)) {
  X <- model.matrix(~ poly(x, d_fail, raw = TRUE))
  # conditionMessage(): the text of the error message
  print(tryCatch(solve(t(X) %*% X), error = function(e) conditionMessage(e)))
}

# What stands out (polynomial computation):
# [degree to be filled in from your run_all.R output]
# - The direct formula fails (solve() reports the system is "computationally
#   singular") at the degree shown in the direct_ok column. The powers of x
#   grow so fast (20^d) that X'X has a huge condition number and cannot be
#   inverted in double precision.
# - lm() keeps working for longer because it uses a QR decomposition of X
#   instead of inverting X'X. When it finally gives up it drops a column (NA
#   coefficient) rather than giving a wrong answer.
# - Before failing outright, the direct method already loses accuracy: see
#   max_diff growing with the degree.


# =============================================================================
# Exercise 4: prostate - 97 men about to receive a radical prostatectomy
#   Start with lpsa ~ lcavol, then add lweight, svi, lbph, age, lcp, pgg45,
#   gleason one at a time, recording the residual SE and R^2 each time.
# =============================================================================
data(prostate)

vars <- c("lcavol", "lweight", "svi", "lbph", "age", "lcp", "pgg45", "gleason")
trend <- data.frame(added = vars, sigma = NA, r.squared = NA)

for (k in seq_along(vars)) {
  # reformulate(): builds the formula lpsa ~ (first k variables) from text
  f <- reformulate(vars[1:k], response = "lpsa")
  fit <- lm(f, data = prostate)
  trend$sigma[k] <- summary(fit)$sigma          # residual standard error
  trend$r.squared[k] <- summary(fit)$r.squared
}
print(trend)

# par(mfrow = c(1, 2)): two plots side by side
par(mfrow = c(1, 2))
# type = "b": points joined by lines; xaxt = "n" suppresses the default x-axis
plot(1:8, trend$sigma, type = "b", xaxt = "n", xlab = "variable added",
     ylab = "residual SE", main = "Residual SE")
# axis(1, ...): draws our own x-axis with the variable names; las = 2 turns the labels
axis(1, at = 1:8, labels = vars, las = 2, cex.axis = 0.8)
plot(1:8, trend$r.squared, type = "b", xaxt = "n", xlab = "variable added",
     ylab = "R-squared", main = "R-squared")
axis(1, at = 1:8, labels = vars, las = 2, cex.axis = 0.8)
par(mfrow = c(1, 1))

# What stands out (prostate, adding variables):
# [numbers to be filled in from your run_all.R output]
# - R^2 can only go up as variables are added, because a bigger model can
#   always do at least as well. That is why R^2 alone cannot choose a model.
# - The residual SE can go up or down: it divides the RSS by n - p, so a
#   variable that hardly reduces the RSS makes it worse.
# - The big gains come from the first few variables (lcavol, then lweight and
#   svi); after that the curves flatten.


# =============================================================================
# Exercise 5: prostate - regression of lpsa on lcavol AND of lcavol on lpsa
# =============================================================================
plot(lpsa ~ lcavol, data = prostate, main = "Two regression lines")

g_y <- lm(lpsa ~ lcavol, data = prostate)   # lpsa = a1 + b1 * lcavol
g_x <- lm(lcavol ~ lpsa, data = prostate)   # lcavol = a2 + b2 * lpsa
# abline(): draws a straight line; given a model it uses its intercept and slope
abline(g_y, col = "red", lwd = 2)

# Rewrite lcavol = a2 + b2 * lpsa as lpsa = -a2/b2 + (1/b2) * lcavol to draw it on these axes
a2 <- coef(g_x)[1]
b2 <- coef(g_x)[2]
abline(-a2 / b2, 1 / b2, col = "blue", lwd = 2)

# Intersection: solve a1 + b1 * x = -a2/b2 + (1/b2) * x for x
a1 <- coef(g_y)[1]
b1 <- coef(g_y)[2]
x_int <- (-a2 / b2 - a1) / (b1 - 1 / b2)
y_int <- a1 + b1 * x_int
print(c(lcavol = unname(x_int), lpsa = unname(y_int)))
# Compare with the two means
print(colMeans(prostate[, c("lcavol", "lpsa")]))
# points(): marks the intersection; pch = 19 is a filled circle
points(x_int, y_int, pch = 19, cex = 1.5)
legend("topleft", legend = c("lpsa on lcavol", "lcavol on lpsa"),
       col = c("red", "blue"), lwd = 2)

# What stands out (two regression lines):
# - The two lines are different: regressing y on x minimises vertical
#   distances, x on y minimises horizontal ones. The x-on-y line is steeper.
# - They intersect at exactly (mean(lcavol), mean(lpsa)), because every least
#   squares line with an intercept passes through the point of the means.
#   [values to be filled in from your run_all.R output]
