# Checking the error assumptions (constant variance, normality), using the savings data

library(faraway)
data(savings)
g <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)


# ---- Constant variance ----
par(mfrow = c(1, 2))
# Residuals vs fitted: look for a funnel shape (non-constant variance) or curvature (wrong model)
plot(fitted(g), residuals(g), xlab = "Fitted", ylab = "Residuals")
abline(h = 0)                     # abline(h = 0): horizontal line at zero
# |Residuals| vs fitted folds the plot in half, making changes in spread easier to see
plot(fitted(g), abs(residuals(g)), xlab = "Fitted", ylab = "|Residuals|")
par(mfrow = c(1, 1))

# A rough numerical check: regress |residuals| on fitted values; a significant slope suggests
# the spread changes with the fitted value
print(summary(lm(abs(residuals(g)) ~ fitted(g))))

# Calibration: what do these plots look like when we KNOW the truth?
# Each line draws nine artificial residual plots (1:50 on the x-axis)
set.seed(1)
par(mfrow = c(3, 3))
for (i in 1:9) plot(1:50, rnorm(50))                         # constant variance
for (i in 1:9) plot(1:50, (1:50) * rnorm(50))                # strong non-constant variance
for (i in 1:9) plot(1:50, sqrt((1:50)) * rnorm(50))          # mild non-constant variance
for (i in 1:9) plot(1:50, cos((1:50) * pi / 25) + rnorm(50)) # non-linearity
par(mfrow = c(1, 1))

# Residuals against the predictors
par(mfrow = c(1, 2))
plot(savings$pop15, residuals(g), xlab = "Population under 15", ylab = "Residuals")
plot(savings$pop75, residuals(g), xlab = "Population over 75", ylab = "Residuals")
par(mfrow = c(1, 1))

# The pop15 plot shows two groups. var.test(): F-test that two groups have equal variance
print(var.test(residuals(g)[savings$pop15 > 35], residuals(g)[savings$pop15 < 35]))


# ---- Normality ----
par(mfrow = c(1, 2))
# qqnorm(): sorted residuals against normal quantiles; qqline(): line through the quartiles
qqnorm(residuals(g), ylab = "Residuals")
qqline(residuals(g))
# hist(): a histogram is a poor way to judge normality (the bins matter too much)
hist(residuals(g))
par(mfrow = c(1, 1))

# Calibration Q-Q plots: normal, lognormal (skewed), Cauchy (long tails), uniform (short tails)
par(mfrow = c(3, 3))
for (i in 1:9) qqnorm(rnorm(50))
for (i in 1:9) qqnorm(exp(rnorm(50)))
for (i in 1:9) qqnorm(rcauchy(50))
for (i in 1:9) qqnorm(runif(50))
par(mfrow = c(1, 1))

# shapiro.test(): formal test of H0 "the residuals are normal"
print(shapiro.test(residuals(g)))

# What it shows (numbers from the book):
# - No clear sign of non-constant variance in the fitted-value plots. The
#   |residual| regression has slope -0.20 with p = 0.093: not significant.
# - The artificial plots show how much pattern pure chance produces, which
#   helps avoid over-reading a real residual plot.
# - The residuals vs pop15 plot shows two groups of countries (pop15 above and
#   below 35%), and their variances differ: F = 2.79 on 22 and 26 df,
#   p = 0.014. The high-pop15 (younger) countries have more variable residuals.
# - The Q-Q plot is close to a straight line, and Shapiro-Wilk gives W = 0.987,
#   p = 0.85: no evidence against normality. The histogram is inconclusive.
# - Formal tests only detect the problem they were designed for, so plots
#   come first; tests are for confirming what a plot suggests.
