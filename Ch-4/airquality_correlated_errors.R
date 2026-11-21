# Checking for correlated errors, using the airquality data (153 consecutive days, New York, 1973)
# Columns: Ozone (ppb), Solar.R (radiation), Wind (mph), Temp (F), Month, Day
# airquality comes with R itself (datasets package), so no library() is needed for it.

data(airquality)
print(head(airquality))       # there are missing values (NA)
# pairs(..., panel = panel.smooth): scatterplot matrix with a smooth trend in each panel
pairs(airquality, panel = panel.smooth)

# na.action = na.exclude: rows with NA are left out of the fit, but residuals() and
# fitted() still return 153 values, with NA in the dropped rows (keeps the time order)
g <- lm(Ozone ~ Solar.R + Wind + Temp, airquality, na.action = na.exclude)
print(summary(g))

par(mfrow = c(1, 2))
plot(fitted(g), residuals(g), xlab = "Fitted", ylab = "Residuals", main = "Ozone")
# log() of the response fixes the fan shape and the curvature
gl <- lm(log(Ozone) ~ Solar.R + Wind + Temp, airquality, na.action = na.exclude)
plot(fitted(gl), residuals(gl), xlab = "Fitted", ylab = "Residuals", main = "log(Ozone)")
par(mfrow = c(1, 1))

# ---- Graphical checks for serial correlation ----
par(mfrow = c(1, 2))
# Index plot: residuals in time order. Long runs above/below 0 suggest positive correlation
plot(residuals(gl), ylab = "Residuals", main = "Residuals in time order")
abline(h = 0)
# Successive residuals: residual i (x-axis) vs residual i+1 (y-axis).
# [-153] drops the last value, [-1] drops the first, so the pairs line up.
# expression(hat(epsilon)[i]): writes the axis label as the maths symbol epsilon-hat_i
plot(residuals(gl)[-153], residuals(gl)[-1],
     xlab = expression(hat(epsilon)[i]), ylab = expression(hat(epsilon)[i + 1]),
     main = "Successive residuals")
par(mfrow = c(1, 1))

# Regression of each residual on the previous one; -1 removes the intercept
# (residuals have mean zero). Pairs with an NA are dropped automatically
print(summary(lm(residuals(gl)[-1] ~ -1 + residuals(gl)[-153])))

# ---- Durbin-Watson test ----
# DW = sum of squared successive differences / sum of squared residuals; about 2 means
# no correlation. By hand, using the complete cases (as the book's test does):
gna <- lm(Ozone ~ Solar.R + Wind + Temp, data = na.omit(airquality))   # na.omit(): drops NA rows
r <- residuals(gna)
# diff(r): r[2] - r[1], r[3] - r[2], ...
print(sum(diff(r)^2) / sum(r^2))
# dwtest() in the lmtest package also gives a p-value.
# requireNamespace(): TRUE if the package is installed (install with install.packages("lmtest"))
if (requireNamespace("lmtest", quietly = TRUE)) {
  print(lmtest::dwtest(Ozone ~ Solar.R + Wind + Temp, data = na.omit(airquality)))
} else {
  cat("Install the lmtest package for the Durbin-Watson p-value\n")
}

# What it shows (numbers from the book):
# - 111 complete days are used, so the model has 107 residual df.
#   R^2 = 0.61; all three predictors are significant.
# - The untransformed fit shows non-constant variance and curvature;
#   log(Ozone) largely fixes both.
# - Neither the index plot nor the successive-residual plot shows obvious
#   serial correlation; there is one outlier (it appears twice in the second plot).
# - The lag-1 regression slope is 0.11 (SE 0.105, p = 0.30): not significant.
# - Durbin-Watson: DW = 1.94, p = 0.33, no evidence of correlation. Treat it
#   with caution: dropping the days with missing values breaks the time
#   sequence the test relies on.
