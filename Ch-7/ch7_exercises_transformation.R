# =============================================================================
# Linear Models with R (Faraway), Chapter 7 exercises 1-5
# Trend, polynomial, broken-stick and spline models for aatemp; transformations
# and goodness of fit for cornnit, ozone, pressure and trees.
# =============================================================================

library(faraway)
library(MASS)      # boxcox()
library(nlme)      # gls()
library(splines)   # splineDesign()


# =============================================================================
# Exercise 1: aatemp - annual mean temperature (F) in Ann Arbor, Michigan
# =============================================================================
data(aatemp)
print(range(aatemp$year))
print(nrow(aatemp))

# (a) Linear trend
g_lin <- lm(temp ~ year, aatemp)
print(summary(g_lin))
plot(temp ~ year, aatemp, main = "Ann Arbor annual mean temperature")
abline(g_lin, lwd = 2)

# (b) AR(1) correlation between successive years, estimated together with the trend
g_ar <- gls(temp ~ year, correlation = corAR1(form = ~ year), data = aatemp)
print(summary(g_ar))
print(intervals(g_ar))

# (c) Degree-10 polynomial, then backward elimination of the top degree.
# poly(year, d) uses orthogonal polynomials, which avoids the numerical trouble raw powers
# of numbers like 2000 would cause (2000^10 = 1e33). The p-value of the last term is the
# test of the highest degree.
d <- 10
repeat {
  g_poly <- lm(temp ~ poly(year, d), aatemp)
  ct <- summary(g_poly)$coefficients
  p_top <- ct[nrow(ct), "Pr(>|t|)"]        # p-value of the degree-d term (last row)
  cat("degree", d, ": p-value of top term =", signif(p_top, 3), "\n")
  if (p_top < 0.05 || d == 1) break         # stop when the top term is significant
  d <- d - 1
}
print(summary(g_poly))
grid <- data.frame(year = seq(min(aatemp$year), 2020, by = 1))
lines(grid$year, predict(g_poly, grid), col = "red", lwd = 2)
# Prediction for 2020 (outside the data: an extrapolation)
print(predict(g_poly, data.frame(year = 2020), interval = "prediction"))

# (d) Constant until 1930, then a linear trend: one hockey-stick term
rhs <- function(x) ifelse(x < 1930, 0, x - 1930)
g_claim <- lm(temp ~ rhs(year), aatemp)
print(summary(g_claim))
# Compare with a broken stick that ALSO allows a slope before 1930: the claim is that
# this before-1930 slope is zero
lhs <- function(x) ifelse(x < 1930, 1930 - x, 0)
g_bs <- lm(temp ~ lhs(year) + rhs(year), aatemp)
print(anova(g_claim, g_bs))
lines(grid$year, predict(g_claim, grid), col = "blue", lwd = 2, lty = 2)

# (e) Cubic B-spline with six basis functions, knots evenly spaced on the range.
# Six cubic basis functions need 6 + 4 = 10 knots: 4 repeated at each end plus 2 interior
# knots at 1/3 and 2/3 of the way along the range
a <- min(aatemp$year); b <- max(aatemp$year)
knots <- c(rep(a, 4), a + (b - a) / 3, a + 2 * (b - a) / 3, rep(b, 4))
bx <- splineDesign(knots, aatemp$year)
print(dim(bx))                       # 115 rows, 6 basis functions
g_spl <- lm(temp ~ bx, aatemp)       # one coefficient will be NA: the basis sums to 1
print(summary(g_spl)$r.squared)
o <- order(aatemp$year)
lines(aatemp$year[o], fitted(g_spl)[o], col = "darkgreen", lwd = 2, lty = 3)
legend("topleft", legend = c("straight line", paste("polynomial, degree", d),
                             "constant then linear", "cubic spline"),
       col = c("black", "red", "blue", "darkgreen"), lty = c(1, 1, 2, 3), lwd = 2, cex = 0.8)
# The straight line is a special case of the cubic spline, so the F-test compares them
print(anova(g_lin, g_spl))

# What stands out (aatemp):
# [numbers to be filled in from your run_all.R output]
# - (a) The slope of temp on year and its p-value answer "is there a linear trend?".
# - (b) phi is the estimated year-to-year correlation; with positive
#   correlation the GLS standard error of the slope is larger, so the trend
#   may become less (or not) significant.
# - (c) The degree reached by backward elimination and the 2020 prediction;
#   a polynomial extrapolated beyond the data can head off in any direction,
#   so the 2020 value should not be trusted.
# - (d) A large anova() p-value means the data are consistent with no trend
#   before 1930; the rhs(year) slope then gives the warming rate after 1930.
# - (e) A small anova() p-value means the spline captures real curvature the
#   straight line misses.


# =============================================================================
# Exercise 2: cornnit - corn yield (bushels/acre) vs nitrogen (lb/acre), Wisconsin 1994
# =============================================================================
data(cornnit)
print(table(cornnit$nitrogen))       # repeated nitrogen levels allow a lack-of-fit test
plot(yield ~ nitrogen, cornnit, main = "cornnit")

# Pure-error model: one mean per nitrogen level
g_pe <- lm(yield ~ factor(nitrogen), cornnit)

# Candidate models; nitrogen has zeros, so log() uses nitrogen + 1
cands <- list(
  linear    = lm(yield ~ nitrogen, cornnit),
  quadratic = lm(yield ~ nitrogen + I(nitrogen^2), cornnit),
  sqrt      = lm(yield ~ sqrt(nitrogen), cornnit),
  log       = lm(yield ~ log(nitrogen + 1), cornnit)
)
# For each candidate: R^2, residual SE and the lack-of-fit p-value against the pure-error model
lof <- t(sapply(cands, function(m) c(r.squared = summary(m)$r.squared,
                                      sigma = summary(m)$sigma,
                                      lack_of_fit_p = anova(m, g_pe)[2, "Pr(>F)"])))
print(round(lof, 4))

grid <- data.frame(nitrogen = seq(0, max(cornnit$nitrogen), length.out = 100))
cols <- c("black", "red", "blue", "darkgreen")
for (k in seq_along(cands)) lines(grid$nitrogen, predict(cands[[k]], grid), col = cols[k], lwd = 2)
points(cornnit$nitrogen, fitted(g_pe), pch = 18, cex = 1.5)    # level means
legend("bottomright", legend = names(cands), col = cols, lwd = 2)

# Residual plot for the best candidate (largest lack-of-fit p-value)
best <- cands[[which.max(lof[, "lack_of_fit_p"])]]
plot(fitted(best), residuals(best), main = paste("Residuals:", names(cands)[which.max(lof[, 3])]),
     xlab = "Fitted", ylab = "Residuals")
abline(h = 0)

# What stands out (cornnit):
# [numbers to be filled in from your run_all.R output]
# - Yield rises steeply with the first nitrogen and then levels off
#   (diminishing returns), so the straight line should fail the lack-of-fit
#   test, while a concave shape (sqrt, log or quadratic) should pass.


# =============================================================================
# Exercise 3: ozone - daily ozone in the Los Angeles basin
#   O3 (response), temp, humidity, ibh (inversion base height), and others
# =============================================================================
data(ozone)
g <- lm(O3 ~ temp + humidity + ibh, ozone)
print(summary(g))
bc <- boxcox(g, lambda = seq(-0.5, 1.5, by = 0.01), plotit = TRUE)
lam <- bc$x[which.max(bc$y)]
print(lam)
print(range(bc$x[bc$y > max(bc$y) - qchisq(0.95, 1) / 2]))   # 95% CI for lambda

# Refit with the nearest interpretable power (here assumed to be the square root; check
# that 0.5 is inside the interval printed above) and compare the residual plots
g_sqrt <- lm(sqrt(O3) ~ temp + humidity + ibh, ozone)
par(mfrow = c(1, 2))
plot(fitted(g), residuals(g), main = "O3", xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
plot(fitted(g_sqrt), residuals(g_sqrt), main = "sqrt(O3)", xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
par(mfrow = c(1, 1))

# What stands out (ozone):
# [numbers to be filled in from your run_all.R output]
# - The best lambda and its 95% CI. If 1 is outside the interval, the
#   response should be transformed; round lambda to a simple power inside
#   the interval.


# =============================================================================
# Exercise 4: pressure - vapour pressure of mercury (mm) vs temperature (deg C), 19 points
#   (comes with R)
# =============================================================================
data(pressure)
g_raw <- lm(pressure ~ temperature, pressure)
# Physics (Clausius-Clapeyron): log(pressure) is roughly linear in 1 / absolute temperature
g_cc <- lm(log(pressure) ~ I(1 / (temperature + 273.15)), pressure)
g_log <- lm(log(pressure) ~ temperature, pressure)
print(summary(g_cc))
print(c(raw = summary(g_raw)$r.squared, log_linear = summary(g_log)$r.squared,
        clausius = summary(g_cc)$r.squared))

par(mfrow = c(2, 2))
plot(pressure ~ temperature, pressure, main = "Raw data")
plot(fitted(g_raw), residuals(g_raw), main = "pressure ~ temperature", xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
plot(fitted(g_log), residuals(g_log), main = "log(pressure) ~ temperature", xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
plot(fitted(g_cc), residuals(g_cc), main = "log(pressure) ~ 1/(T + 273.15)", xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
par(mfrow = c(1, 1))

# What stands out (pressure):
# [numbers to be filled in from your run_all.R output]
# - The raw relationship is extremely curved (pressure grows roughly
#   exponentially). log(pressure) on temperature is still curved.
#   log(pressure) on 1 / absolute temperature is almost exactly linear, as the
#   physics predicts: the residuals should be tiny and patternless.


# =============================================================================
# Exercise 5: trees - Volume, Girth (actually diameter, in) and Height (ft) of 31 black
#   cherry trees (comes with R)
# =============================================================================
data(trees)
g_lin <- lm(Volume ~ Girth + Height, trees)
bc <- boxcox(g_lin, lambda = seq(-0.5, 1.5, by = 0.01), plotit = TRUE)
print(bc$x[which.max(bc$y)])
print(range(bc$x[bc$y > max(bc$y) - qchisq(0.95, 1) / 2]))

# Geometry: a tree is roughly a cone or cylinder, Volume ~ c * Girth^2 * Height,
# so log(Volume) = log(c) + 2 log(Girth) + 1 log(Height)
g_log <- lm(log(Volume) ~ log(Girth) + log(Height), trees)
print(summary(g_log))
print(confint(g_log))     # do the intervals include 2 and 1?

par(mfrow = c(1, 2))
plot(fitted(g_lin), residuals(g_lin), main = "Volume ~ Girth + Height", xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
plot(fitted(g_log), residuals(g_log), main = "log-log model", xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
par(mfrow = c(1, 1))

# What stands out (trees):
# [numbers to be filled in from your run_all.R output]
# - The linear model's residuals are curved. Box-Cox should point to a power
#   near 1/3 (cube root), which makes sense for a volume.
# - The log-log model has coefficients close to 2 for Girth and 1 for Height,
#   exactly what the cone/cylinder formula predicts, and patternless residuals.
