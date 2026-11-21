# Multiple outliers hiding each other, using the star data
# 47 stars in the CYG OB1 cluster: temp = log(surface temperature), light = log(light intensity)

library(faraway)
data(star)

plot(star$temp, star$light, xlab = "log(Temperature)", ylab = "log(Light Intensity)")

# Least squares line through all 47 stars
ga <- lm(light ~ temp, star)
abline(ga)

# Does the outlier test find the four stars on the left?
print(range(rstudent(ga)))
# Bonferroni critical value for comparison
print(qt(0.05 / (2 * nrow(star)), df.residual(ga) - 1))

# Refit without the four giants (temp below 3.6) and draw it dashed
ga2 <- lm(light ~ temp, data = star, subset = (temp > 3.6))
abline(ga2, lty = 2)
legend("bottomleft", legend = c("all stars", "without the 4 giants"), lty = c(1, 2))
print(coef(ga))
print(coef(ga2))

# What it shows (numbers from the book):
# - The four stars at the upper left (red giants) pull the least squares line
#   so far towards themselves that the line no longer follows the main group.
# - The jackknife residuals only range from -2.05 to 1.91, so the outlier test
#   finds nothing, although the outliers are obvious in the plot. Several
#   outliers together can hide each other.
# - Without the giants the slope changes sign and the line follows the bulk of
#   the data. In more than two dimensions such points are much harder to see,
#   which is the case for robust regression (Chapter 6).
