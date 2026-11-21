# Errors in the predictors and the SIMEX method, using the cars data (50 cars from the 1920s)
# Columns: speed (mph), dist (stopping distance, ft). cars comes with R itself.

data(cars)
plot(dist ~ speed, cars, ylab = "distance", main = "Stopping distance vs speed")
g <- lm(dist ~ speed, cars)
print(summary(g))
abline(g)

# Add extra random error to speed and refit: the slope shrinks towards zero.
# I(speed + rnorm(50)): speed plus N(0, 1) noise, computed inside the formula
set.seed(1)
ge1 <- lm(dist ~ I(speed + rnorm(50)), cars)
print(coef(ge1))
abline(ge1, lty = 2)
ge2 <- lm(dist ~ I(speed + 2 * rnorm(50)), cars)      # noise SD 2
print(coef(ge2))
abline(ge2, lty = 3)
ge5 <- lm(dist ~ I(speed + 5 * rnorm(50)), cars)      # noise SD 5
print(coef(ge5))
abline(ge5, lty = 4)
legend("topleft", legend = c("original", "+ SD 1 noise", "+ SD 2 noise", "+ SD 5 noise"),
       lty = 1:4)


# ---- SIMEX: SIMulate the effect of extra error, then EXtrapolate back to no error ----
# Suppose speed was measured with error variance 0.5. Add extra error with variance
# 0.1, 0.2, ..., 0.5, 1000 times each, and record the slope every time.
# rep(1:5/10, each = 1000): 0.1 repeated 1000 times, then 0.2 1000 times, ...
vv <- rep(1:5 / 10, each = 1000)
slopes <- numeric(5000)
for (i in 1:5000) {
  # sqrt(vv[i]) * rnorm(50): noise with variance vv[i]; $coefficients[2] is the slope
  slopes[i] <- lm(dist ~ I(speed + sqrt(vv[i]) * rnorm(50)), cars)$coefficients[2]
}

# matrix(slopes, nrow = 1000): one column per variance setting; colMeans(): mean slope per column
betas <- c(coef(g)[2], colMeans(matrix(slopes, nrow = 1000)))
# Total error variance at each setting: the assumed 0.5 plus what we added
variances <- c(0, 1:5 / 10) + 0.5
print(cbind(variances, betas))
plot(variances, betas, xlim = c(0, 1), ylim = c(3.86, 4),
     main = "SIMEX extrapolation", ylab = "mean slope", xlab = "error variance")

# Fit a straight line through the points and read off its value at variance 0
gv <- lm(betas ~ variances)
print(coef(gv))
abline(gv, lty = 2)
points(0, coef(gv)[1], pch = 3, cex = 1.5)   # the SIMEX estimate of the slope

# What it shows (numbers from the book; the noisy fits are random):
# - Least squares: dist = -17.58 + 3.93 speed, R^2 = 0.65.
# - Adding noise to speed flattens the line: slopes of about 3.76, 3.17 and
#   1.87 for noise SD 1, 2 and 5 in the book. Measurement error in x biases
#   the slope towards zero (attenuation).
# - SIMEX: the mean slope falls steadily as the error variance grows.
#   Extrapolating the straight line back to variance 0 gives about 4.0, a
#   little steeper than the least squares 3.93, as an estimate of the slope
#   without measurement error.
# - A better extrapolation model than a straight line could be used (Cook and
#   Stefanski, 1994).
