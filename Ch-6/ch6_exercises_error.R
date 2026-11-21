# =============================================================================
# Linear Models with R (Faraway), Chapter 6 exercises 1-5
# Weighted least squares and transformation (pipeline), GLS with AR(1) errors
# (divusa), lack-of-fit tests (salmonella, cars) and robust regression (stackloss).
# Packages: nlme and MASS (come with R), quantreg (install.packages("quantreg")).
# =============================================================================

library(faraway)
library(nlme)
library(MASS)


# =============================================================================
# Exercise 1: pipeline - Alaska pipeline defect depths (NIST)
#   Field = in-field ultrasonic measurement, Lab = laboratory measurement (more accurate),
#   Batch = measurement batch (ignored, as the question says)
# =============================================================================
data(pipeline)

# (a) Lab ~ Field and a check for non-constant variance
g <- lm(Lab ~ Field, pipeline)
print(summary(g))
par(mfrow = c(1, 2))
plot(Lab ~ Field, pipeline, main = "Lab vs Field")
abline(g)
plot(fitted(g), residuals(g), xlab = "Fitted", ylab = "Residuals", main = "Residuals vs fitted")
abline(h = 0)
par(mfrow = c(1, 1))
print(summary(lm(abs(residuals(g)) ~ fitted(g)))$coefficients)

# (b) Estimate how the variance depends on Field, using groups of 9 sorted cases
i <- order(pipeline$Field)                 # order(): the row order that sorts Field
npipe <- pipeline[i, ]
# gl(12, 9): a factor 1,1,...,1,2,...,12 (12 levels, each repeated 9 times) = 108 values;
# [-108] drops the last one because there are only 107 rows
ff <- gl(12, 9)[-108]
# split(): divides the values into the 12 groups; lapply() + unlist() gives one mean / variance per group
meanfield <- unlist(lapply(split(npipe$Field, ff), mean))
varlab <- unlist(lapply(split(npipe$Lab, ff), var))
print(cbind(meanfield, varlab))

# Model: var(Lab) = a0 * Field^a1, i.e. log(varlab) = log(a0) + a1 * log(meanfield).
# The last group has only 8 cases, so it is left out ([-12])
gv <- lm(log(varlab[-12]) ~ log(meanfield[-12]))
print(coef(gv))
a1 <- coef(gv)[2]
plot(log(meanfield), log(varlab), main = "log variance vs log mean")
abline(gv)

# Weights are 1 / variance, i.e. proportional to 1 / Field^a1 (a0 does not matter)
gw <- lm(Lab ~ Field, pipeline, weights = 1 / Field^a1)
print(summary(gw))

# (c) Transformations instead of weights: try every combination of
# identity / sqrt / log / inverse on Lab and Field
trans <- list(identity = function(z) z, sqrt = sqrt, log = log, inverse = function(z) 1 / z)
tab <- expand.grid(y = names(trans), x = names(trans), stringsAsFactors = FALSE)
tab$r.squared <- NA
tab$spread_p <- NA     # p-value of the |residual| ~ fitted slope (large = constant variance)
for (k in 1:nrow(tab)) {
  yy <- trans[[tab$y[k]]](pipeline$Lab)
  xx <- trans[[tab$x[k]]](pipeline$Field)
  fit <- lm(yy ~ xx)
  tab$r.squared[k] <- summary(fit)$r.squared
  tab$spread_p[k] <- summary(lm(abs(residuals(fit)) ~ fitted(fit)))$coefficients[2, 4]
}
print(tab[order(-tab$spread_p), ], digits = 3)

# log-log is the natural candidate to look at (multiplicative error)
glog <- lm(log(Lab) ~ log(Field), pipeline)
print(summary(glog))
par(mfrow = c(1, 2))
plot(log(Lab) ~ log(Field), pipeline, main = "log-log")
abline(glog)
plot(fitted(glog), residuals(glog), xlab = "Fitted", ylab = "Residuals", main = "log-log residuals")
abline(h = 0)
par(mfrow = c(1, 1))

# What stands out (pipeline):
# [numbers to be filled in from your run_all.R output]
# - (a) The residuals fan out as Field grows: deeper defects are measured
#   less precisely, so the variance is not constant.
# - (b) The slope a1 says how fast var(Lab) grows with Field; the WLS fit
#   gives the deep (noisy) defects less weight.
# - (c) The table ranks the transformation pairs by how constant the spread
#   is; the chosen pair should also look linear.


# =============================================================================
# Exercise 2: divusa - divorce on unemployed, femlab, marriage, birth, military
# =============================================================================
data(divusa)
g <- lm(divorce ~ unemployed + femlab + marriage + birth + military, divusa)
print(summary(g))
r <- residuals(g)
nr <- length(r)

# (a) Two graphical checks for correlated errors
par(mfrow = c(1, 2))
plot(divusa$year, r, type = "b", xlab = "year", ylab = "Residuals", main = "Residuals over time")
abline(h = 0)
plot(r[-nr], r[-1], xlab = expression(hat(epsilon)[i]),
     ylab = expression(hat(epsilon)[i + 1]), main = "Successive residuals")
par(mfrow = c(1, 1))

# (b) GLS with AR(1) errors, fitted by maximum likelihood (method = "ML")
gg <- gls(divorce ~ unemployed + femlab + marriage + birth + military,
          correlation = corAR1(form = ~ year), data = divusa, method = "ML")
print(summary(gg))
print(intervals(gg))
# Compare the p-values: OLS vs GLS
print(round(cbind(OLS_p = summary(g)$coefficients[, 4],
                  GLS_p = summary(gg)$tTable[, 4]), 4))

# What stands out (divusa, correlated errors):
# [numbers to be filled in from your run_all.R output]
# - (a) Long runs of same-sign residuals and a strong upward trend in the
#   successive-residual plot: clear positive serial correlation.
# - (b) phi is the estimated lag-1 correlation of the errors; if its 95% CI
#   excludes 0, it is significant. Compare the p-value columns to see which
#   predictors stop (or start) being significant once the correlation is
#   allowed for.
# - (c) Why: divorce rates move slowly. Social attitudes, divorce laws
#   (no-fault divorce spread through the 1970s) and other unmeasured
#   influences persist from one year to the next, so this year's error
#   carries into next year's.


# =============================================================================
# Exercise 3: salmonella - revertant colonies at 6 doses of quinoline (3 plates each)
# =============================================================================
data(salmonella)
print(salmonella)
g <- lm(colonies ~ log(dose + 1), salmonella)
print(summary(g))
# Lack of fit: one mean per dose (the replicates give a model-free error estimate)
ga <- lm(colonies ~ factor(dose), salmonella)
print(anova(g, ga))
plot(colonies ~ log(dose + 1), salmonella, main = "salmonella")
abline(g)
points(log(salmonella$dose + 1), fitted(ga), pch = 18, col = "red")   # group means

# What stands out (salmonella):
# [numbers to be filled in from your run_all.R output]
# - A small anova() p-value means the straight line in log(dose + 1) does not
#   fit the dose means. The plot shows how: the counts rise and then fall
#   again at the highest doses, which a straight line cannot follow.


# =============================================================================
# Exercise 4: cars - stopping distance vs speed (comes with R)
# =============================================================================
data(cars)
g <- lm(dist ~ speed, cars)
print(summary(g))
# table(): how many cars share each speed - repeated speeds are the replicates
print(table(cars$speed))
ga <- lm(dist ~ factor(speed), cars)
print(anova(g, ga))

# What stands out (cars):
# [numbers to be filled in from your run_all.R output]
# - A large p-value: no evidence of lack of fit for the straight line, given
#   the scatter between cars at the same speed. This does not prove the line
#   is right (physics suggests distance grows like speed^2); the replicates
#   are just too noisy to show the curvature.


# =============================================================================
# Exercise 5: stackloss - 21 days of a plant oxidising ammonia (comes with R)
#   Air.Flow, Water.Temp, Acid.Conc.; stack.loss = ammonia lost (the response)
# =============================================================================
data(stackloss)

# (a) Least squares
g_ls <- lm(stack.loss ~ ., stackloss)
print(summary(g_ls))
# (b) Least absolute deviations
if (requireNamespace("quantreg", quietly = TRUE)) {
  g_lad <- quantreg::rq(stack.loss ~ ., data = stackloss)
  print(summary(g_lad))
  lad_coef <- coef(g_lad)
} else {
  cat("Install quantreg for the LAD fit\n")
  lad_coef <- rep(NA, 4)
}
# (c) Huber
g_hub <- rlm(stack.loss ~ ., stackloss)
print(summary(g_hub))
# (d) Least trimmed squares, exact search (choose(21, 4) = 5985 subsets, quick)
g_lts <- ltsreg(stack.loss ~ ., stackloss, nsamp = "exact")
print(coef(g_lts))

print(round(rbind(LS = coef(g_ls), LAD = lad_coef, Huber = coef(g_hub), LTS = coef(g_lts)), 3))
# The Huber weights show which days rlm() down-weighted
print(round(g_hub$w, 2))

# Diagnostics on the least squares fit
n <- nrow(stackloss); p <- 4
diag_tab <- data.frame(day = 1:n,
                       leverage = round(hatvalues(g_ls), 3),
                       jackknife = round(rstudent(g_ls), 2),
                       cook = round(cooks.distance(g_ls), 3))
print(diag_tab[order(-diag_tab$cook), ][1:6, ])
print(c(leverage_cutoff = 2 * p / n, bonferroni = abs(qt(0.05 / (2 * n), n - p - 1)),
        cook_cutoff = 4 / n))
par(mfrow = c(1, 2))
halfnorm(cooks.distance(g_ls), 3, ylab = "Cook's distances")
halfnorm(hatvalues(g_ls), 3, ylab = "Leverages")
par(mfrow = c(1, 1))

# Remove the cases flagged as influential (Cook's distance above 4/n) and refit by least squares
keep_days <- cooks.distance(g_ls) <= 4 / n     # TRUE for the days we keep
print(which(!keep_days))                         # the days removed
g_ls2 <- lm(stack.loss ~ ., stackloss[keep_days, ])
print(summary(g_ls2))
print(round(rbind(LS_all = coef(g_ls), LS_dropped = coef(g_ls2), LTS = coef(g_lts)), 3))

# What stands out (stackloss):
# [numbers to be filled in from your run_all.R output]
# - Compare the four rows: the robust fits usually agree with each other and
#   differ from least squares in the Acid.Conc. and Water.Temp coefficients.
# - The Huber weights below 1 point to the days the robust fit distrusts;
#   days 1, 3, 4 and 21 are the classic suspects for these data.
# - After removing the influential days, least squares moves towards the
#   robust estimates.
