# =============================================================================
# Linear Models with R (Faraway), Chapter 4 exercises 1-5
# Regression diagnostics for sat, teengamb, prostate and swiss (questions a-f),
# and a check for serial correlation in divusa.
# The same six checks are needed four times, so they are written once as a
# function, check_model(), and then called for each dataset.
# =============================================================================

library(faraway)

# -----------------------------------------------------------------------------
# check_model(g, data, labs): runs checks (a)-(f) on a fitted lm() model g.
#   data = the data frame g was fitted to; labs = labels for the cases (e.g. state or country names) used in the plots.
# It prints the key numbers and draws the plots.
# -----------------------------------------------------------------------------
check_model <- function(g, data, labs, title = "") {
  n <- length(residuals(g))
  p <- length(coef(g))

  # (a) Constant variance: residuals vs fitted, plus the |residual| regression
  par(mfrow = c(1, 2))
  plot(fitted(g), residuals(g), xlab = "Fitted", ylab = "Residuals",
       main = paste(title, "(a) residuals vs fitted"))
  abline(h = 0)
  # (b) Normality: Q-Q plot of the residuals
  qqnorm(residuals(g), main = paste(title, "(b) normal Q-Q"))
  qqline(residuals(g))
  par(mfrow = c(1, 1))
  cat("\n(a) |residual| ~ fitted slope test:\n")
  print(summary(lm(abs(residuals(g)) ~ fitted(g)))$coefficients[2, ])
  cat("\n(b) Shapiro-Wilk:\n")
  print(shapiro.test(residuals(g)))

  # (c) Leverage: hatvalues() gives the leverages; flag those above 2p/n
  h <- hatvalues(g)
  names(h) <- labs
  cat("\n(c) Leverages above 2p/n =", round(2 * p / n, 3), ":\n")
  print(round(sort(h[h > 2 * p / n], decreasing = TRUE), 3))

  # (d) Outliers: jackknife residuals against the Bonferroni critical value
  jack <- rstudent(g)
  names(jack) <- labs
  crit <- abs(qt(0.05 / (2 * n), n - p - 1))
  cat("\n(d) Largest |jackknife residual|:", round(max(abs(jack)), 3),
      "at", labs[which.max(abs(jack))], "; Bonferroni critical value:", round(crit, 3), "\n")
  # Bonferroni p-value of the largest: n times its ordinary two-sided p-value (capped at 1)
  cat("    Bonferroni p-value:", signif(min(1, n * 2 * pt(-max(abs(jack)), n - p - 1)), 3), "\n")

  # (e) Influence: Cook's distances, and the change in the fit when the largest is dropped
  cook <- cooks.distance(g)
  names(cook) <- labs
  cat("\n(e) Three largest Cook's distances:\n")
  print(round(sort(cook, decreasing = TRUE)[1:3], 3))
  par(mfrow = c(1, 2))
  halfnorm(hatvalues(g), 3, labs = labs, ylab = "Leverages", main = paste(title, "(c)"))
  halfnorm(cook, 3, labs = labs, ylab = "Cook's distances", main = paste(title, "(e)"))
  par(mfrow = c(1, 1))
  # Refit the same formula without the most influential case.
  # formula(g): the model formula; data[cook < max(cook), ] keeps every other row
  g_drop <- lm(formula(g), data = data[cook < max(cook), ])
  cat("    Coefficients with and without", labs[which.max(cook)], ":\n")
  print(round(cbind(all = coef(g), dropped = coef(g_drop)), 4))

  # (f) Structure: termplot() draws a partial residual plot for every predictor;
  # smooth = panel.smooth adds a smooth curve to show any bend
  k <- p - 1                               # number of predictors
  par(mfrow = c(ceiling(k / 3), min(k, 3)))
  termplot(g, partial.resid = TRUE, smooth = panel.smooth, ask = FALSE)
  par(mfrow = c(1, 1))

  invisible(NULL)   # invisible(): return nothing without printing
}


# =============================================================================
# Exercise 1: sat - total SAT on expend, salary, ratio, takers
# =============================================================================
data(sat)
g_sat <- lm(total ~ expend + salary + ratio + takers, data = sat)
print(summary(g_sat))
check_model(g_sat, data = sat, labs = row.names(sat), title = "sat")

# What stands out (sat):
# [numbers to be filled in from your run_all.R output]
# - (a) constant variance, (b) normality: read off the slope test and Shapiro-Wilk.
# - (c) high-leverage states, (d) outlier test, (e) the state whose removal moves
#   the coefficients most.
# - (f) The takers partial residual plot is expected to be curved: the score
#   drops steeply as takers goes from about 5% to 20% and then levels off,
#   which suggests a transformation of takers (e.g. log) - Chapter 7.


# =============================================================================
# Exercise 2: teengamb - gamble on sex, status, income, verbal
# =============================================================================
data(teengamb)
g_tg <- lm(gamble ~ ., data = teengamb)
print(summary(g_tg))
check_model(g_tg, data = teengamb, labs = as.character(1:nrow(teengamb)), title = "teengamb")

# What stands out (teengamb):
# [numbers to be filled in from your run_all.R output]
# - Chapter 1 found the spread of gamble grows with income; expect a funnel in
#   the residuals vs fitted plot and a significant |residual| slope.
# - The large positive residual found in Chapter 2 exercise 1 is the case to
#   watch in the outlier and Cook's distance checks.


# =============================================================================
# Exercise 3: prostate - lpsa on all other variables
# =============================================================================
data(prostate)
g_pr <- lm(lpsa ~ ., data = prostate)
print(summary(g_pr))
check_model(g_pr, data = prostate, labs = as.character(1:nrow(prostate)), title = "prostate")

# What stands out (prostate):
# [numbers to be filled in from your run_all.R output]
# - Chapter 1 found case 32 has an extreme lweight (6.11 vs a next-largest of
#   4.78); expect it to show up with high leverage.


# =============================================================================
# Exercise 4: swiss - Fertility in 47 French-speaking Swiss provinces, 1888
#   Agriculture, Examination, Education, Catholic (%), Infant.Mortality
#   (swiss comes with R itself)
# =============================================================================
data(swiss)
g_sw <- lm(Fertility ~ ., data = swiss)
print(summary(g_sw))
check_model(g_sw, data = swiss, labs = row.names(swiss), title = "swiss")

# What stands out (swiss):
# [numbers to be filled in from your run_all.R output]
# - Catholic is bimodal (provinces are mostly Catholic or mostly Protestant),
#   so its partial residual plot shows two clusters rather than a spread.


# =============================================================================
# Exercise 5: divusa - divorce on unemployed, femlab, marriage, birth, military
#   (year left out of the predictors); check for serial correlation
# =============================================================================
data(divusa)
g_dv <- lm(divorce ~ unemployed + femlab + marriage + birth + military, data = divusa)
print(summary(g_dv))
r <- residuals(g_dv)
nr <- length(r)

par(mfrow = c(1, 2))
# Residuals in time order: long runs above or below 0 mean positive correlation
plot(divusa$year, r, type = "l", xlab = "year", ylab = "Residuals",
     main = "Residuals over time")
abline(h = 0)
# Each residual against the next one
plot(r[-nr], r[-1], xlab = expression(hat(epsilon)[i]),
     ylab = expression(hat(epsilon)[i + 1]), main = "Successive residuals")
par(mfrow = c(1, 1))

# Lag-1 correlation and regression of each residual on the previous one
print(cor(r[-nr], r[-1]))
print(summary(lm(r[-1] ~ -1 + r[-nr]))$coefficients)
# Durbin-Watson statistic by hand (about 2 = no correlation, near 0 = strong positive)
print(sum(diff(r)^2) / sum(r^2))
if (requireNamespace("lmtest", quietly = TRUE)) {
  print(lmtest::dwtest(g_dv))
}

# What stands out (divusa):
# [numbers to be filled in from your run_all.R output]
# - Expect strong positive serial correlation: long runs of same-sign
#   residuals, a clear positive slope in the successive-residual plot, and a
#   Durbin-Watson statistic far below 2. Divorce rates change slowly from year
#   to year, so this year's error carries into next year. The usual standard
#   errors and p-values are then too optimistic (see GLS, Chapter 6).
