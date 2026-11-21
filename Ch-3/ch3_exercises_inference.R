# =============================================================================
# Linear Models with R (Faraway), Chapter 3 exercises 1-5
# Confidence intervals and regions, predictions, permutation tests and F-tests
# on prostate, teengamb and sat.
# =============================================================================

library(faraway)


# =============================================================================
# Exercise 1: prostate - lpsa on all other variables
# =============================================================================
data(prostate)
g <- lm(lpsa ~ ., data = prostate)   # . = every other column
print(summary(g))

# (a) 90% and 95% CIs for the age coefficient
# confint(model, parm, level): CI for the named coefficient
print(confint(g, "age", level = 0.90))
print(confint(g, "age", level = 0.95))
print(summary(g)$coefficients["age", ])   # the p-value, to check the deduction

# (b) 95% joint confidence region for (age, lbph), with the origin marked
b <- coef(g)[c("age", "lbph")]
V <- vcov(g)[c("age", "lbph"), c("age", "lbph")]          # their 2 x 2 covariance matrix
radius <- sqrt(2 * qf(0.95, 2, df.residual(g)))
theta <- seq(0, 2 * pi, length.out = 200)
# Stretch a unit circle by chol(V) and shift it to b (see savings_confidence_intervals.R)
ell <- sweep(radius * cbind(cos(theta), sin(theta)) %*% chol(V), 2, b, "+")
# range(): smallest and largest value; including 0 makes sure the origin is on the plot
plot(ell, type = "l", xlim = range(ell[, 1], 0), ylim = range(ell[, 2], 0),
     xlab = "age coefficient", ylab = "lbph coefficient",
     main = "95% joint confidence region")
points(0, 0, pch = 3, cex = 1.5)    # origin
points(b[1], b[2], pch = 18)        # estimates
abline(v = confint(g)["age", ], lty = 2)
abline(h = confint(g)["lbph", ], lty = 2)
# The test the picture answers: H0: beta_age = beta_lbph = 0 (with the other predictors in)
print(anova(lm(lpsa ~ . - age - lbph, data = prostate), g))

# (c) Prediction for a new patient
x0 <- data.frame(lcavol = 1.44692, lweight = 3.62301, age = 65, lbph = 0.30010,
                 svi = 0, lcp = -0.79851, gleason = 7, pgg45 = 15)
# For one new patient the "prediction" interval is the right one;
# the "confidence" interval is for the mean lpsa of all such patients
print(predict(g, x0, interval = "prediction"))
print(predict(g, x0, interval = "confidence"))

# (d) Same patient but aged 20
x20 <- x0
x20$age <- 20
print(predict(g, x20, interval = "prediction"))
print(range(prostate$age))   # the ages actually in the data

# (e) Permutation test for age, matching its t-test
tobs <- summary(g)$coefficients["age", "t value"]
set.seed(123)
nperm <- 4000
tstats <- numeric(nperm)
pp <- prostate                    # a copy we can shuffle
for (i in 1:nperm) {
  pp$age <- sample(prostate$age)  # shuffle age, keep everything else fixed
  tstats[i] <- summary(lm(lpsa ~ ., data = pp))$coefficients["age", "t value"]
}
print(c(observed_t = tobs,
        perm_p = mean(abs(tstats) >= abs(tobs)),
        normal_p = summary(g)$coefficients["age", "Pr(>|t|)"]))

# What stands out (prostate full model):
# [numbers to be filled in from your run_all.R output]
# - (a) If 0 is inside the 95% interval but outside the 90% interval, the
#   age p-value must lie between 0.05 and 0.10 (the printed p-value confirms).
# - (b) The ellipse tests H0: beta_age = beta_lbph = 0 jointly. Origin outside
#   the ellipse = reject at 5%; inside = do not reject. The anova() F-test
#   gives the same verdict.
# - (c) Use the prediction interval: it covers one patient's lpsa, so it
#   includes the patient-to-patient error and is much wider than the CI for
#   the mean.
# - (d) The age-20 interval is wider because age 20 is far outside the observed
#   ages (41-79): x0'(X'X)^-1 x0 grows with the distance from the centre of the
#   data. It is an extrapolation and should not be trusted.
# - (e) The permutation p-value should be close to the normal-theory p-value.


# =============================================================================
# Exercise 2: drop the predictors not significant at 5% and compare
# =============================================================================
coefs <- summary(g)$coefficients
# Rows with p < 0.05, without the intercept
keep <- setdiff(rownames(coefs)[coefs[, "Pr(>|t|)"] < 0.05], "(Intercept)")
print(keep)
g2 <- lm(reformulate(keep, response = "lpsa"), data = prostate)
print(summary(g2))

# (a) Predictions for both patients from the smaller model
print(predict(g2, x0, interval = "prediction"))
print(predict(g2, x20, interval = "prediction"))
# Width of each interval, full vs reduced model
# function(p): a small helper that takes a predict() result and returns upper - lower
width <- function(p) unname(p[, "upr"] - p[, "lwr"])
print(c(full_65 = width(predict(g, x0, interval = "prediction")),
        full_20 = width(predict(g, x20, interval = "prediction")),
        small_65 = width(predict(g2, x0, interval = "prediction")),
        small_20 = width(predict(g2, x20, interval = "prediction"))))

# (b) F-test of the small model against the full one
print(anova(g2, g))

# What stands out (reduced prostate model):
# [numbers to be filled in from your run_all.R output]
# - If age is not among the kept predictors, the two patients get the same
#   prediction, and the age-20 interval is no longer inflated.
# - (b) A large anova() p-value means the dropped variables add nothing
#   significant together, so the smaller model is preferred: it is simpler
#   and gives more stable predictions.


# =============================================================================
# Exercise 3: teengamb - gamble on sex, status, income, verbal
# =============================================================================
data(teengamb)
g <- lm(gamble ~ ., data = teengamb)
print(summary(g))

# (c) A male with average status, income and verbal, then one with maximal values
xm  <- data.frame(sex = 0, status = mean(teengamb$status),
                  income = mean(teengamb$income), verbal = mean(teengamb$verbal))
xmx <- data.frame(sex = 0, status = max(teengamb$status),
                  income = max(teengamb$income), verbal = max(teengamb$verbal))
print(predict(g, xm,  interval = "prediction"))
print(predict(g, xmx, interval = "prediction"))
print(predict(g, xm,  interval = "confidence"))
print(predict(g, xmx, interval = "confidence"))

# (d) income alone vs the full model
gi <- lm(gamble ~ income, data = teengamb)
print(anova(gi, g))

# What stands out (teengamb):
# [numbers to be filled in from your run_all.R output]
# - (a) The variables with p < 0.05 in the summary are the significant ones.
# - (b) sex = 1 for females, so its coefficient is how much less (if negative)
#   a female is predicted to gamble than a male with the same status, income
#   and verbal score.
# - (c) The interval for the "maximal" male is wider: his predictor values
#   are far from the averages, where x0'(X'X)^-1 x0 is smallest.
#   The average male sits near the centre of the data, so his interval is
#   the narrowest possible.
# - (d) A small anova() p-value means status, verbal and sex together add
#   significantly to income alone.


# =============================================================================
# Exercise 4: sat - SAT total on school spending variables
# =============================================================================
data(sat)

# (a) expend, ratio and salary
g1 <- lm(total ~ expend + ratio + salary, data = sat)
print(summary(g1))
# H0: beta_salary = 0 -> the salary row of the table
print(summary(g1)$coefficients["salary", ])
# H0: all three = 0 -> the overall F-test at the bottom of the summary
f1 <- summary(g1)$fstatistic
print(c(f1, p = unname(1 - pf(f1["value"], f1["numdf"], f1["dendf"]))))

# (b) add takers
g2 <- lm(total ~ expend + ratio + salary + takers, data = sat)
print(summary(g2))
print(summary(g2)$coefficients["salary", ])
print(anova(g1, g2))
# The F for adding one variable equals the square of its t value
t_takers <- summary(g2)$coefficients["takers", "t value"]
print(c(t_squared = t_takers^2, F = anova(g1, g2)$F[2]))

# What stands out (sat):
# [numbers to be filled in from your run_all.R output]
# - Without takers, the spending variables look weakly related (or wrongly
#   signed) to total. Chapter 1 showed why: takers is the dominant driver
#   (r = -0.89 with total) and is itself correlated with spending.
# - Adding takers changes the picture: takers is hugely significant and the
#   salary test changes.
# - The anova() F-statistic equals t^2 for takers and the two p-values are
#   identical: for a single added variable the F-test and t-test are the same test.


# =============================================================================
# Exercise 5: a formula relating R^2 and the overall F-statistic
#   RSS = (1 - R^2) TSS and TSS - RSS = R^2 TSS, so
#   F = ((TSS - RSS) / (p - 1)) / (RSS / (n - p)) = (R^2 / (p - 1)) / ((1 - R^2) / (n - p))
# Check it numerically on the sat model g2:
# =============================================================================
r2 <- summary(g2)$r.squared
n <- nrow(sat)
p <- length(coef(g2))    # number of parameters, including the intercept
print(c(from_R2 = (r2 / (p - 1)) / ((1 - r2) / (n - p)),
        from_summary = unname(summary(g2)$fstatistic["value"])))

# What stands out (R^2 and F):
# - The two numbers agree exactly. F is a monotone function of R^2 for fixed
#   n and p, so the overall F-test is really a test of whether R^2 is larger
#   than chance would give.
