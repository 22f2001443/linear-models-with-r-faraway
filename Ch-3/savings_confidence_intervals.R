# Confidence intervals for the coefficients and a joint confidence region, using the savings data

library(faraway)
data(savings)
g <- lm(sr ~ ., savings)   # sr ~ . : all other columns are predictors
print(summary(g))

# qt(0.975, 45): the t critical value for a 95% interval with 45 df
tcrit <- qt(0.975, df.residual(g))
print(tcrit)

# Coefficient table: column 1 = estimate, column 2 = standard error
ct <- summary(g)$coefficients
# 95% CI for beta_pop75: estimate +- t * SE
print(ct["pop75", 1] + c(-1, 1) * tcrit * ct["pop75", 2])
# 95% CI for beta_ddpi
print(ct["ddpi", 1] + c(-1, 1) * tcrit * ct["ddpi", 2])

# confint(): all the one-at-a-time 95% intervals in one go
print(confint(g))


# ---- Joint 95% confidence region for (beta_pop15, beta_pop75) ----
# The region is the ellipse of beta values with
#   (b - beta)' V^-1 (b - beta) <= 2 * F(0.95; 2, n - p),
# where b holds the two estimates and V is their 2 x 2 covariance matrix.
# (The book uses the ellipse package: plot(ellipse(g, c(2, 3)), type = "l").
#  The base-R version below draws the same curve without installing anything.)
b <- coef(g)[c("pop15", "pop75")]
# vcov(): estimated covariance matrix of the coefficients; [rows, cols] keeps the 2 x 2 block
V <- vcov(g)[c("pop15", "pop75"), c("pop15", "pop75")]
radius <- sqrt(2 * qf(0.95, 2, df.residual(g)))   # qf(): F quantile
# seq(): 200 angles round a circle; cbind(): joins cos and sin as two columns
theta <- seq(0, 2 * pi, length.out = 200)
circle <- cbind(cos(theta), sin(theta))
# chol(V): a matrix R with t(R) %*% R = V. Multiplying the unit circle by it
# stretches the circle into the right ellipse; sweep(..., "+") then shifts it to the centre b
ell <- sweep(radius * circle %*% chol(V), 2, b, "+")

plot(ell, type = "l", xlim = c(-1, 0), xlab = "pop15 coefficient",
     ylab = "pop75 coefficient", main = "95% joint confidence region")
points(0, 0)                          # the origin
points(b[1], b[2], pch = 18)          # the estimates
# abline(v = ...) / abline(h = ...): vertical / horizontal lines at the one-at-a-time CIs
abline(v = confint(g)["pop15", ], lty = 2)
abline(h = confint(g)["pop75", ], lty = 2)

# The origin is outside the ellipse exactly when the F-test of
# H0: beta_pop15 = beta_pop75 = 0 rejects at 5%. Check it numerically:
# solve(V): the inverse of V; drop(): turns the 1 x 1 matrix into a plain number
q0 <- drop(t(b) %*% solve(V) %*% b) / 2
print(c(F = q0, critical = qf(0.95, 2, df.residual(g))))
print(anova(lm(sr ~ dpi + ddpi, savings), g))   # same F-test via anova()

# Correlation between the predictors vs correlation between their coefficients
print(cor(savings$pop15, savings$pop75))
# summary(..., correlation = TRUE)$correlation: correlation matrix of the estimates
print(round(summary(g, correlation = TRUE)$correlation, 3))

# What it shows (numbers from the book):
# - 95% CI for pop75: (-3.87, 0.49). It contains 0, matching p = 0.125 > 0.05.
# - 95% CI for ddpi: (0.015, 0.805). It excludes 0, so H0: beta_ddpi = 0 is
#   rejected at 5% (p = 0.042). But it is wide: the upper end is about 50 times
#   the lower, so the size of the growth effect is very uncertain.
# - The origin lies outside the ellipse, so the joint hypothesis
#   beta_pop15 = beta_pop75 = 0 is rejected, even though the pop75 interval
#   alone includes zero. Prefer the joint region when coefficients are correlated.
# - The dashed lines are not tangent to the ellipse: the rectangle of the two
#   separate CIs is not a valid joint 95% region.
# - pop15 and pop75 are strongly negatively correlated (r = -0.908), but their
#   coefficient estimates are positively correlated (0.765). The sign often
#   flips like this.
