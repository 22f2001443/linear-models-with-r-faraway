# Checking the structure of the model: added variable and partial residual plots (savings data)

library(faraway)
data(savings)
g <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)

par(mfrow = c(1, 2))

# ---- Added variable (partial regression) plot for pop15 ----
# d: sr with the effect of the other predictors removed
d <- residuals(lm(sr ~ pop75 + dpi + ddpi, savings))
# m: pop15 with the effect of the other predictors removed
m <- residuals(lm(pop15 ~ pop75 + dpi + ddpi, savings))
plot(m, d, xlab = "pop15 residuals", ylab = "Savings residuals",
     main = "Added variable plot")
# The slope of d on m equals the pop15 coefficient in the full model
print(coef(lm(d ~ m)))
print(coef(g))
abline(0, coef(g)["pop15"])

# ---- Partial residual plot for pop15 ----
# residuals + beta_pop15 * pop15 against pop15
plot(savings$pop15, residuals(g) + coef(g)["pop15"] * savings$pop15,
     xlab = "pop'n under 15", ylab = "Savings (adjusted)", main = "Partial residual plot")
abline(0, coef(g)["pop15"])
par(mfrow = c(1, 1))

# prplot(): (faraway) the same partial residual plot in one line; 1 = first predictor
prplot(g, 1)

# The plot shows two groups, so fit the model separately in each
g1 <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings, subset = (pop15 > 35))
g2 <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings, subset = (pop15 < 35))
print(summary(g1))
print(summary(g2))

# What it shows (numbers from the book):
# - The slope in the added variable plot is -0.4612, exactly the pop15
#   coefficient in the full model. The plot shows the pop15 relationship
#   after adjusting for the other predictors.
# - Added variable plots are better for spotting outliers and influential
#   points; partial residual plots are better for spotting non-linearity.
# - The partial residual plot shows two distinct groups of countries.
# - In the high pop15 group (> 35%, mostly less developed countries) there is
#   no relationship at all: F p-value = 0.52.
# - In the low pop15 group (developed countries) the relationship is strong
#   (p = 0.0027, R^2 = 0.51), mainly through ddpi (p = 0.0067), with pop15
#   borderline (p = 0.061). A purely numerical analysis could easily miss this.
