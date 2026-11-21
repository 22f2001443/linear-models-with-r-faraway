# Testing for lack of fit using replicates, with the corrosion data
# 13 specimens of 90/10 Cu-Ni alloy: Fe = iron content (%), loss = weight loss from corrosion
# (mg per dm^2 per day) after 60 days in seawater. Some Fe values are repeated (replicates).

library(faraway)
data(corrosion)

plot(loss ~ Fe, corrosion, xlab = "Iron content", ylab = "Weight loss")
g <- lm(loss ~ Fe, corrosion)
print(summary(g))
abline(coef(g))

# factor(Fe): treats each distinct Fe value as its own group, so the model has one mean per
# group. Its residual variance is the "pure error": model-free, from the replicates only
ga <- lm(loss ~ factor(Fe), corrosion)
points(corrosion$Fe, fitted(ga), pch = 18)     # the group means
# Lack-of-fit test: straight line vs one-mean-per-group
print(anova(g, ga))
# Pure error SD vs the straight line's residual SE
print(c(pure_error_sd = summary(ga)$sigma, line_sd = summary(g)$sigma))
print(summary(ga)$r.squared)                   # the best R^2 any model could reach here

# Overfitting: a 6th-degree polynomial goes through all the group means
gp <- lm(loss ~ Fe + I(Fe^2) + I(Fe^3) + I(Fe^4) + I(Fe^5) + I(Fe^6), corrosion)
plot(loss ~ Fe, data = corrosion, ylim = c(60, 130), main = "6th-degree polynomial")
points(corrosion$Fe, fitted(ga), pch = 18)
grid <- seq(0, 2, len = 50)
# lines(): adds the fitted curve, predicted on a fine grid of Fe values
lines(grid, predict(gp, data.frame(Fe = grid)))
print(summary(gp)$r.squared)

# What it shows (numbers from the book):
# - The straight line looks excellent: R^2 = 0.97, residual SE 3.06.
# - Yet the lack-of-fit test rejects: F = 9.28 on 5 and 6 df, p = 0.0086. The
#   pure error SD is only about 1.4, much smaller than 3.06, so the line misses
#   the group means by more than the replicates' scatter allows.
# - Before blaming the model, check the replicates are genuine (independent
#   specimens, not repeat measurements of one), or whether a third variable is
#   involved.
# - Even the saturated one-mean-per-group model only reaches R^2 = 0.997, not
#   1: with replicates the maximum achievable R^2 is below 100%.
# - The 6th-degree polynomial has R^2 = 0.9965 but is absurd: it dips at
#   Fe = 1.7 and then shoots up. Better fit statistics are not better models.
