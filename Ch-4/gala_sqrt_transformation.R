# A variance-stabilising transformation, using the gala data
# Species is a count. For counts the variance tends to grow with the mean, so a square root helps.

library(faraway)
data(gala)

par(mfrow = c(1, 2))
# Untransformed response
gg <- lm(Species ~ Area + Elevation + Scruz + Nearest + Adjacent, gala)
plot(fitted(gg), residuals(gg), xlab = "Fitted", ylab = "Residuals",
     main = "Untransformed response")
abline(h = 0)

# sqrt(): square root of the response
gs <- lm(sqrt(Species) ~ Area + Elevation + Scruz + Nearest + Adjacent, gala)
plot(fitted(gs), residuals(gs), xlab = "Fitted", ylab = "Residuals",
     main = "Square root response")
abline(h = 0)
par(mfrow = c(1, 1))

# The same |residual| regression check as for the savings data, for both models
print(summary(lm(abs(residuals(gg)) ~ fitted(gg)))$coefficients)
print(summary(lm(abs(residuals(gs)) ~ fitted(gs)))$coefficients)

# What it shows:
# - On the original scale the residuals fan out as the fitted values grow:
#   clear non-constant variance.
# - After the square root the spread is roughly constant across the fitted values.
# - The square root suits counts because a Poisson count has variance equal
#   to its mean, and the square root makes that variance roughly constant.
