# Confidence intervals for predictions, using the gala data
# Two different questions:
#   - the mean number of species for islands like x0  -> confidence interval
#   - the number of species on one new island like x0 -> prediction interval (wider)

library(faraway)
data(gala)
g <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala)

# A new island: 1 for the intercept, then Area, Elevation, Nearest, Scruz, Adjacent
x0 <- c(1, 0.08, 93, 6.0, 12.0, 0.34)
# Point prediction = x0' beta-hat; sum(x0 * coef(g)) does the multiply-and-add
y0 <- sum(x0 * coef(g))
print(y0)

# t critical value with n - p = 24 df
tcrit <- qt(0.975, df.residual(g))
print(tcrit)
x <- model.matrix(g)
xtxi <- solve(t(x) %*% x)
sig <- summary(g)$sigma

# Mean response: half-width = t * sigma * sqrt(x0' (X'X)^-1 x0)
bm <- drop(sqrt(x0 %*% xtxi %*% x0) * tcrit * sig)
print(c(y0 - bm, y0 + bm))
# Single future island: add 1 inside the square root for the new error
bp <- drop(sqrt(1 + x0 %*% xtxi %*% x0) * tcrit * sig)
print(c(y0 - bp, y0 + bp))

# predict() does this directly. It needs a data frame with the same column names
x0df <- data.frame(Area = 0.08, Elevation = 93, Nearest = 6.0, Scruz = 12, Adjacent = 0.34)
# se.fit = TRUE also returns the standard error of the fitted mean.
# str(): shows the structure of the returned list compactly
str(predict(g, x0df, se.fit = TRUE))
print(predict(g, x0df, interval = "confidence"))
print(predict(g, x0df, interval = "prediction"))


# ---- Extrapolation: bands get wider away from the data ----
grid <- seq(0, 100, 1)   # Nearest from 0 to 100 km, other predictors held at x0
p <- predict(g, data.frame(Area = 0.08, Elevation = 93, Nearest = grid,
                           Scruz = 12, Adjacent = 0.34),
             se.fit = TRUE, interval = "confidence")
# matplot(): plots each column of a matrix as its own line (fit, lower, upper)
matplot(grid, p$fit, lty = c(1, 2, 2), type = "l", col = 1,
        xlab = "Nearest", ylab = "Species",
        main = "Predicted Species with 95% bands for the mean")
# rug(): small tick marks on the x-axis at the observed Nearest values
rug(gala$Nearest)

# What it shows (numbers from the book):
# - Predicted species for the new island: 33.92 (about 34).
# - 95% CI for the mean response: (1.03, 66.81).
# - 95% prediction interval for one new island: (-96.15, 163.99). It is much
#   wider because it includes the island-to-island variation (sigma = 61).
# - The negative lower limit is impossible for a count: the normal model is a
#   poor fit for counts. A transformed response (Chapter 7) or a Poisson model
#   would avoid this.
# - The bands widen as Nearest moves away from the observed values (the rug).
#   They only allow for uncertainty in the estimates, not for the model being
#   wrong outside the range of the data.
