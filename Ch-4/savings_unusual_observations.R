# Leverage, outliers and influential points, using the savings data

library(faraway)
data(savings)
g <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)


# ---- Leverage ----
# influence(): leave-one-out quantities for every case: $hat (leverages),
# $coefficients (change in each coefficient when the case is dropped), $sigma
ginf <- influence(g)
print(head(ginf$hat))
print(sum(ginf$hat))      # the leverages always add up to p, the number of parameters (5)

# row.names(): the country names, used to label the plots
countries <- row.names(savings)
# Rule of thumb: leverages above 2p/n deserve a closer look
p <- length(coef(g))
n <- nrow(savings)
print(2 * p / n)
# sort(..., decreasing = TRUE): largest first; [ ... > 2p/n] keeps only the large ones
print(sort(ginf$hat[ginf$hat > 2 * p / n], decreasing = TRUE))

par(mfrow = c(1, 2))
# halfnorm(): (faraway) half-normal plot; unusually large values stand out at the top right.
# The 2 largest are labelled with their country names
halfnorm(ginf$hat, labs = countries, ylab = "Leverages")

# Studentized residuals: residual / (sigma * sqrt(1 - h)), so they all have variance about 1
gs <- summary(g)
print(gs$sigma)
stud <- residuals(g) / (gs$sigma * sqrt(1 - ginf$hat))
# rstandard() computes the same thing; all.equal() checks the two agree
print(all.equal(stud, rstandard(g)))
qqnorm(stud)
abline(0, 1)       # the y = x line: standardized residuals should lie close to it
par(mfrow = c(1, 1))


# ---- Outliers ----
# rstudent(): jackknife (externally studentized) residuals, each case left out of its own fit.
# Under the model each one follows a t distribution with n - p - 1 df
jack <- rstudent(g)
print(jack[which.max(abs(jack))])     # the largest in absolute value
# Bonferroni critical value: level 0.05 shared out over n two-sided tests
print(qt(0.05 / (n * 2), n - p - 1))


# ---- Influential points ----
# cooks.distance(): Cook's distance, combining residual size and leverage
cook <- cooks.distance(g)
par(mfrow = c(1, 2))
halfnorm(cook, 3, labs = countries, ylab = "Cook's distances")   # label the 3 largest

# Refit without the largest (Libya): subset = (cook < max(cook)) keeps every other country
gl <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings, subset = (cook < max(cook)))
print(summary(gl))
print(summary(g))

# Change in the pop15 coefficient when each country is left out
plot(ginf$coefficients[, 2], ylab = "Change in pop15 coef")
# The book uses identify() to click on points. In a script we label the largest
# instead: which.max(abs()) finds it and text() writes its name next to the point
big <- which.max(abs(ginf$coefficients[, 2]))
text(big, ginf$coefficients[big, 2], countries[big], pos = 4)
par(mfrow = c(1, 1))

# Refit without Japan
gj <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings, subset = (countries != "Japan"))
print(summary(gj))

# What it shows (numbers from the book):
# - The leverages sum to 5 = p. The two largest leverages (labelled in the
#   half-normal plot) are Libya and the United States, the countries with
#   unusual predictor values.
# - Studentized residuals follow the y = x line in the Q-Q plot: no problem.
# - The largest jackknife residual is Zambia, 2.85. The Bonferroni critical
#   value is 3.53 in absolute value, so Zambia is NOT an outlier. Unless a
#   jackknife residual exceeds about 3.5, an outlier test is not worth doing.
# - Libya has the largest Cook's distance. Dropping it changes the ddpi
#   coefficient by about 50% (0.41 -> 0.61), too much for one country.
# - Japan stands out for the change in the pop15 coefficient. Without Japan,
#   ddpi is no longer significant (p = 0.099) and R^2 drops from 0.338 to 0.277.
