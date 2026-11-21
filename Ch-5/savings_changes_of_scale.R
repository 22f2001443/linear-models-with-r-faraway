# Changes of scale, using the savings data
# Rescaling a predictor changes its coefficient and SE but not the t-tests, F-test or R^2.

library(faraway)
data(savings)

g <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)
print(summary(g))

# Income in thousands of dollars instead of dollars. I() does the division inside the formula
g <- lm(sr ~ pop15 + pop75 + I(dpi / 1000) + ddpi, savings)
print(summary(g))

# Standardise every variable: scale() subtracts the mean and divides by the SD of each column.
# scale() returns a matrix, so data.frame() turns it back into a data frame
scsav <- data.frame(scale(savings))
g <- lm(sr ~ ., scsav)
print(summary(g))

# What it shows (numbers from the book):
# - Dividing dpi by 1000 multiplies its coefficient and SE by 1000
#   (-0.000337 -> -0.337, SE 0.000931 -> 0.931). The t value (-0.36), p-value,
#   residual SE (3.8), R^2 (0.338) and F (5.76) are unchanged.
# - With everything standardised, the intercept is 0 (the fitted plane passes
#   through the means, now the origin) and the coefficients are on a common
#   scale: pop15 -0.94, pop75 -0.49, dpi -0.07, ddpi 0.26. Each is the change in
#   sr, in SDs, for a 1 SD change in that predictor.
# - The t-tests and R^2 are again unchanged. Standardised coefficients make
#   predictors easier to compare, but are less natural to interpret.
