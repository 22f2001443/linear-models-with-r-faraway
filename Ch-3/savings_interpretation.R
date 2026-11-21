# Interpreting coefficients in observational data, using the savings data
# The same predictor (pop75) is fitted in four models with different companions.

library(faraway)
data(savings)

# Four models, all containing pop75
g  <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)
g2 <- lm(sr ~ pop75 + dpi + ddpi, savings)
g3 <- lm(sr ~ pop75 + ddpi, savings)
g4 <- lm(sr ~ pop75, savings)

print(summary(g))
print(summary(g2))
print(summary(g3))
print(summary(g4))

# The pop75 row from each model side by side.
# lapply(): applies a function to each model in the list; rbind with do.call() stacks the rows
models <- list(four = g, three = g2, two = g3, one = g4)
pop75_rows <- do.call(rbind, lapply(models, function(m) summary(m)$coefficients["pop75", ]))
print(round(pop75_rows, 4))

# pop15 and pop75 measure the same thing (the age structure), so they are strongly correlated
print(cor(savings$pop15, savings$pop75))

# Predictions are far more stable than the coefficients
x0 <- data.frame(pop15 = 32, pop75 = 3, dpi = 700, ddpi = 3)
# sapply(): applies predict() to each model and returns a vector
print(sapply(models, predict, newdata = x0))

# What it shows (numbers from the book):
# - pop75 coefficient: -1.69 (p = 0.13) with all four predictors, +0.95
#   (p = 0.22) without pop15, +1.07 (p = 0.023) with only ddpi, and +1.10
#   (p = 0.025) alone. Both the sign and the significance depend on what else
#   is in the model.
# - The flip happens because pop75 and pop15 are strongly negatively
#   correlated (r = -0.91): once pop15 is in the model, pop75 has little extra
#   to explain.
# - Predictions at the same x0 barely move: 9.73, 9.91, 10.08, 10.45.
#   Prediction is more stable than estimation of individual effects.
# - In observational data a coefficient is not "the effect of changing x by
#   one". It compares individuals with different x who share the other
#   predictor values, and no causal claim follows without more evidence.
