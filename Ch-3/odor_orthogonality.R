# Designed experiments and orthogonality, using the odor data (15 runs of a chemical process)
# Columns: odor (odor score), temp, gas (gas/liquid ratio), pack (packing height),
# each predictor coded -1, 0, 1 (e.g. temp = (Fahrenheit - 80) / 40). A central composite design.

library(faraway)
data(odor)
print(odor)

# as.matrix(cbind(1, odor[, -1])): the X matrix, a column of 1s plus the three predictors
# (odor[, -1] drops column 1, the response)
x <- as.matrix(cbind(1, odor[, -1]))
# X'X is diagonal: the columns of X are orthogonal
print(t(x) %*% x)

g <- lm(odor ~ temp + gas + pack, odor)
# summary(..., correlation = TRUE): also prints the correlations between the estimates
print(summary(g, correlation = TRUE))

# Drop temp: with an orthogonal design the other estimates do not change at all
g <- lm(odor ~ gas + pack, odor)
print(summary(g))


# ---- Compare with observational data (savings), where X is not orthogonal ----
data(savings)
g <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)
print(summary(g))
# update(g, . ~ . - pop15): refits the same model with pop15 removed
g <- update(g, . ~ . - pop15)
print(summary(g))

# What it shows (numbers from the book):
# - X'X = diag(15, 8, 8, 8), so the correlations between the coefficient
#   estimates are all exactly 0, and the three slopes have equal standard
#   errors (12.7) because the design is balanced.
# - Dropping temp leaves the gas (-17.0) and pack (-21.4) estimates exactly
#   the same. Only sigma changes a little (36.0 -> 35.9), which nudges the
#   SEs and p-values. The effect of each factor can be read off on its own.
# - In the observational savings data, dropping pop15 flips the sign of the
#   pop75 coefficient from -1.69 to +0.95. Neither is significant here, but it
#   shows that in non-orthogonal data a coefficient's meaning depends on which
#   other predictors are in the model.
