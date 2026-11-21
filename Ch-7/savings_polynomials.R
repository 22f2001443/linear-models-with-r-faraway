# Polynomial regression, using savings: sr on ddpi (growth rate of income)

library(faraway)
data(savings)

# Add terms one at a time until the newest one is not significant
print(summary(lm(sr ~ ddpi, savings)))
print(summary(lm(sr ~ ddpi + I(ddpi^2), savings)))               # I(ddpi^2): ddpi squared
print(summary(lm(sr ~ ddpi + I(ddpi^2) + I(ddpi^3), savings)))

# Why lower-order terms stay in: shift ddpi by 10 and refit the quadratic
savings2 <- data.frame(savings, mddpi = savings$ddpi - 10)
print(summary(lm(sr ~ mddpi + I(mddpi^2), savings2)))

# poly(ddpi, 4): orthogonal polynomials up to degree 4. Each term is uncorrelated with the
# others, so one fit shows which degrees matter
g <- lm(sr ~ poly(ddpi, 4), savings)
print(summary(g))

# A second-degree polynomial in two predictors (a response surface):
# polym(): all terms pop15, ddpi, pop15^2, ddpi^2 and pop15*ddpi, as orthogonal polynomials
g2 <- lm(sr ~ polym(pop15, ddpi, degree = 2), savings)
print(summary(g2))

# Plot the chosen quadratic
gq <- lm(sr ~ ddpi + I(ddpi^2), savings)
plot(sr ~ ddpi, savings, main = "Quadratic fit")
grid <- seq(min(savings$ddpi), max(savings$ddpi), length.out = 100)
lines(grid, predict(gq, data.frame(ddpi = grid)))

# What it shows (numbers from the book):
# - Linear: ddpi significant (p = 0.031). Quadratic: the squared term is
#   significant (p = 0.013). Cubic: the cubed term is not (p = 0.99), so stop
#   at the quadratic. Starting from degree 4 and working down gives the same answer.
# - Shifting ddpi by 10 leaves the quadratic term unchanged (-0.093, p = 0.013)
#   but makes the linear term insignificant (p = 0.74). An arbitrary shift of
#   origin should not decide the model, so lower-order terms are kept
#   whenever a higher-order term is.
# - With poly(ddpi, 4) one summary is enough: degrees 1 and 2 are
#   significant (p = 0.025, 0.015), degrees 3 and 4 are not.
