# Identifiability, using the gala data
# If one predictor is an exact linear combination of others, X'X cannot be inverted
# and the coefficients are not identifiable. Being "almost" a combination is nearly as bad.

library(faraway)
data(gala)

# Create a new variable: island area minus the area of the adjacent island.
# It is an exact linear combination of Area and Adjacent
gala$Adiff <- gala$Area - gala$Adjacent

# lm() with the redundant variable added
g <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent + Adiff, gala)
# summary() reports "1 not defined because of singularities": Adiff gets NA.
# R drops variables in reverse order of the formula, so the last one (Adiff) goes
print(summary(g))

# eigen(): eigenvalues and eigenvectors of X'X. An eigenvalue of (almost) zero shows a
# linear dependence; its eigenvector shows which columns are involved
X <- model.matrix(g)
e <- eigen(crossprod(X))
print(signif(e$values, 3))
# The last eigenvector belongs to the smallest eigenvalue. Its non-zero entries
# (Area, Adjacent and Adiff, with weights about +-0.58) show that
# Area - Adjacent - Adiff = 0
print(round(e$vectors[, ncol(X)], 3))


# ---- Near-identifiability ----

# set.seed(): fixes the random numbers so the result is the same every time you run it
set.seed(1)
# runif(30): 30 random numbers from Uniform(0, 1). 0.001 * (runif(30) - 0.5) is a tiny
# random shift between -0.0005 and +0.0005 added to Adiff
gala$Adiffe <- gala$Adiff + 0.001 * (runif(30) - 0.5)

# Now X'X can be inverted, but only just
g2 <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent + Adiffe, gala)
print(summary(g2))

# What it shows (numbers from the book; the near-identifiable fit is random):
# - With Adiff added, the rank of X is 6 but it has 7 columns, so R fits the
#   largest identifiable model: the output is identical to the original model
#   and Adiff is NA.
# - The eigenvector of the zero eigenvalue points to Area, Adjacent and Adiff:
#   the eigen-decomposition finds the culprit when the cause is not obvious.
# - After the tiny perturbation every coefficient is estimated, but Area,
#   Adjacent and Adiffe get estimates around +-2.4e+04 with standard errors of
#   about 4.7e+04 (book's values; yours will differ with the random numbers).
#   The fit itself (R^2 about 0.77) is unchanged: the data cannot tell these
#   three coefficients apart, so they are estimated very unstably.
#   Section 5.3 (collinearity) shows how to detect this.
