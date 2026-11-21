# Weighted least squares, using the fpe data: 1981 French presidential election,
# every fourth department. All counts in thousands.
#   EI = registered voters; A (Mitterrand), B (Giscard), C-K other candidates: first-round votes
#   A2, B2: Mitterrand's and Giscard's second-round votes
#   N = second-round total minus first-round total (extra voters in round 2)

library(faraway)
data(fpe)
print(head(fpe))

# If your version of the data has no N column, compute it as described in the book
if (is.null(fpe$N)) {
  fpe$N <- (fpe$A2 + fpe$B2) - rowSums(fpe[, c("A", "B", "C", "D", "E", "F", "G", "H", "J", "K")])
}

# Model: A2 = beta_A A + beta_B B + ... + beta_N N, where beta_i is the share of
# candidate i's voters who moved to Mitterrand. No intercept: "- 1".
# The error is like a sum over voters, so its variance grows with EI: weights = 1/EI
g <- lm(A2 ~ A + B + C + D + E + F + G + H + J + K + N - 1, fpe, weights = 1 / EI)
print(coef(g))

# Without weights, some of the small candidates' coefficients change a lot
print(coef(lm(A2 ~ A + B + C + D + E + F + G + H + J + K + N - 1, fpe)))

# Only the relative weights matter: multiplying them all by 53 changes nothing
print(coef(lm(A2 ~ A + B + C + D + E + F + G + H + J + K + N - 1, fpe, weights = 53 / EI)))

# Proportions should lie between 0 and 1. Crude fix: fix the coefficients that are above 1
# (A, G, K) at exactly 1 with offset(), and drop the negative ones (B, H) - i.e. set them to 0.
# Note: the book's model also leaves out J (setting its share to 0), even though J's
# estimate (0.61) was inside [0, 1]. This line copies the book; add + J to keep it.
g2 <- lm(A2 ~ offset(A + G + K) + C + D + E + F + N - 1, fpe, weights = 1 / EI)
print(coef(g2))

# What it shows (numbers from the book):
# - Weighted estimates: A 1.07, B -0.11, C 0.25, D 0.93, E 0.25, F 0.76,
#   G 1.97, H -0.57, J 0.61, K 1.21, N 0.53.
# - Unweighted, E jumps from 0.25 to 0.67 and K from 1.21 to 0.52: the weights
#   matter for the small candidates. Rescaling all weights by 53 changes nothing.
# - After the truncation fix: almost all Communist voters (D, 0.97) went to
#   Mitterrand. Surprisingly, about 23% of Gaullist voters (C) did too, which
#   appears to have been decisive. Ecology voters (E, 0.39) and new round-two
#   voters (N, 0.61) split roughly evenly.
