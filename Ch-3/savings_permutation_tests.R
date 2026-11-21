# Permutation tests, using the savings data
# Idea: if sr has no relation to the predictors, shuffling sr among the countries
# should give F-statistics as large as the real one quite often.

library(faraway)
data(savings)

# A model with just pop75 and dpi (its p-value is not tiny, which makes it more interesting)
g <- lm(sr ~ pop75 + dpi, savings)
print(summary(g))

# $fstatistic holds the F value and its two degrees of freedom
gs <- summary(g)
print(gs$fstatistic)
fobs <- gs$fstatistic["value"]   # the observed F-statistic

set.seed(123)        # set.seed(): makes the random permutations reproducible
nperm <- 4000
# numeric(n): a vector of n zeros, to be filled in by the loop
fstats <- numeric(nperm)
# for loop: repeat the fit 4000 times
for (i in 1:nperm) {
  # sample(sr): a random reordering (permutation) of the sr values
  ge <- lm(sample(sr) ~ pop75 + dpi, data = savings)
  fstats[i] <- summary(ge)$fstatistic["value"]
}
# Proportion of permuted F-statistics at least as large as the observed one = permutation p-value
# (mean() of TRUE/FALSE values gives the proportion of TRUEs)
print(mean(fstats > fobs))
# Normal-theory p-value from the F distribution, for comparison
print(1 - pf(fobs, gs$fstatistic["numdf"], gs$fstatistic["dendf"]))

# hist(): the permutation distribution, with the observed F marked in red
hist(fstats, breaks = 50, main = "Permutation distribution of F", xlab = "F")
abline(v = fobs, col = "red", lwd = 2)


# ---- Permutation test for one predictor (pop75) ----
# Now we shuffle the predictor instead of the response
tobs <- summary(g)$coefficients["pop75", "t value"]
print(tobs)
tstats <- numeric(nperm)
for (i in 1:nperm) {
  ge <- lm(sr ~ sample(pop75) + dpi, savings)
  # row 2 of the coefficient table is the shuffled pop75; column 3 is the t value
  tstats[i] <- summary(ge)$coefficients[2, 3]
}
# Two-sided: count permuted |t| larger than the observed |t|
print(mean(abs(tstats) > abs(tobs)))

# What it shows (numbers from the book; yours will vary slightly with the random draws):
# - Observed F = 2.68 on 2 and 47 df, normal-theory p = 0.079.
# - The permutation p-value was 0.074 in the book, very close to 0.079.
#   The permutation test does not assume normal errors, so the close agreement
#   is reassuring for the normal-theory result.
# - For pop75 alone: observed t = 1.68; about 10.5% of the permuted |t| values
#   are larger, matching the normal-based p-value of 0.10.
# - Rerunning with a different seed changes the third decimal; more
#   permutations would make the estimate more precise.
