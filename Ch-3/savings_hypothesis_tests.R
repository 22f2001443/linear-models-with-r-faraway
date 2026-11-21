# F-tests and t-tests for comparing nested models, using the savings data (50 countries, 1960-70)
# Columns: sr (savings rate, % of disposable income), pop15 (% under 15), pop75 (% over 75),
# dpi (per-capita disposable income, US$), ddpi (% growth rate of dpi)

library(faraway)
data(savings)
print(head(savings))

# Full model with all four predictors
g <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)
print(summary(g))


# ---- Test of all the predictors: H0: all four slopes = 0 ----
# By hand: F = ((TSS - RSS) / (p - 1)) / (RSS / (n - p))
tss <- sum((savings$sr - mean(savings$sr))^2)   # total sum of squares
rss <- deviance(g)                               # residual sum of squares
print(c(TSS = tss, RSS = rss, df = df.residual(g)))   # df = n - p = 50 - 5 = 45
fstat <- ((tss - rss) / 4) / (rss / df.residual(g))
print(fstat)
# pf(): the F distribution function, so 1 - pf() is the upper-tail p-value
print(1 - pf(fstat, 4, df.residual(g)))


# ---- Testing one predictor: H0: beta_pop15 = 0 ----
# Fit the model without pop15 (the null model) and compare RSS
g2 <- lm(sr ~ pop75 + dpi + ddpi, savings)
rss2 <- deviance(g2)
print(rss2)
fstat <- (deviance(g2) - deviance(g)) / (deviance(g) / df.residual(g))
print(fstat)
print(1 - pf(fstat, 1, df.residual(g)))

# The t-statistic from the summary table is the square root of this F
print(sqrt(fstat))
# summary(g)$coefficients is the coefficient table; row "pop15", column "t value"
tstat <- summary(g)$coefficients["pop15", "t value"]
print(tstat)
# pt(): t distribution function; 2 * (1 - pt(|t|)) is the two-sided p-value
print(2 * (1 - pt(sqrt(fstat), 45)))

# anova(small, big): does the whole comparison in one line
print(anova(g2, g))


# ---- Testing a pair of predictors: H0: beta_pop75 = beta_ddpi = 0 ----
g3 <- lm(sr ~ pop15 + dpi, savings)
print(anova(g3, g))


# ---- Testing a subspace: H0: beta_pop15 = beta_pop75 ----
# The dot in sr ~ . means "all other columns of the data frame"
g <- lm(sr ~ ., savings)
# I(): makes R add pop15 + pop75 arithmetically, instead of reading + as "another term".
# Under H0 both age groups share one coefficient, so the model uses their sum
gr <- lm(sr ~ I(pop15 + pop75) + dpi + ddpi, savings)
print(anova(gr, g))


# ---- Testing a fixed value: H0: beta_ddpi = 0.5 ----
# offset(): a term with its coefficient fixed (here at 0.5) instead of estimated
gr <- lm(sr ~ pop15 + pop75 + dpi + offset(0.5 * ddpi), savings)
print(anova(gr, g))

# Simpler: t = (estimate - 0.5) / SE
est <- summary(g)$coefficients["ddpi", "Estimate"]
se  <- summary(g)$coefficients["ddpi", "Std. Error"]
tstat <- (est - 0.5) / se
print(tstat)
print(2 * pt(-abs(tstat), df.residual(g)))   # two-sided p-value
tstat^2                               # equals the F value from anova(gr, g)

# What it shows (numbers from the book):
# - Overall F = 5.756 on 4 and 45 df, p = 0.00079: the predictors together
#   clearly help explain sr. TSS = 983.63 and RSS = 650.71.
# - Dropping pop15: RSS rises to 797.72, F = 10.17, p = 0.0026, which is exactly
#   the t-test in the summary (t = -3.19, and 3.19^2 = 10.17). A t-test for one
#   coefficient is an F-test comparing the models with and without it.
# - pop75 and ddpi together: F = 3.23 on 2 and 45 df, p = 0.049, just
#   significant at 5%. Testing a pair jointly can give a different answer
#   from looking at the two separate t-tests.
# - beta_pop15 = beta_pop75: p = 0.21, so there is no evidence the young and
#   old proportions need separate coefficients, in this model.
# - beta_ddpi = 0.5: t = -0.46, p = 0.65, the same p-value as the offset
#   F-test, and t^2 = 0.21 = F. Not rejected.
# - Every test is relative to the other predictors in the model: change them
#   and the result can change.
