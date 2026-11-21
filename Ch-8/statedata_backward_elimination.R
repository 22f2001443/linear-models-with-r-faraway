# Testing-based variable selection (backward elimination), using the 1970s US state data
# state.x77 (comes with R): Population, Income, Illiteracy, Life.Exp, Murder, HS.Grad,
# Frost (days below freezing), Area. Response: Life.Exp (life expectancy, years)

data(state)
# data.frame(state.x77, row.names = state.abb): one row per state, labelled by its abbreviation
statedata <- data.frame(state.x77, row.names = state.abb)
g <- lm(Life.Exp ~ ., data = statedata)
print(summary(g))

# Backward elimination: remove the predictor with the largest p-value above 0.05, refit, repeat.
# update(g, . ~ . - Area): the same model without Area
g <- update(g, . ~ . - Area)
print(summary(g))
g <- update(g, . ~ . - Illiteracy)
print(summary(g))
g <- update(g, . ~ . - Income)
print(summary(g))
g <- update(g, . ~ . - Population)
print(summary(g))

# The same procedure written as a loop, so it works for any model.
# drop1(g, test = "F"): the F-test for dropping each term in turn
backward <- function(g, data, alpha = 0.05) {
  repeat {
    d <- drop1(g, test = "F")
    pv <- d[-1, "Pr(>F)"]                # [-1] skips the "<none>" row
    if (length(pv) == 0 || max(pv) <= alpha) return(g)
    worst <- rownames(d)[-1][which.max(pv)]
    cat("dropping", worst, "(p =", signif(max(pv), 3), ")\n")
    # update(formula(g), ...): the current formula with that term removed; refit on data
    g <- lm(update(formula(g), as.formula(paste(". ~ . -", worst))), data = data)
  }
}
g_back <- backward(lm(Life.Exp ~ ., data = statedata), data = statedata)
print(formula(g_back))

# Variables left out can still be related to the response
print(summary(lm(Life.Exp ~ Illiteracy + Murder + Frost, statedata)))

# What it shows (numbers from the book):
# - Full model: only Murder (p = 8.7e-08) and HS.Grad (p = 0.042) are
#   significant at 5%; R^2 = 0.736. Income is not significant, surprisingly.
# - Removal order: Area, Illiteracy, Income, Population (the last is a close
#   call, p = 0.052). Final model: Murder, HS.Grad, Frost with R^2 = 0.713,
#   barely below the full model's 0.736.
# - Illiteracy, dropped along the way, is significant (p = 0.05) in a model
#   with Murder and Frost: being dropped means "adds nothing once HS.Grad is
#   in", not "unrelated to life expectancy".
# - p-values after selection are too optimistic: removing the weaker
#   predictors inflates the apparent significance of those left.
