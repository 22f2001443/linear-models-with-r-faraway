# =============================================================================
# Linear Models with R (Faraway), Chapter 8 exercises 1-5
# Backward elimination, AIC, adjusted R^2 and Mallows Cp on prostate, teengamb
# and divusa; simplifying a quadratic model for trees; selection with and
# without influential points for stackloss.
# Packages: leaps (install.packages("leaps")).
# =============================================================================

library(faraway)
library(leaps)

# -----------------------------------------------------------------------------
# backward(g, data, alpha): backward elimination - drop the term with the largest
# p-value above alpha, refit, repeat (same function as in the chapter script).
# data = the data frame g was fitted to (passed in so the refit can find it)
# -----------------------------------------------------------------------------
backward <- function(g, data, alpha = 0.05) {
  repeat {
    d <- drop1(g, test = "F")
    pv <- d[-1, "Pr(>F)"]
    if (length(pv) == 0 || max(pv) <= alpha) return(g)
    worst <- rownames(d)[-1][which.max(pv)]
    cat("  dropping", worst, "(p =", signif(max(pv), 3), ")\n")
    # update(formula(g), ...): the current formula with that term removed; refit on data
    g <- lm(update(formula(g), as.formula(paste(". ~ . -", worst))), data = data)
  }
}

# -----------------------------------------------------------------------------
# select_all(formula, data, title): runs all four methods and prints the chosen models
# -----------------------------------------------------------------------------
select_all <- function(formula, data, title = "") {
  cat("\n==========", title, "==========\n")
  g <- lm(formula, data = data)

  cat("(a) Backward elimination (alpha = 0.05):\n")
  gb <- backward(g, data)
  print(formula(gb))

  cat("(b) AIC (step):\n")
  # trace = 0 hides the step-by-step printout
  ga <- step(g, trace = 0)
  print(formula(ga))

  # nvmax = number of predictors, so models of every size are searched
  k <- ncol(model.matrix(g)) - 1
  rs <- summary(regsubsets(formula, data = data, nvmax = k))
  tab <- data.frame(p = 2:(k + 1), RSS = rs$rss, adjR2 = rs$adjr2, Cp = rs$cp)
  print(round(tab, 3))

  cat("(c) Adjusted R^2 - best model:\n")
  w <- rs$which[which.max(rs$adjr2), -1]          # -1 drops the intercept column
  print(names(w)[w])

  cat("(d) Mallows Cp - smallest model with Cp <= p:\n")
  ok <- which(rs$cp <= tab$p)
  pick <- if (length(ok) > 0) ok[1] else which.min(rs$cp)
  w <- rs$which[pick, -1]
  print(names(w)[w])

  par(mfrow = c(1, 2))
  plot(tab$p, tab$Cp, xlab = "No. of parameters", ylab = "Cp", main = paste(title, "Cp"))
  abline(0, 1)
  plot(tab$p, tab$adjR2, xlab = "No. of parameters", ylab = "Adjusted R^2",
       main = paste(title, "adjusted R^2"))
  par(mfrow = c(1, 1))
  invisible(list(backward = gb, aic = ga, subsets = rs))
}


# =============================================================================
# Exercise 1: prostate - lpsa on all other variables
# =============================================================================
data(prostate)
sel_pr <- select_all(lpsa ~ ., prostate, "prostate")

# What stands out (prostate):
# [numbers to be filled in from your run_all.R output]
# - Compare the four chosen sets. Backward elimination at 5% usually gives the
#   smallest model; AIC and adjusted R^2 tend to keep one or two more variables
#   that help prediction without being individually significant.


# =============================================================================
# Exercise 2: teengamb - gamble on sex, status, income, verbal
# =============================================================================
data(teengamb)
sel_tg <- select_all(gamble ~ ., teengamb, "teengamb")

# What stands out (teengamb):
# [numbers to be filled in from your run_all.R output]


# =============================================================================
# Exercise 3: divusa - divorce on all other variables (including year)
# =============================================================================
data(divusa)
sel_dv <- select_all(divorce ~ ., divusa, "divusa")

# What stands out (divusa):
# [numbers to be filled in from your run_all.R output]
# - Remember Chapter 4 exercise 5: the errors are strongly serially
#   correlated, so all these p-value and RSS based criteria overstate how
#   well determined the model is.


# =============================================================================
# Exercise 4: trees - log(Volume) on a second-order polynomial in Girth and Height
# =============================================================================
data(trees)
g_full <- lm(log(Volume) ~ Girth + Height + I(Girth^2) + I(Height^2) + Girth:Height, trees)
print(summary(g_full))

# Respect the hierarchy: only the second-order terms (squares and the interaction) may be
# removed while they are present; the linear terms come after
second <- c("I(Girth^2)", "I(Height^2)", "Girth:Height")
g <- g_full
repeat {
  # attr(terms(g), "term.labels"): the names of the terms currently in the model
  present <- intersect(second, attr(terms(g), "term.labels"))
  if (length(present) == 0) break
  # drop1(..., scope = ...): only consider dropping the listed terms
  d <- drop1(g, scope = reformulate(present), test = "F")
  pv <- d[-1, "Pr(>F)"]
  if (max(pv) <= 0.05) break
  worst <- rownames(d)[-1][which.max(pv)]
  cat("dropping", worst, "(p =", signif(max(pv), 3), ")\n")
  g <- update(g, as.formula(paste(". ~ . -", worst)))
}
print(summary(g))
# Joint test: the full quadratic surface against the plane (all second-order terms at once)
g_plane <- lm(log(Volume) ~ Girth + Height, trees)
print(anova(g_plane, g_full))

# What stands out (trees):
# [numbers to be filled in from your run_all.R output]
# - The joint F-test says whether the curved surface is needed at all.
# - Girth^2 is the term most likely to survive: log(Volume) grows like
#   2 log(Girth), which is curved in Girth itself. Compare with the log-log
#   model from Chapter 7 exercise 5, which needs no second-order terms.


# =============================================================================
# Exercise 5: stackloss - stack.loss on Air.Flow, Water.Temp, Acid.Conc.
#   (The question says "stack.loss as the predictor"; it means the response.)
# =============================================================================
data(stackloss)
sel_sl <- select_all(stack.loss ~ ., stackloss, "stackloss (all 21 days)")

# Outliers and influential points in the full model
g_sl <- lm(stack.loss ~ ., stackloss)
n <- nrow(stackloss); p <- 4
print(data.frame(day = 1:n, leverage = round(hatvalues(g_sl), 3),
                 jackknife = round(rstudent(g_sl), 2),
                 cook = round(cooks.distance(g_sl), 3))[order(-cooks.distance(g_sl)), ][1:5, ])
print(c(bonferroni = abs(qt(0.05 / (2 * n), n - p - 1)), cook_cutoff = 4 / n))
halfnorm(cooks.distance(g_sl), 3, ylab = "Cook's distances")

# Remove the influential days (Cook's distance above 4/n) and repeat the selection
keep_days <- cooks.distance(g_sl) <= 4 / n
print(which(!keep_days))
sel_sl2 <- select_all(stack.loss ~ ., stackloss[keep_days, ], "stackloss (influential days removed)")

# What stands out (stackloss):
# [numbers to be filled in from your run_all.R output]
# - With all days, Acid.Conc. is the predictor the methods are most likely to
#   drop. Day 21 is the classic influential point (with 1, 3 and 4 often
#   flagged as outliers).
# - Compare the selections with and without the influential days: variable
#   selection is sensitive to a handful of points, so diagnostics should come
#   before (and after) selection.
