# Criterion-based variable selection (AIC, adjusted R^2, Mallows Cp), using the state data
# Packages: leaps for regsubsets() - install once with install.packages("leaps")

data(state)
statedata <- data.frame(state.x77, row.names = state.abb)
g <- lm(Life.Exp ~ ., data = statedata)

# ---- AIC ----
# step(): starting from the full model, repeatedly drops (or adds back) the term that
# lowers AIC = n log(RSS/n) + 2p the most, until no change helps
g_aic <- step(g)
print(formula(g_aic))
print(coef(g_aic))


# ---- Best subsets: adjusted R^2 and Mallows Cp ----
library(leaps)
# regsubsets(): for each model size, finds the subset of predictors with the smallest RSS
b <- regsubsets(Life.Exp ~ ., data = statedata)
rs <- summary(b)
print(rs)
# rs$which: TRUE/FALSE table of which variables are in the best model of each size
print(rs$which)

par(mfrow = c(1, 2))
# Cp against the number of parameters p; good models lie on or below the line Cp = p
plot(2:8, rs$cp, xlab = "No. of Parameters", ylab = "Cp Statistic")
abline(0, 1)
# Adjusted R^2 against p: pick the largest
plot(2:8, rs$adjr2, xlab = "No. of Parameters", ylab = "Adjusted R-square")
par(mfrow = c(1, 1))
print(rs$which[which.max(rs$adjr2), ])
print(round(cbind(p = 2:8, Cp = rs$cp, adjR2 = rs$adjr2), 3))


# ---- Sensitivity to unusual points ----
# lm.influence()$hat: the leverages; rev(sort()) lists them from largest down
h <- lm.influence(g)$hat
names(h) <- state.abb
print(round(rev(sort(h))[1:7], 3))

# Without Alaska (the highest-leverage state)
b <- regsubsets(Life.Exp ~ ., data = statedata[state.abb != "AK", ])
rs <- summary(b)
print(rs$which[which.max(rs$adjr2), ])

# ---- Sensitivity to transformations ----
# stripchart(): one column of jittered points per (standardised) variable, to spot skewness
stripchart(data.frame(scale(statedata)), vertical = TRUE, method = "jitter", las = 2)
# Population and Area are skewed, so try their logs
b <- regsubsets(Life.Exp ~ log(Population) + Income + Illiteracy + Murder + HS.Grad +
                  Frost + log(Area), data = statedata)
rs <- summary(b)
print(rs$which[which.max(rs$adjr2), ])
print(max(rs$adjr2))

# What it shows (numbers from the book):
# - AIC removes Area, Illiteracy and Income (the same order as backward
#   elimination) but keeps Population: Population + Murder + HS.Grad + Frost.
# - Best subsets: Murder for one predictor, + HS.Grad for two, + Frost for
#   three, + Population for four. The three- and four-predictor models both
#   lie on or below Cp = p. Adjusted R^2 is largest for the four-predictor
#   model (Population, Murder, HS.Grad, Frost).
# - Alaska has leverage 0.81. Without it, the adjusted R^2 choice also brings
#   in Area: the selected model depends on one state.
# - Using log(Population) and log(Area) changes the best model to
#   log(Population), Murder, HS.Grad, Frost, with adjusted R^2 = 0.717, the
#   highest so far.
# - Several models fit about equally well: use the criteria as a guide, then
#   compare their predictions, diagnostics and interpretability.
