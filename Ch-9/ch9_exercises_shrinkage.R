# =============================================================================
# Linear Models with R (Faraway), Chapter 9 exercises 1-4
# PCR, PLS and ridge regression on seatpos, and a comparison of five methods
# on the body fat data (fat) with a held-out test sample.
# Packages: pls (install.packages("pls")) and MASS (comes with R).
# =============================================================================

library(faraway)
library(pls)
library(MASS)

rmse <- function(x, y) sqrt(mean((x - y)^2))

data(seatpos)
# The new driver from exercise 1 (all eight predictors, in the data's column order)
x0 <- data.frame(Age = 64.800, Weight = 263.700, HtShoes = 181.080, Ht = 178.560,
                 Seated = 91.440, Arm = 35.640, Thigh = 40.950, Leg = 38.790)


# =============================================================================
# Exercise 1: PCR on seatpos
# =============================================================================
lengths <- c("HtShoes", "Ht", "Seated", "Arm", "Thigh", "Leg")

# --- (i) the six length variables (all in cm, so centre only, no scaling) ---
pc6 <- prcomp(seatpos[, lengths])
print(summary(pc6))                       # SD and proportion of variance for each PC
print(round(pc6$rotation[, 1:2], 3))      # loadings of the first two PCs
plot(pc6$sdev, type = "b", xlab = "PC number", ylab = "SD of PC", main = "Scree: 6 lengths")

set.seed(123)
pcr6 <- pcr(hipcenter ~ HtShoes + Ht + Seated + Arm + Thigh + Leg, data = seatpos,
            validation = "CV")
cv6 <- RMSEP(pcr6, estimate = "CV")
print(cv6)
# selectNcomp(..., "onesigma"): smallest number of components whose CV error is within
# one standard error of the best - a common way to avoid over-fitting the CV curve
k6 <- selectNcomp(pcr6, method = "onesigma")
print(k6)
k6 <- max(k6, 1)                          # use at least one component
print(drop(predict(pcr6, x0, ncomp = k6)))

# --- (ii) add Age and Weight: different units, so scale = TRUE (standardise) ---
pc8 <- prcomp(seatpos[, -9], scale. = TRUE)          # column 9 is hipcenter
print(summary(pc8))
print(round(pc8$rotation[, 1:3], 3))
set.seed(123)
pcr8 <- pcr(hipcenter ~ ., data = seatpos, scale = TRUE, validation = "CV")
print(RMSEP(pcr8, estimate = "CV"))
k8 <- max(selectNcomp(pcr8, method = "onesigma"), 1)
print(k8)
print(drop(predict(pcr8, x0, ncomp = k8)))
# The new driver's values next to the range of each predictor in the data
print(rbind(new_driver = unlist(x0), min = sapply(seatpos[, -9], min), max = sapply(seatpos[, -9], max)))

# What stands out (PCR):
# [numbers to be filled in from your run_all.R output]
# - PC1 of the lengths has loadings of the same sign on all six: it is overall
#   body size and carries most of the variance. Later PCs are contrasts
#   between body parts (e.g. leg vs trunk length).
# - One component is likely to be enough: tall drivers sit further back.
# - With Age and Weight added, PC1 is still "size" (Weight joins it); Age
#   probably gets a component of its own, being almost unrelated to the lengths.
# - Before trusting either prediction, check the new driver's values against
#   the range of the data (printed just above).


# =============================================================================
# Exercise 2: PLS on seatpos with all predictors
# =============================================================================
set.seed(123)
plsmod <- plsr(hipcenter ~ ., data = seatpos, scale = TRUE, validation = "CV")
print(RMSEP(plsmod, estimate = "CV"))
plot(RMSEP(plsmod, estimate = "CV"), main = "PLS: CV RMSE")
kp <- max(selectNcomp(plsmod, method = "onesigma"), 1)
print(kp)
print(drop(predict(plsmod, x0, ncomp = kp)))

# What stands out (PLS):
# [numbers to be filled in from your run_all.R output]
# - PLS should need no more components than PCR, because it uses hipcenter
#   when building them.


# =============================================================================
# Exercise 3: ridge regression on seatpos with all predictors
# =============================================================================
rg <- lm.ridge(hipcenter ~ ., data = seatpos, lambda = seq(0, 50, by = 0.1))
matplot(rg$lambda, t(rg$coef), type = "l", lty = 1, xlab = expression(lambda),
        ylab = "standardised coefficient", main = "Ridge trace: seatpos")
MASS::select(rg)
best <- which.min(rg$GCV)
print(rg$lambda[best])
abline(v = rg$lambda[best], lty = 2)
# coef(): coefficients on the original scale, with the intercept, one row per lambda
b <- coef(rg)[best, ]
print(round(b, 4))
# Prediction: intercept + sum of coefficient x value (x0 is in the same column order)
print(sum(b * c(1, unlist(x0))))
# Least squares prediction for comparison
print(predict(lm(hipcenter ~ ., seatpos), x0, interval = "prediction"))

# What stands out (ridge):
# [numbers to be filled in from your run_all.R output]
# - The trace shows the wildly unstable least squares coefficients (lambda = 0)
#   settling quickly; GCV picks the amount of shrinkage.
# - Compare the ridge, PCR, PLS and least squares predictions for the new driver.


# =============================================================================
# Exercise 4: fat - percentage body fat (Brozek) for 252 men
#   Response: brozek. siri (another body-fat formula), density (from which brozek is
#   computed) and free (fat-free weight, computed from brozek) are left out: they are
#   functions of the response, not genuine predictors.
#   Test sample = every 10th row, training sample = the rest.
# =============================================================================
data(fat)
fat2 <- fat[, !(names(fat) %in% c("siri", "density", "free"))]
test_i <- seq(10, nrow(fat2), by = 10)
train <- fat2[-test_i, ]
test <- fat2[test_i, ]
print(c(train = nrow(train), test = nrow(test)))

results <- c()

# (a) Linear regression with all predictors
m_all <- lm(brozek ~ ., train)
results["all predictors"] <- rmse(predict(m_all, test), test$brozek)

# (b) Variables selected by AIC
m_aic <- step(m_all, trace = 0)
print(formula(m_aic))
results["AIC selection"] <- rmse(predict(m_aic, test), test$brozek)

# (c) Principal component regression (standardised, CV choice of components)
set.seed(123)
m_pcr <- pcr(brozek ~ ., data = train, scale = TRUE, validation = "CV")
k_pcr <- max(selectNcomp(m_pcr, method = "onesigma"), 1)
results["PCR"] <- rmse(drop(predict(m_pcr, test, ncomp = k_pcr)), test$brozek)

# (d) Partial least squares
set.seed(123)
m_pls <- plsr(brozek ~ ., data = train, scale = TRUE, validation = "CV")
k_pls <- max(selectNcomp(m_pls, method = "onesigma"), 1)
results["PLS"] <- rmse(drop(predict(m_pls, test, ncomp = k_pls)), test$brozek)

# (e) Ridge regression, lambda chosen by GCV
m_rg <- lm.ridge(brozek ~ ., data = train, lambda = seq(0, 20, by = 0.05))
b_rg <- coef(m_rg)[which.min(m_rg$GCV), ]
# Test predictions: intercept + X_test %*% slopes (columns in the same order as the fit)
xt <- as.matrix(test[, names(b_rg)[-1]])
results["ridge"] <- rmse(b_rg[1] + xt %*% b_rg[-1], test$brozek)

print(c(PCR_components = k_pcr, PLS_components = k_pls,
        ridge_lambda = m_rg$lambda[which.min(m_rg$GCV)]))
print(round(sort(results), 3))
# barplot(): one bar per method
barplot(sort(results), ylab = "Test RMSE", main = "Body fat: test-sample RMSE", las = 2)

# What stands out (fat):
# [numbers to be filled in from your run_all.R output]
# - The test RMSEs are the fair comparison: none of the 25 test men was used
#   to build any model.
# - With only 14 predictors and 227 training cases, the methods are likely to
#   be close to each other. Shrinkage helps most when there are many
#   collinear predictors (like meatspec), less so here.
# - With a test sample of only 25 men, small differences in RMSE are not
#   meaningful.
