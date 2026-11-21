# Principal components regression (PCR), using the meatspec data
# 215 meat samples: V1-V100 = absorbances at 100 near-infrared wavelengths, fat = fat content.
# Training sample = first 172 rows, test sample = last 43 rows.
# Packages: pls for the cross-validation part (install once with install.packages("pls")).
# The book used mva and pls.pcr, which no longer exist: prcomp() is now in base R
# and pls.pcr has been replaced by the pls package.

library(faraway)
data(meatspec)
train <- 1:172
test <- 173:215

# ---- Least squares with all 100 predictors ----
model1 <- lm(fat ~ ., meatspec[train, ])
print(summary(model1)$r.squared)
# rmse(): root mean squared error between predictions x and true values y
rmse <- function(x, y) sqrt(mean((x - y)^2))
print(rmse(fitted(model1), meatspec$fat[train]))                          # training error
print(rmse(predict(model1, meatspec[test, ]), meatspec$fat[test]))        # test error

# Stepwise AIC selection among the 100 predictors (takes a few seconds; trace = 0 hides the steps)
model2 <- step(model1, trace = 0)
print(length(coef(model2)) - 1)                                           # predictors kept
print(rmse(fitted(model2), meatspec$fat[train]))
print(rmse(predict(model2, meatspec[test, ]), meatspec$fat[test]))


# ---- Principal components of the training predictors ----
# prcomp(): principal components; it centres each column (column 101, fat, is left out)
meatpca <- prcomp(meatspec[train, -101])
# $sdev: SD of each component (square roots of the eigenvalues)
print(round(meatpca$sdev, 3))
# $rotation: the eigenvectors (loadings); plot the first three against wavelength
matplot(1:100, meatpca$rotation[, 1:3], type = "l", xlab = "Frequency", ylab = "",
        main = "First three PC loadings")

# Regression on the first four PCs. $x holds the PC scores
model3 <- lm(fat ~ meatpca$x[, 1:4], meatspec[train, ])
print(rmse(fitted(model3), meatspec$fat[train]))

# Shrinkage: the implied coefficients on the original 100 predictors
par(mfrow = c(1, 2))
plot(model1$coefficients[-1], ylab = "Coefficient", main = "Least squares")
svb <- meatpca$rotation[, 1:4] %*% model3$coefficients[-1]
plot(svb, ylab = "Coefficient", main = "PCR, 4 components")
par(mfrow = c(1, 1))

# Scree plot: SD of the first 10 PCs; look for an "elbow"
plot(meatpca$sdev[1:10], type = "l", ylab = "SD of PC", xlab = "PC number")

# ---- Predicting the test sample ----
# Centre the test predictors with the TRAINING means, then project onto the eigenvectors
mm <- apply(meatspec[train, -101], 2, mean)            # apply(..., 2, mean): column means
# sweep(x, 2, mm): subtracts mm from every row of x
tx <- as.matrix(sweep(meatspec[test, -101], 2, mm))
nx <- tx %*% meatpca$rotation[, 1:4]
pv <- cbind(1, nx) %*% model3$coefficients
print(rmse(pv, meatspec$fat[test]))

# How many PCs would be best for the test sample? (Cheating: in practice we don't know the test data)
rmsmeat <- numeric(50)
for (i in 1:50) {
  nx <- tx %*% meatpca$rotation[, 1:i, drop = FALSE]
  model3 <- lm(fat ~ meatpca$x[, 1:i], meatspec[train, ])
  pv <- cbind(1, nx) %*% model3$coefficients
  rmsmeat[i] <- rmse(pv, meatspec$fat[test])
}
plot(rmsmeat, ylab = "Test RMS", xlab = "No. of components")
print(which.min(rmsmeat))
print(min(rmsmeat))


# ---- Choosing the number of components honestly: cross-validation ----
library(pls)
set.seed(123)       # CV splits the data at random
# pcr(): principal components regression with up to 50 components;
# validation = "CV" uses 10-fold cross-validation
pcrmod <- pcr(fat ~ ., data = meatspec[train, ], ncomp = 50, validation = "CV")
# RMSEP(): cross-validated RMSE for 0, 1, ..., 50 components
cv <- RMSEP(pcrmod, estimate = "CV")
plot(cv, main = "CV RMSE, PCR")
# $val[1, 1, ] holds the values; the first is the 0-component model, so subtract 1
ncomp_cv <- which.min(cv$val[1, 1, ]) - 1
print(ncomp_cv)
# drop(): turns the 3-d prediction array into a plain vector
print(rmse(drop(predict(pcrmod, meatspec[test, ], ncomp = ncomp_cv)), meatspec$fat[test]))

# What it shows (numbers from the book; the CV part is random):
# - Least squares: R^2 = 0.997 and training RMSE 0.69, but test RMSE 3.81, about
#   five times worse. The fit to the training data is over-optimistic.
# - Stepwise AIC removes 28 variables: test RMSE improves a little, to 3.59.
# - The first PC has SD 5.06, the second 0.51: about ten times smaller, then a
#   sharp drop. PC1 loads almost equally on every frequency (overall
#   absorbance level), PC2 contrasts high and low frequencies.
# - Four PCs: training RMSE 4.06, test RMSE 4.53 - not good. But the implied
#   coefficients range only from about -10 to 10 instead of thousands, and
#   vary smoothly across frequencies: that is the shrinkage.
# - Using the test data to choose (not allowed in practice), 27 PCs give the
#   best test RMSE, 1.85.
# - Cross-validation in the book chose 21 components, with test RMSE 2.21,
#   close to that best value. Your CV choice may differ slightly with the
#   random folds and the newer pls package.
