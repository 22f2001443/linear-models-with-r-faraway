# Partial least squares (PLS), using the meatspec data (first 172 rows train, last 43 test)
# Package: pls (install once with install.packages("pls")). The book's pls.pcr package
# no longer exists; plsr() from pls does the same job.

library(faraway)
library(pls)
data(meatspec)
train <- 1:172
test <- 173:215
rmse <- function(x, y) sqrt(mean((x - y)^2))

set.seed(123)
# plsr(): PLS regression with up to 50 components and 10-fold cross-validation.
# Unlike PCR, PLS chooses each component using the response as well as the predictors
plsmod <- plsr(fat ~ ., data = meatspec[train, ], ncomp = 50, validation = "CV")

par(mfrow = c(1, 2))
# coef(..., ncomp = 4): the implied coefficients on the 100 frequencies for 4 components
plot(drop(coef(plsmod, ncomp = 4)), ylab = "Coefficient", main = "PLS, 4 components")
cv <- RMSEP(plsmod, estimate = "CV")
plot(cv, main = "CV RMSE, PLS")
par(mfrow = c(1, 1))
ncomp_cv <- which.min(cv$val[1, 1, ]) - 1
print(ncomp_cv)

# Training and test performance, for the CV choice and for the book's 14 components
# fitted(plsmod)[, 1, k]: fitted values using k components
for (k in unique(c(ncomp_cv, 14))) {
  cat("components:", k,
      " train RMSE:", round(rmse(fitted(plsmod)[, 1, k], meatspec$fat[train]), 4),
      " test RMSE:", round(rmse(drop(predict(plsmod, meatspec[test, ], ncomp = k)),
                                meatspec$fat[test]), 4), "\n")
}

# What it shows (numbers from the book; CV is random):
# - The 4-component coefficients look much like the PCR ones: a smooth,
#   shrunken combination of the frequencies.
# - CV suggests about 14 components, half as many as PCR needed (21), because
#   PLS uses the response to build its components.
# - With 14 components: training RMSE 1.95 and test RMSE 2.01, slightly better
#   than PCR's 2.21.
# - PLS is good for prediction but nearly useless for explanation, and PLS and
#   PCR still need the usual diagnostics.
