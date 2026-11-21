# Ridge regression, using the meatspec data (first 172 rows train, last 43 test)
# Ridge shrinks the (standardised) coefficients towards zero by adding lambda to the
# diagonal of X'X: beta = (X'X + lambda I)^-1 X'y

library(faraway)
library(MASS)      # lm.ridge(), select()
data(meatspec)
train <- 1:172
test <- 173:215
rmse <- function(x, y) sqrt(mean((x - y)^2))

# Centre the predictors with the training means, and centre the response
mm <- apply(meatspec[train, -101], 2, mean)
trainx <- as.matrix(sweep(meatspec[train, -101], 2, mm))
testx <- as.matrix(sweep(meatspec[test, -101], 2, mm))
yc <- meatspec$fat[train] - mean(meatspec$fat[train])

# lm.ridge(): ridge fits for a whole grid of lambda values (here 0 to 5e-8 in steps of 1e-9)
gridge <- lm.ridge(yc ~ trainx, lambda = seq(0, 5e-8, 1e-9))
# Ridge trace: each coefficient as a function of lambda
matplot(gridge$lambda, t(gridge$coef), type = "l", lty = 1,
        xlab = expression(lambda), ylab = expression(hat(beta)), main = "Ridge trace")

# MASS::select(): automatic choices of lambda (HKB, L-W and generalised cross-validation).
# Written MASS::select because other packages (e.g. dplyr) also have a select()
MASS::select(gridge)
best <- which.min(gridge$GCV)        # the grid position with the smallest GCV
print(gridge$lambda[best])
abline(v = gridge$lambda[best])

# lm.ridge scales the predictors internally, so do the same (scale = gridge$scales)
# before multiplying by the coefficients, and add back the mean of the response
ypredg <- scale(trainx, center = FALSE, scale = gridge$scales) %*% gridge$coef[, best] +
  mean(meatspec$fat[train])
print(rmse(ypredg, meatspec$fat[train]))
ytpredg <- scale(testx, center = FALSE, scale = gridge$scales) %*% gridge$coef[, best] +
  mean(meatspec$fat[train])
print(rmse(ytpredg, meatspec$fat[test]))

# One test case is predicted very badly; compare it with the truth
bad <- which.max(abs(ytpredg - meatspec$fat[test]))
print(c(case = 172 + bad, ridge = ytpredg[bad], truth = meatspec$fat[172 + bad]))
# Test RMSE without that case
print(rmse(ytpredg[-bad], meatspec$fat[test][-bad]))

# What it shows (numbers from the book):
# - GCV chooses lambda = 1.8e-08 (grid position 19). The trace shows the
#   coefficients settling down quickly as lambda increases from 0.
# - Training RMSE 0.80, but test RMSE 4.10: dismayingly poor.
# - Almost all of that comes from one sample, case 185: ridge predicts 11.2
#   when the truth is 34.8 (PLS predicted 35.7). Without it the test RMSE is
#   1.98, as good as PLS.
# - Ridge trades a little bias for a big reduction in variance. Which of
#   PCR, PLS and ridge is best varies from dataset to dataset.
