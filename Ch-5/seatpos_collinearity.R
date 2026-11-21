# Collinearity, using the seatpos data (38 car drivers, University of Michigan HuMoSim lab)
# Columns: Age (years), Weight (lb), HtShoes / Ht (height with / without shoes, cm), Seated
# (seated height), Arm, Thigh, Leg (lengths, cm), hipcenter (horizontal hip position, mm)

library(faraway)
data(seatpos)

g <- lm(hipcenter ~ ., seatpos)
print(summary(g))

# Pairwise correlations, rounded to 3 decimals
print(round(cor(seatpos), 3))


# ---- Condition numbers ----
# model.matrix(g)[, -1]: the X matrix without the intercept column
x <- model.matrix(g)[, -1]
# eigen(): eigenvalues of X'X. A big spread between the largest and smallest means collinearity
e <- eigen(t(x) %*% x)
print(e$values)
# Condition numbers kappa_j = sqrt(lambda_1 / lambda_j); above 30 is considered large
print(sqrt(e$values[1] / e$values))


# ---- Variance inflation factors ----
# VIF for Age by hand: R^2 from regressing Age on the other predictors, then 1 / (1 - R^2)
r2 <- summary(lm(x[, 1] ~ x[, -1]))$r.squared
print(r2)
print(1 / (1 - r2))
# vif(): (faraway) all the VIFs at once
print(vif(x))
# sqrt(VIF) = how many times larger the SE is than with uncorrelated predictors
print(sqrt(vif(x)["HtShoes"]))


# ---- Instability: add a small amount of noise to the response ----
set.seed(1)
# hipcenter + 10 * rnorm(38): measurement error with SD 10 mm
g_noise <- lm(hipcenter + 10 * rnorm(38) ~ ., seatpos)
print(summary(g_noise))
# Compare the coefficients side by side
print(round(cbind(original = coef(g), with_noise = coef(g_noise)), 3))


# ---- The length variables are all measuring body size ----
print(round(cor(x[, 3:8]), 2))
# Keep just one of them (height, the simplest to measure)
g2 <- lm(hipcenter ~ Age + Weight + Ht, seatpos)
print(summary(g2))

# What it shows (numbers from the book):
# - Classic signs of collinearity: R^2 = 0.687 and the overall F-test is very
#   significant (p = 1.3e-05), yet no single predictor is significant.
# - HtShoes and Ht are almost the same variable (r = 0.998), and the six
#   length variables are all strongly correlated with each other.
# - The eigenvalues range from 3.65e+06 down to 7.3, giving condition numbers
#   up to 708. Several are above 30, so more than one linear combination is
#   causing trouble.
# - VIFs: Age 2.0, Weight 3.6, HtShoes 307, Ht 333, Seated 9.0, Arm 4.5,
#   Thigh 2.8, Leg 6.7. sqrt(307) = 17.5, so the HtShoes SE is about 17.5
#   times what it would be with uncorrelated predictors.
# - Adding noise with SD 10 mm barely changes R^2, but several coefficients
#   change a lot (in the book, Seated flipped sign from 0.53 to -0.75; your
#   random noise will give different changes): the estimates are unstable.
# - Keeping only Age, Weight and Ht: R^2 = 0.656, almost as good as with all
#   eight, and Ht is now highly significant (p = 0.00017). The dropped
#   variables are not unrelated to hipcenter; they are just redundant once Ht is in.
