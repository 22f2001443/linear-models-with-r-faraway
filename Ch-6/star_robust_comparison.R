# Least squares vs Huber vs least trimmed squares on the star data (see Chapter 4)

library(faraway)
library(MASS)
data(star)

plot(light ~ temp, star, main = "Three fits to the star data")
gs1 <- lm(light ~ temp, star)                         # least squares
abline(coef(gs1))
gs2 <- rlm(light ~ temp, star)                        # Huber M-estimate
abline(coef(gs2), lty = 2)
gs3 <- ltsreg(light ~ temp, star, nsamp = "exact")    # least trimmed squares, exact search
abline(coef(gs3), lty = 5)
legend("bottomleft", legend = c("least squares", "Huber", "LTS"), lty = c(1, 2, 5))

print(rbind(LS = coef(gs1), Huber = coef(gs2), LTS = coef(gs3)))

# What it shows:
# - Only LTS follows the main group of stars. Huber is almost the same as
#   least squares: an M-estimate down-weights large residuals, but the four
#   giants have high leverage, pull the line towards themselves, and so do
#   not get large residuals.
# - LTS is "resistant": it fits roughly half the data and ignores the rest,
#   so a cluster of bad high-leverage points cannot drag it away.
