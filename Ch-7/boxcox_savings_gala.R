# The Box-Cox method for choosing a response transformation, using the savings and gala data
# Box-Cox looks at y^lambda (log y when lambda = 0) and picks lambda by maximum likelihood.

library(faraway)
library(MASS)      # boxcox()

# ---- savings ----
data(savings)
g <- lm(sr ~ pop15 + pop75 + dpi + ddpi, savings)
par(mfrow = c(1, 2))
# boxcox(): plots the profile log-likelihood of lambda with a 95% confidence interval
boxcox(g, plotit = TRUE)
# lambda = seq(0.5, 1.5, by = 0.1): zoom in to read the interval more easily
bc <- boxcox(g, plotit = TRUE, lambda = seq(0.5, 1.5, by = 0.1))
par(mfrow = c(1, 1))

# Read the numbers off instead of the plot: boxcox() returns $x (lambda) and $y (log-likelihood).
# The 95% CI is every lambda whose log-likelihood is within qchisq(0.95, 1) / 2 of the maximum
bc <- boxcox(g, plotit = FALSE, lambda = seq(0, 2, by = 0.01))
print(bc$x[which.max(bc$y)])                                  # best lambda
print(range(bc$x[bc$y > max(bc$y) - qchisq(0.95, 1) / 2]))    # 95% CI


# ---- gala ----
data(gala)
g <- lm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent, gala)
par(mfrow = c(1, 2))
boxcox(g, plotit = TRUE)
boxcox(g, lambda = seq(0.0, 1.0, by = 0.05), plotit = TRUE)
par(mfrow = c(1, 1))
bc <- boxcox(g, plotit = FALSE, lambda = seq(-0.5, 1.5, by = 0.01))
print(bc$x[which.max(bc$y)])
print(range(bc$x[bc$y > max(bc$y) - qchisq(0.95, 1) / 2]))

# What it shows:
# - savings: the 95% CI for lambda runs from about 0.6 to 1.4 and includes
#   1 (no transformation), so there is no good reason to transform.
# - gala: the CI excludes 1, so the response clearly needs transforming. The
#   best lambda is near 1/3 (a cube root); 0.5 (square root, as used in
#   Chapter 4) is just inside the interval.
# - Round lambda to an interpretable value inside the interval. Box-Cox is
#   upset by outliers, needs y > 0, and has little effect when max(y)/min(y)
#   is small.
