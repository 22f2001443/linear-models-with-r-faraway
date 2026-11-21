# Subset regression vs broken stick (segmented) regression, using savings: sr on pop15 only

library(faraway)
data(savings)

# Separate straight lines for countries below and above pop15 = 35
g1 <- lm(sr ~ pop15, savings, subset = (pop15 < 35))
g2 <- lm(sr ~ pop15, savings, subset = (pop15 > 35))
plot(sr ~ pop15, savings, xlab = "Pop'n under 15", ylab = "Savings Rate")
abline(v = 35, lty = 5)              # the dividing line
# segments(x0, y0, x1, y1): draws a line segment from (x0, y0) to (x1, y1)
segments(20, g1$coef[1] + g1$coef[2] * 20, 35, g1$coef[1] + g1$coef[2] * 35)
segments(48, g2$coef[1] + g2$coef[2] * 48, 35, g2$coef[1] + g2$coef[2] * 35)

# Broken stick: two "hockey stick" basis functions that meet at c = 35
# ifelse(test, a, b): a where the test is TRUE, b where it is FALSE
lhs <- function(x) ifelse(x < 35, 35 - x, 0)   # left arm: slopes down to 0 at 35
rhs <- function(x) ifelse(x < 35, 0, x - 35)   # right arm: 0 until 35, then rises
gb <- lm(sr ~ lhs(pop15) + rhs(pop15), savings)
print(summary(gb))
x <- seq(20, 48, by = 1)
py <- gb$coef[1] + gb$coef[2] * lhs(x) + gb$coef[3] * rhs(x)
lines(x, py, lty = 2)
legend("topright", legend = c("separate fits", "broken stick"), lty = c(1, 2))

# Slopes of the two separate fits vs the broken stick (the lhs slope is minus the coefficient)
print(rbind(separate = c(left = unname(coef(g1)[2]), right = unname(coef(g2)[2])),
            broken_stick = c(left = unname(-coef(gb)[2]), right = unname(coef(gb)[3]))))

# What it shows:
# - The two separate lines do not meet at 35. The broken stick forces them to
#   join there and uses 3 parameters instead of 4. Its intercept is the
#   fitted savings rate at the join (pop15 = 35).
# - Imposing continuity changes the sign of the slope for the high-pop15
#   countries. The two groups are so different, with so few countries near
#   35, that continuity is not obviously a sensible assumption here.
