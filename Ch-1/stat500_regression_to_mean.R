# Regression to the mean, using the stat500 exam scores (55 students)
# Columns: midterm, final, hw (homework) and total course score

# library(): loads the faraway package, which contains the stat500 dataset
library(faraway)
# data(): loads the stat500 dataset into the workspace
data(stat500)

# scale(): standardises every column (subtracts the mean, divides by the standard deviation),
# so each column has mean 0 and sd 1. scale() returns a matrix, so data.frame() turns it back into a data frame
stat500 <- data.frame(scale(stat500))

# plot() with a formula (y ~ x): scatter plot of final (y-axis) against midterm (x-axis)
plot(final ~ midterm, stat500)
# abline(0, 1): adds the line y = x (intercept 0, slope 1).
# A student on this line did equally well, relative to the class, on both exams
abline(0, 1)

# lm(): fits the least-squares regression line of final on midterm
g <- lm(final ~ midterm, stat500)
# coef(): pulls out the fitted intercept and slope from g
# abline(): draws that fitted line; lty = 5 makes it dashed
abline(coef(g), lty = 5)

# cor(): correlation matrix of all four columns; print() shows it in the console
print(cor(stat500))

# What it shows:
# - Because both scores are standardised, the slope of the regression line equals
#   the correlation between midterm and final (0.545), so it is flatter than y = x.
# - A student 1 sd above average on the midterm is predicted to be only about
#   0.55 sd above average on the final. Equally, a student well below average is
#   predicted to do less badly. Extreme scores tend to move back towards the mean:
#   this is "regression to the mean", where the word "regression" comes from.
# - Homework is barely related to the final (r = 0.09) and only weakly to the
#   midterm (r = 0.27), but all three feed into total.
