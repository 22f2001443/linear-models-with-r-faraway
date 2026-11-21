# library(): loads the faraway package, which contains the pima dataset
library(faraway)
# data(): loads the pima dataset from the package into the workspace
data(pima)

# print(): shows the result in the console
# names(): returns the column names of the data frame
print(names(pima))     # all column names

# head(): returns the first 6 rows of the data frame
print(head(pima))      # first 6 rows

# factor(): converts the 0/1 test column into a categorical variable
pima$test <- factor(pima$test)

# summary(): gives min, quartiles, mean and max for numeric columns, and counts for factors
print(summary(pima))   # summary statistics for each column
pima                   # typing the object name alone prints the whole data frame



# plot(): with one vector, draws each glucose value against its row index
plot(pima$glucose)

# sort(): orders the glucose values from smallest to largest (the default)
print(sort(pima$glucose))
# sort(decreasing = TRUE): orders the diastolic values from largest to smallest
print(sort(pima$diastolic, decreasing = TRUE))

# plot() with a formula (y ~ x): scatter plot of glucose (y-axis) against diastolic (x-axis)
# xlab / ylab set the axis labels, main sets the plot title
plot(glucose ~ diastolic, data = pima,
     xlab = "Diastolic blood pressure (mm Hg)",
     ylab = "Plasma glucose",
     main = "Glucose vs diastolic")


# Replace impossible 0 readings with NA (missing) so they don't distort plots and summaries
pima$diastolic [pima$diastolic == 0] <- NA
pima$glucose [pima$glucose == 0] <- NA
pima$triceps [pima$triceps == 0] <- NA
pima$insulin [pima$insulin == 0] <- NA
pima$bmi [pima$bmi == 0] <- NA

# View(): opens the data frame in a spreadsheet-style viewer (RStudio / VS Code)
View(pima)

# hist(): draws a histogram of diastolic pressure
# freq = TRUE puts counts on the y-axis, col sets the bar colour
hist (pima$diastolic, freq = TRUE,
      main = "Histogram of diastolic blood pressure",
      xlab = "Diastolic blood pressure (mm Hg)",
      col = "lightblue", )
# density(): estimates a smooth curve of the distribution (na.rm = TRUE skips the NA values)
# lines(): adds that curve to the existing plot as a red line, lwd = 2 makes it thicker
lines(density(pima$diastolic, na.rm = TRUE), col = "red", lwd = 2)


# plot() on a density() result: draws the smooth density curve on its own axes
plot(density(pima$diastolic, na.rm = TRUE),
     main = "Density of diastolic blood pressure",
     xlab = "Diastolic blood pressure (mm Hg)")

# First block (histogram plus curve): draws the bars, then adds the curve on top as a red line. Use it to see the raw counts and the smoothed shape together.

# Second block (curve only): draws just the curve on its own axes, with no bars. It gives a cleaner view of the distribution's shape.

# Cumulative frequency histogram: each bar = number of people at or below that value
# hist(plot = FALSE): computes the bins and counts without drawing anything, and stores them in h
h <- hist(pima$diastolic, plot = FALSE)
# cumsum(): replaces each bin's count with the running total up to that bin
h$counts <- cumsum(h$counts)
# plot() on a histogram object: draws the bars using the cumulative counts
plot(h,
     main = "Cumulative frequency of diastolic blood pressure",
     xlab = "Diastolic blood pressure (mm Hg)",
     ylab = "Cumulative frequency",
     col = "lightgreen")
# c(): puts 0 in front of the counts so there is one y value for every bin edge in h$breaks
# lines(): joins those points with a red line on top of the bars
lines(h$breaks, c(0, h$counts), col = "red", lwd = 2)   # cumulative frequency curve (ogive)
# sort() + plot(): plots the sorted diastolic values against their rank, pch = "." draws each point as a small dot
plot (sort (pima$diastolic), pch=".")

# plot() with a formula: scatter plot of the diabetes pedigree function against diastolic pressure
plot (diabetes  ~ diastolic,pima)
# plot() with a factor on the x-axis: draws a boxplot of diabetes for each test group (0 and 1)
plot (diabetes  ~ test, pima)

# pairs(): draws a grid of scatter plots, one for every pair of columns in pima
pairs (pima)
