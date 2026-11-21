# =============================================================================
# Linear Models with R (Faraway), Chapter 1 exercises 1-5
# Numerical and graphical summaries of teengamb, uswages, prostate, sat, divusa.
# Each section: load and clean the data, a short numerical summary, one panel
# of plots, and comments on the features worth noticing.
# =============================================================================

library(faraway)   # library(): loads the faraway package, which contains all five datasets


# =============================================================================
# Exercise 1: teengamb - teenage gambling in Britain (47 teenagers)
#   sex     0 = male, 1 = female
#   status  socio-economic status score based on parents' occupation
#   income  pounds per week
#   verbal  verbal score: words out of 12 correctly defined
#   gamble  expenditure on gambling, pounds per year
# =============================================================================
data(teengamb)   # data(): loads the dataset into the workspace

# cor(): correlation matrix of all columns; round() keeps 2 decimals
# (run before sex becomes a factor, because cor() needs numbers)
print(round(cor(teengamb), 2))

# factor(): turns the 0/1 sex code into a labelled categorical variable
teengamb$sex <- factor(teengamb$sex, labels = c("male", "female"))

# summary(): min, quartiles, mean and max for numeric columns, counts for factors
print(summary(teengamb))
# tapply(): runs summary() on gamble separately for males and females
print(tapply(teengamb$gamble, teengamb$sex, summary))

# par(mfrow): splits the plot window into a 2 x 2 grid
par(mfrow = c(2, 2))
# hist(): distribution of gambling expenditure
hist(teengamb$gamble, main = "Gambling expenditure",
     xlab = "gamble (pounds/year)", col = "lightblue")
# boxplot() with a formula: one box of gamble for each sex
boxplot(gamble ~ sex, data = teengamb, main = "Gambling by sex",
        ylab = "gamble (pounds/year)", col = c("lightblue", "pink"))
# plot() with a formula: scatter plot; as.numeric(sex) gives colour 1 (black) for males, 2 (red) for females
plot(gamble ~ income, data = teengamb, col = as.numeric(sex), pch = 19,
     main = "Gambling vs income", xlab = "income (pounds/week)")
# legend(): adds a key for the two colours
legend("topleft", legend = levels(teengamb$sex), col = 1:2, pch = 19)
plot(gamble ~ verbal, data = teengamb, col = as.numeric(sex), pch = 19,
     main = "Gambling vs verbal score")
par(mfrow = c(1, 1))   # reset to one plot per window

# What stands out (teengamb):
# - 28 males and 19 females. gamble is very right-skewed: median 6, mean 19.3,
#   maximum 156, and 4 teenagers spent nothing at all.
# - Sex makes a big difference: male median 14.25 vs female median 1.70.
#   No female spent more than 19.6 a year, while a quarter of males spent over 42.
# - income is the variable most related to gamble (r = 0.62). The spread of
#   gamble also grows with income (sd about 10 when income <= 3, about 47 when
#   income > 6): a funnel shape to keep in mind for any regression model.
# - status and verbal are only weakly related to gamble (r = -0.05 and -0.22),
#   but are related to each other (r = 0.53).


# =============================================================================
# Exercise 2: uswages - weekly wages of 2000 US male workers, CPS 1988
#   wage          real weekly wage in dollars
#   educ, exper   years of education and of experience
#   race          1 = Black, 0 = White
#   smsa          1 = lives in a metropolitan area
#   ne, mw, so, we  region dummies (north-east, mid-west, south, west)
#   pt            1 = works part time
# =============================================================================
data(uswages)

# Negative experience is impossible (exper is worked out as age - educ - 6), so mark it as missing
uswages$exper[uswages$exper < 0] <- NA
# factor(): label the 0/1 indicator columns
uswages$race <- factor(uswages$race, labels = c("White", "Black"))
uswages$smsa <- factor(uswages$smsa, labels = c("No", "Yes"))
uswages$pt   <- factor(uswages$pt,   labels = c("No", "Yes"))
# The four region columns hold exactly one 1 per row.
# max.col(): finds which column holds the 1, giving one region factor instead of four dummies
regions <- c("ne", "mw", "so", "we")
uswages$region <- factor(regions[max.col(uswages[, regions], ties.method = "first")],
                         levels = regions)
# Drop the four dummy columns; %in% tests which column names are in the regions vector
uswages <- uswages[, !(names(uswages) %in% regions)]

print(summary(uswages))
# table(): how many workers report each number of years of education
print(table(uswages$educ))
# tapply(..., median): median weekly wage in each group
print(round(tapply(uswages$wage, uswages$pt, median)))       # part time vs full time
print(round(tapply(uswages$wage, uswages$race, median)))     # White vs Black
print(round(tapply(uswages$wage, uswages$region, median)))   # by region

par(mfrow = c(2, 2))
# log(): wage is very skewed, so it is plotted on the log scale
hist(log(uswages$wage), main = "log(weekly wage)", xlab = "log(wage)", col = "lightblue")
# boxplot(): one box of log(wage) for each number of years of education
boxplot(log(wage) ~ educ, data = uswages, main = "log(wage) by education",
        xlab = "years of education", ylab = "log(wage)")
# plot(): scatter of log(wage) vs experience; pch = "." draws small dots because there are 2000 points
plot(log(wage) ~ exper, data = uswages, pch = ".",
     main = "log(wage) vs experience", xlab = "years of experience")
# !is.na(): keeps the rows where exper is not missing
ok <- !is.na(uswages$exper)
# lowess(): smooth trend through the cloud; lines() draws it in red
lines(lowess(uswages$exper[ok], log(uswages$wage[ok])), col = "red", lwd = 2)
boxplot(log(wage) ~ pt, data = uswages, main = "log(wage) by part-time",
        xlab = "part time", ylab = "log(wage)", col = c("lightblue", "orange"))
par(mfrow = c(1, 1))

# What stands out (uswages):
# - wage is strongly right-skewed: median $522, mean $608, maximum $7716 a
#   week (4 workers above $3000). log(wage) is roughly symmetric, so the log
#   scale is the natural one to work on.
# - 33 workers had negative experience, which is impossible; they are set to NA.
#   8 report 0 years of education, which is possible but worth checking.
# - Education piles up at 12 years (719 workers, high school), 16 (280, college)
#   and 18 (197). 18 is the maximum, so it may mean "18 or more".
# - Wages rise with education: median about $475 at 12 years vs $678 at 16.
# - Wages rise with experience up to about 20-30 years (median $707), then fall
#   (median $499 beyond 40 years). The relationship is curved, not straight.
# - Part-time workers earn far less per week (median $148 vs $558), which is
#   expected because this is a weekly, not hourly, wage.
# - Black workers (7.8% of the sample) have a lower median wage ($398 vs $522).
#   The south has the lowest regional median ($475); the others are $522-542.


# =============================================================================
# Exercise 3: prostate - 97 men about to receive a radical prostatectomy
#   lcavol   log(cancer volume)            lweight  log(prostate weight)
#   age      age in years                  lbph     log(benign prostatic hyperplasia)
#   svi      seminal vesicle invasion 0/1  lcp      log(capsular penetration)
#   gleason  Gleason score                 pgg45    % of Gleason scores 4 or 5
#   lpsa     log(prostate specific antigen), the usual response variable
# =============================================================================
data(prostate)

# Correlation of every variable with lpsa (svi is still 0/1 here, so cor() can use it)
print(round(cor(prostate)[, "lpsa"], 2))

prostate$svi <- factor(prostate$svi, labels = c("No", "Yes"))
print(summary(prostate))
print(table(prostate$gleason))   # how many men have each Gleason score

# par(mfrow = c(2, 3)): a 2 x 3 grid for six plots
par(mfrow = c(2, 3))
hist(prostate$lpsa, main = "lpsa", xlab = "log(PSA)", col = "lightblue")
plot(lpsa ~ lcavol, data = prostate, main = "lpsa vs log cancer volume")
boxplot(lpsa ~ svi, data = prostate, main = "lpsa by seminal vesicle invasion",
        xlab = "svi", col = c("lightblue", "orange"))
# factor(gleason) makes boxplot() draw one box per score
boxplot(lpsa ~ factor(gleason), data = prostate, main = "lpsa by Gleason score",
        xlab = "gleason")
# breaks = 20 uses narrow bins so the spike at the lower limit is visible
hist(prostate$lbph, breaks = 20, main = "lbph: spike at log(0.25)",
     xlab = "lbph", col = "lightgreen")
# boxplot() of a single variable: shows the one very large lweight value as a lone point
boxplot(prostate$lweight, main = "lweight: one extreme value", ylab = "lweight")
par(mfrow = c(1, 1))

# What stands out (prostate):
# - 97 men aged 41-79 (median 65). Most variables are already on the log scale
#   (names starting with "l").
# - lpsa is roughly symmetric, from -0.43 to 5.58.
# - lbph and lcp have a pile of values at exactly -1.386 = log(0.25): 43 men
#   for lbph and 45 for lcp. Values below 0.25 seem to have been recorded as
#   0.25, so these variables have a spike at their minimum instead of a smooth
#   shape. pgg45 is similar, with 35 men at 0.
# - gleason only takes 6, 7, 8 or 9, and 8 appears just once, so it behaves
#   more like a category than a continuous measurement.
# - 21 of the 97 men (22%) have seminal vesicle invasion, and their lpsa is
#   clearly higher.
# - lpsa is most strongly related to lcavol (r = 0.73), then svi (0.57) and
#   lcp (0.55). Age hardly matters (r = 0.17).
# - Man 32 has lweight 6.11 when the next largest is 4.78: a prostate about
#   3.8 times heavier than any other. It could be a recording error and is
#   worth checking before modelling.


# =============================================================================
# Exercise 4: sat - school spending and SAT scores, 50 US states, 1994-95
#   expend  spending per pupil, thousands of dollars
#   ratio   pupils per teacher
#   salary  average teacher salary, thousands of dollars
#   takers  % of eligible students who took the SAT
#   verbal, math, total  average scores (total = verbal + math)
# =============================================================================
data(sat)

print(summary(sat))
print(round(cor(sat), 2))
# order(): sorts the states by total score; head() / tail() show the 3 lowest and 3 highest
print(head(sat[order(sat$total), c("expend", "takers", "total")], 3))
print(tail(sat[order(sat$total), c("expend", "takers", "total")], 3))

# Split the states into those where few (< 20%) and many students take the SAT
lowtake <- sat$takers < 20
# ifelse(): picks "low" or "high" for each state; tapply() then averages total within each group
print(round(tapply(sat$total, ifelse(lowtake, "low takers", "high takers"), mean)))

par(mfrow = c(2, 2))
# breaks = 10 asks for about 10 bins so the two clusters of takers show up
hist(sat$takers, breaks = 10, main = "% of students taking the SAT",
     xlab = "takers (%)", col = "lightblue")
plot(total ~ expend, data = sat, main = "SAT total vs spending",
     xlab = "spending per pupil ($1000s)")
# lm(): fits a straight line of total on expend; abline() draws it in red
abline(lm(total ~ expend, data = sat), col = "red", lwd = 2)
plot(total ~ takers, data = sat, main = "SAT total vs % taking the test",
     xlab = "takers (%)")
# Same spending plot again, coloured by takers group: blue = few takers, red = many
plot(total ~ expend, data = sat, pch = 19, col = ifelse(lowtake, "blue", "red"),
     main = "Spending vs total, by takers group", xlab = "spending per pupil ($1000s)")
legend("topright", legend = c("takers < 20%", "takers >= 20%"),
       col = c("blue", "red"), pch = 19, cex = 0.8)
par(mfrow = c(1, 1))

# What stands out (sat):
# - One row per state. total is exactly verbal + math, and both have r = 0.99
#   with total, so total alone is enough to summarise the scores.
# - takers falls into two groups: 22 states where fewer than 20% of students
#   take the SAT (where most students take the ACT instead) and 28 states
#   where 21-81% do.
# - Surprise: states that spend more get LOWER scores (r = -0.38 for expend,
#   -0.44 for teacher salary).
# - The explanation is takers (r = -0.89 with total). Where few students take
#   the SAT, only the strongest ones do, so the state average is high (mean
#   total 1038 vs 910). Higher-spending states also tend to have more takers
#   (r = 0.59), which drags their averages down.
# - Within each takers group, spending and score are unrelated (r = 0.03 and
#   -0.02). Any model of total on expend has to include takers.
# - Lowest totals: South Carolina 844, Georgia 854, North Carolina 865 (58-65%
#   takers). Highest: North Dakota 1107, Iowa 1099, Minnesota 1085 (5-9% takers).


# =============================================================================
# Exercise 5: divusa - divorce in the USA, one row per year 1920-1996
#   divorce     divorces per 1000 married women aged 15+
#   unemployed  unemployment rate (%)
#   femlab      female labour-force participation (%)
#   marriage    marriages per 1000 unmarried women aged 16+
#   birth       births per 1000 women aged 15-44
#   military    military personnel per 1000 population
# =============================================================================
data(divusa)

print(summary(divusa[, -1]))                  # [, -1] leaves out the year column
print(round(cor(divusa)[, "divorce"], 2))     # correlation of each variable with divorce

# This is a time series, so every variable is plotted against year
par(mfrow = c(2, 3))
# for loop: one plot per variable; names(divusa)[-1] is every column except year
for (v in names(divusa)[-1]) {
  # divusa[[v]] picks the column whose name is stored in v; type = "l" joins the points into a line
  plot(divusa$year, divusa[[v]], type = "l", lwd = 2, main = v, xlab = "year", ylab = v)
  # abline(v = ...): dashed vertical lines at 1933 (Great Depression) and 1945 (end of WWII)
  abline(v = c(1933, 1945), lty = 2, col = "grey")
}
par(mfrow = c(1, 1))

# What stands out (divusa):
# - 77 consecutive years, so this is a time series: plot against year, and
#   remember that neighbouring years are not independent observations.
# - divorce rose from 8.0 in 1920 to a peak of 22.8 in 1979, then eased to
#   19.5 by 1996. The fastest rise was 1965-1979 (10.6 to 22.8). There is a
#   short spike after WWII: 17.9 in 1946, back to 10.3 by 1950.
# - History is easy to see: unemployment peaks at 24.9% in 1933 (Great
#   Depression), military jumps to 86.6 per 1000 in 1945 (WWII), marriage
#   spikes to 118.1 in 1946 as troops returned, and births peak at 122.9 in
#   1957 (baby boom) before dropping to about 65 from the mid-1970s.
# - The early 1930s show up as dips in divorce, marriage and birth too:
#   people married, divorced and had children less during the Depression.
# - femlab is a perfectly straight line from 1920 to 1940, so those years were
#   probably filled in between census values rather than measured each year.
# - femlab climbs steadily from 22.7% to 59.3%. Its correlation with divorce
#   is high (0.91), but both mostly trend upward with time (divorce vs year
#   r = 0.88), so the correlation alone doesn't show that one causes the other.
# - divorce moves opposite to birth (r = -0.72) and marriage (r = -0.53), and
#   is barely related to unemployment (-0.21) or military (0.02).
