# RStudio: open ptcompare.Rproj, then Source this file.
#
# Provisional example (added 2026-10-27 on user request): it shows how the
# output and the plots look with one variable and with three variables. It has
# not been through the review the other examples had, and it may be revised,
# renamed or dropped.
source("scripts/load-local.R")

# ---- One variable -----------------------------------------------------------
# Both groups have mean 0 and variance 1; group 2 is a two-component mixture.
# A plain numeric vector is accepted and becomes a single column named X1.
set.seed(101)
n <- 800L
s <- sqrt(1 - 0.9^2)
x1 <- rnorm(n)
y1 <- ifelse(runif(n) < 0.5, rnorm(n, -0.9, s), rnorm(n, 0.9, s))

set.seed(102)
fit1 <- compare_pt(x1, y1)
print(summary(fit1))
print(regions(fit1, prob = 0.9, min_n = 30))

# With one variable every node divides the same variable, so the tree is a
# sequence of thresholds, and "regions" draws the distributions of the two
# groups inside each region instead of a scatter plot.
plot(fit1)
plot(fit1, type = "regions")
# plot(fit1, type = "levels") is an error: "levels" needs exactly two variables.

# Independent samples from the same distribution.
set.seed(103)
null_fit1 <- compare_pt(rnorm(600), rnorm(600))
print(null_fit1$prob_difference)

# ---- Three variables --------------------------------------------------------
# Group 2 differs only where X < 0 and Y > 0, where Z is shifted upwards.
set.seed(301)
n <- 900L
draw <- function(n) cbind(X = runif(n, -1, 1), Y = runif(n, -1, 1), Z = rnorm(n))
x3 <- draw(n)
y3 <- draw(n)
affected <- y3[, "X"] < 0 & y3[, "Y"] > 0
y3[affected, "Z"] <- y3[affected, "Z"] + 1.6

set.seed(302)
fit3 <- compare_pt(x3, y3)
print(summary(fit3))
region_table3 <- regions(fit3, prob = 0.9, min_n = 40)
print(region_table3)

# Each panel draws two of the three variables, so its box is the region seen in
# those two only. Observations that the third variable puts outside the region
# stay gray inside the box; the symbols say what is in the region.
plot(fit3, max_depth = 3)
plot(fit3, type = "regions")

# The division that was built into the data need not appear on the same
# variable: here the shift in Z is reported as a difference in the conditional
# probability of Y within regions of Z and X.

# Identical N(0,1) margins in all three variables; only the correlations differ.
set.seed(303)
x3b <- matrix(rnorm(3*n), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
L <- chol(matrix(c(1, .7, .7,  .7, 1, .7,  .7, .7, 1), 3))
y3b <- matrix(rnorm(3*n), ncol = 3) %*% L
colnames(y3b) <- c("A", "B", "C")

set.seed(304)
dependence_fit3 <- compare_pt(x3b, y3b)
print(summary(dependence_fit3))
plot(dependence_fit3, type = "regions")
