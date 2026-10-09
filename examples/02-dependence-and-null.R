source("scripts/load-local.R")

# Identical N(0,1) marginal distributions; only the correlation differs.
set.seed(61)
n <- 700L
x <- cbind(A = rnorm(n), B = rnorm(n))
z <- rnorm(n)
y <- cbind(A = z, B = 0.8*z + sqrt(1-0.8^2)*rnorm(n))
set.seed(62)
dependence_fit <- compare_pt(x, y, transform = "range")
print(summary(dependence_fit))
plot(dependence_fit)
plot(dependence_fit, type = "regions")

# Independent samples from the same distribution.
set.seed(63)
null_x <- matrix(rnorm(2*n), ncol = 2, dimnames = list(NULL, c("A", "B")))
null_y <- matrix(rnorm(2*n), ncol = 2, dimnames = list(NULL, c("A", "B")))
set.seed(64)
null_fit <- compare_pt(null_x, null_y, transform = "range")
print(summary(null_fit))
plot(null_fit)

# P(different) is a posterior probability under the finite model, not a p-value.
