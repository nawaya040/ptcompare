source("scripts/load-local.R")
set.seed(501)
n <- 500
d <- 20
x <- matrix(rnorm(n*d), ncol=d)
y <- matrix(rnorm(n*d), ncol=d)
y[,1] <- y[,1]+0.6
colnames(x) <- colnames(y) <- paste0("X",seq_len(d))

# Pooled rank coordinates reduce competition from unrelated marginal shapes.
# This is a different coordinate prior; it conditions on the estimated transform.
# Cutpoints and region conditions are still reported in the original units.
set.seed(502)
rank_fit <- compare_pt(x,y,transform="rank")
print(summary(rank_fit))
plot(rank_fit)
plot(rank_fit,type="region")

# Hold data fixed and inspect sensitivity to the coordinate prior.
set.seed(502)
range_fit <- compare_pt(x,y,transform="range")
print(data.frame(transform=c("range","rank"),
  prob_difference=c(range_fit$prob_difference,rank_fit$prob_difference)))
# See evaluation/REPORT.md for repeated trials, including null behavior.
