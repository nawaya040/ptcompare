source("scripts/load-local.R")
set.seed(501)
n <- 500
d <- 20
x <- matrix(rnorm(n*d), ncol=d)
y <- matrix(rnorm(n*d), ncol=d)
y[,1] <- y[,1]+0.6
colnames(x) <- colnames(y) <- paste0("X",seq_len(d))

# Default: pooled rank coordinates. Cutpoints and region conditions are
# reported in the original units.
set.seed(502)
rank_fit <- compare_pt(x,y,transform="rank")
print(summary(rank_fit))
plot(rank_fit)
plot(rank_fit,type="regions")

# Linear rescaling to the unit cube, as in Awaya and Ma (2024).
set.seed(502)
range_fit <- compare_pt(x,y,transform="range")
print(summary(range_fit))
