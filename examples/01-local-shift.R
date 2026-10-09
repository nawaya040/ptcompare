# RStudio: open ptcompare.Rproj, then Source this file.
# The package has already been built into .local-library on this machine.
source("scripts/load-local.R")

set.seed(20260915)
n <- 1000L
x <- cbind(marker_A = runif(n), marker_B = runif(n))
y <- cbind(marker_A = runif(n), marker_B = runif(n))
# Group 2 differs only in the subpopulation marker_A < 0.35.
affected <- y[, "marker_A"] < 0.35
y[affected, "marker_B"] <- rbeta(sum(affected), 5, 1)

set.seed(2026)
fit <- compare_pt(x, y, support = cbind(c(0, 0), c(1, 1)), transform = "range")
print(summary(fit))
region_table <- regions(fit, prob = 0.9, min_n = 40)
print(region_table)

# The RStudio Plots pane retains previous plots (use its back arrow).
plot(fit, max_depth = 3)
plot(fit, type = "regions")
plot(fit, type = "diagnostics")

# Inspect the generated objects in Environment:
# fit, region_table, x, y, affected.
# An internal node ID from region_table can be used for another region plot:
# plot(fit, type = "regions", node = region_table$node_id[2])
# Local probabilities/intervals condition on the representative tree.
# Positive effect means group 1 has higher odds of entering its left child.
