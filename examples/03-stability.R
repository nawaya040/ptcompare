source("scripts/load-local.R")
set.seed(71)
n <- 600
x <- matrix(rnorm(n*10), ncol=10)
y <- matrix(rnorm(n*10), ncol=10)
y[, 1] <- y[, 1] + 0.45
colnames(x) <- colnames(y) <- paste0("X", 1:10)

# Hold the model fixed; vary only particle count and random seed.
settings <- expand.grid(particles = c(64L, 256L, 1024L), seed = 1:3)
results <- vector("list", nrow(settings))
for (i in seq_len(nrow(settings))) {
  set.seed(settings$seed[i])
  f <- compare_pt(x, y, control = pt_control(particles = settings$particles[i]))
  r <- regions(f, prob = 0)
  results[[i]] <- data.frame(settings[i, ], prob_difference = f$prob_difference,
    seconds = diagnostics(f)$elapsed_seconds, top_variable = r$variable[1],
    top_threshold = r$threshold[1], root_ancestors = tail(diagnostics(f)$trace$root_ancestors, 1))
}
stability_table <- do.call(rbind, results)
print(stability_table)
# Particle count and seed affect only the Monte Carlo approximation. Changing
# max_depth or split_grid in pt_prior() changes the finite model itself.
