#' Settings for sequential Monte Carlo
#' @param particles Number of particles. Increase this to assess Monte Carlo stability.
#' @param ess_threshold Resample when ESS is below this fraction of particles.
#' @param keep_particles Keep every final particle tree (can use substantial memory).
#' @param keep_data Store original data for region plots.
#' @param effect_draws Deterministic quadrature size for absolute effects and intervals.
#' @return A validated list of computational settings.
#' @export
pt_control <- function(particles = 256L, ess_threshold = 0.5,
                       keep_particles = FALSE, keep_data = TRUE,
                       effect_draws = 4096L) {
  .pt_integer(particles, "particles", 2, 100000)
  .pt_number(ess_threshold, "ess_threshold", 0, 1)
  .pt_flag(keep_particles, "keep_particles")
  .pt_flag(keep_data, "keep_data")
  .pt_integer(effect_draws, "effect_draws", 256, 1048576)
  structure(list(particles = as.integer(particles), ess_threshold = ess_threshold,
                 keep_particles = keep_particles, keep_data = keep_data,
                 effect_draws = as.integer(effect_draws)), class = "pt_control")
}

#' Finite hidden Markov Polya tree specification
#' @param max_depth Maximum split depth; the root has depth zero.
#' @param split_grid Odd number of equally spaced relative split candidates.
#' @param min_node Nodes with fewer pooled observations are uniform leaves.
#' @param gamma Probability parameter for entering the difference state.
#' @param beta Depth decay of persistence in the difference state (at least one).
#' @param prune Transition probability into absorbing sharing.
#' @param precision Total beta precision (between 1e-6 and 1e6), centered at the geometric split fraction.
#' @param eta Concentration of the split prior around the midpoint. Candidate
#'   weights are proportional to exp(-eta * pooled_node_size * abs(fraction-0.5)).
#' @param midpoint_spike Use the spike-and-slab location prior of Awaya and Ma
#'   (2024, Appendix A). The midpoint receives its prior mass from the split
#'   prior above; once a node is divided at the midpoint, all its descendants
#'   are divided at their midpoints. Set to FALSE for the location prior (2)
#'   at every node.
#' @return A validated model specification. These settings change the model,
#'   including the finite resolution and pooled-data stopping rule.
#' @export
pt_prior <- function(max_depth = 8L, split_grid = 15L, min_node = 20L,
                     gamma = 0.3, beta = 1, prune = 0.3,
                     precision = 1, eta = 0.1, midpoint_spike = TRUE) {
  .pt_integer(max_depth, "max_depth", 1, 20)
  .pt_integer(split_grid, "split_grid", 1, 127)
  if (split_grid %% 2 != 1) stop("'split_grid' must be odd (to include the midpoint).", call. = FALSE)
  .pt_integer(min_node, "min_node", 2, .Machine$integer.max)
  .pt_number(gamma, "gamma", 0, 1, open = TRUE)
  .pt_number(prune, "prune", 0, 1, open = TRUE)
  .pt_number(beta, "beta", 1, 100)
  .pt_number(precision, "precision", 1e-6, 1e6)
  .pt_number(eta, "eta", 0, 100)
  .pt_flag(midpoint_spike, "midpoint_spike")
  structure(list(max_depth = as.integer(max_depth), split_grid = as.integer(split_grid),
                 min_node = as.integer(min_node), gamma = gamma, beta = beta,
                 prune = prune, precision = precision, eta = eta,
                 midpoint_spike = midpoint_spike), class = "pt_prior")
}

#' Compare two multivariate distributions
#'
#' Fits a finite adaptive hidden Markov Polya tree comparison. Global evidence
#' is averaged over SMC trees. Local summaries are conditional on the sampled
#' tree with the largest prior-times-likelihood score, not a guaranteed global
#' maximizer. Observations must be independent within and across the two groups.
#'
#' @param x Numeric matrix or numeric data frame for group 1, or pooled data.
#' @param y Numeric matrix or data frame for group 2. Use either y or group.
#' @param group Two-level labels for rows of pooled x. Factor levels set group
#'   order; other labels use order of first appearance.
#' @param support Optional matrix with one row per original variable and two
#'   columns (lower, upper). For range coordinates, omitted support uses pooled
#'   ranges padded by 1 percent.
#' @param transform Common coordinate transformation. 'range' uses linear
#'   scaling. 'rank' maps pooled midranks to [0,1], preserves ties and stores
#'   an inverse interpolation for explanations. It changes the partition prior
#'   on the original scale and conditions on an estimated marginal transform.
#'   Explicit support is supported only with 'range'.
#'   The default is 'rank' when support is omitted and 'range' when supplied.
#' @param prior Model settings from [pt_prior()].
#' @param control Computational settings from [pt_control()].
#' @return An object of class `ptcompare_fit`, including the representative tree,
#'   global posterior probability of equality, preprocessing and diagnostics.
#' @details Terminal regions have common uniform conditional distributions.
#'   Global equality concerns this finite model. Data-derived support and node
#'   stopping are conditioned on the pooled data. No group-specific scaling is
#'   applied. Constant columns are removed with a warning. Missing and non-finite
#'   observations are rejected. Ties go to the left child at a split.
#' @references Awaya, N. and Ma, L. (2024). Hidden Markov Polya Trees for
#'   High-Dimensional Distributions. JASA 119(545), 189--201.
#'   doi:10.1080/01621459.2022.2105223.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(120), ncol = 2)
#' y <- matrix(rnorm(120, 0.6), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32),
#'                   prior = pt_prior(max_depth = 3, split_grid = 3))
#' summary(fit)
#' regions(fit, prob = 0)
#' @export
compare_pt <- function(x, y = NULL, group = NULL, support = NULL,
                       prior = pt_prior(), control = pt_control(),
                       transform = if (is.null(support)) "rank" else "range") {
  start <- proc.time()[[3L]]
  transform <- match.arg(transform, c("rank", "range"))
  # Revalidate lists so editing a settings object cannot bypass checks.
  if (!inherits(prior, "pt_prior") || !inherits(control, "pt_control"))
    stop("Use pt_prior() and pt_control() to specify settings.", call. = FALSE)
  prior <- do.call(pt_prior, unclass(prior))
  control <- do.call(pt_control, unclass(control))
  input <- .pt_prepare(x, y, group, support, transform)
  raw <- pt_smc_cpp(input$unit, input$group, unclass(prior), unclass(control))
  tree <- raw$tree
  internal <- tree$left > 0L
  tree$variable <- NA_character_
  tree$variable[internal] <- input$variables[tree$variable_id[internal]]
  tree$threshold <- NA_real_
  j <- tree$variable_id[internal]
  tree$threshold[internal] <- vapply(which(internal), function(i) {
    j <- tree$variable_id[i]
    if (transform == "range") input$support[j, 1] + tree$cut_unit[i] *
      (input$support[j, 2] - input$support[j, 1])
    else stats::approx(input$knots[[j]]$unit, input$knots[[j]]$original,
                       xout = tree$cut_unit[i], ties = "ordered")$y
  }, numeric(1))
  tree$condition <- .pt_conditions(tree)
  effect <- .pt_effects(tree, control$effect_draws)
  tree <- cbind(tree, effect)
  stopped <- tree$left == 0 & tree$depth == prior$max_depth &
    tree$n1 + tree$n2 >= prior$min_node
  out <- list(call = match.call(), prob_null = raw$prob_null,
              prob_difference = 1 - raw$prob_null, tree = tree,
              group_names = input$group_names, variables = input$variables,
              sample_sizes = tabulate(input$group + 1L, 2L), prior = prior,
              control = control, support = input$support,
              dropped_variables = input$dropped,
              support_source = input$support_source,
              transform = transform, transform_knots = input$knots,
              data = if (control$keep_data) input$original else NULL,
              group = if (control$keep_data) input$group else NULL,
              diagnostics = list(trace = raw$trace, weights = raw$weights,
                particle_null = raw$particle_null, particle_scores = raw$particle_scores,
                final_ess = 1 / sum(raw$weights^2),
                resamples = raw$resamples, unique_proposals = raw$proposals,
                peak_particle_payload_bytes = raw$peak_particle_payload_bytes,
                log_evidence = raw$log_evidence,
                depth_limited_leaves = sum(stopped),
                numerical_stop_leaves = sum(tree$stop_reason == "numerical precision"),
                depth_limited_fraction = sum((tree$n1 + tree$n2)[stopped]) / sum(input$group >= 0),
                elapsed_seconds = proc.time()[[3L]] - start),
              particles = raw$particles)
  class(out) <- "ptcompare_fit"
  out
}

.pt_number <- function(x, name, lower, upper, open = FALSE) {
  valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid) valid <- if (open) x > lower && x < upper else x >= lower && x <= upper
  if (!valid) stop(sprintf("'%s' must be a finite number %s %s and %s.",
                          name, if (open) "strictly between" else "between", lower, upper), call. = FALSE)
}
.pt_integer <- function(x, name, lower, upper) {
  .pt_number(x, name, lower, upper)
  if (x != floor(x)) stop(sprintf("'%s' must be an integer.", name), call. = FALSE)
}
.pt_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x))
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
}
.pt_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, function(z) is.numeric(z) && !is.complex(z), logical(1))))
      stop(sprintf("'%s' must contain only numeric columns.", name), call. = FALSE)
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || nrow(x) < 2L || ncol(x) < 1L)
    stop(sprintf("'%s' must be a numeric matrix/data frame with at least two rows.", name), call. = FALSE)
  if (any(!is.finite(x))) stop(sprintf("'%s' contains missing or non-finite values.", name), call. = FALSE)
  storage.mode(x) <- "double"
  x
}
.pt_prepare <- function(x, y, group, support, transform = "range") {
  x <- .pt_matrix(x, "x")
  if (is.null(y) == is.null(group)) stop("Supply exactly one of 'y' and 'group'.", call. = FALSE)
  if (!is.null(y)) {
    y <- .pt_matrix(y, "y")
    if (ncol(x) != ncol(y)) stop("'x' and 'y' must have the same columns.", call. = FALSE)
    if (!is.null(colnames(x)) && !is.null(colnames(y)) && !identical(colnames(x), colnames(y)))
      stop("Column names and order must match between groups.", call. = FALSE)
    g <- c(rep.int(0L, nrow(x)), rep.int(1L, nrow(y)))
    labels <- c("group 1", "group 2")
    x <- rbind(x, y)
  } else {
    if (!is.atomic(group) || !is.null(dim(group)) || length(group) != nrow(x) || anyNA(group))
      stop("'group' must be a nonmissing vector with one label per row.", call. = FALSE)
    labels <- if (is.factor(group)) levels(droplevels(group)) else unique(as.character(group))
    if (length(labels) != 2L) stop("'group' must contain exactly two groups.", call. = FALSE)
    g <- match(as.character(group), labels) - 1L
    if (any(tabulate(g + 1L, 2L) < 2L)) stop("Each group needs at least two observations.", call. = FALSE)
  }
  if (is.null(colnames(x))) colnames(x) <- paste0("X", seq_len(ncol(x)))
  if (anyNA(colnames(x)) || any(!nzchar(colnames(x))) || anyDuplicated(colnames(x)))
    stop("Variable names must be nonempty and unique.", call. = FALSE)
  ranges <- t(apply(x, 2, range))
  keep <- ranges[, 1] < ranges[, 2]
  if (!any(keep)) stop("All variables are constant in the pooled data.", call. = FALSE)
  dropped <- colnames(x)[!keep]
  if (length(dropped)) warning("Removed pooled constant variables: ", paste(dropped, collapse = ", "), call. = FALSE)
  source <- if (is.null(support)) "pooled padded range" else "user supplied"
  if (transform == "rank") {
    if (!is.null(support)) stop("Explicit 'support' is available only with transform = 'range'.", call. = FALSE)
    support <- ranges; source <- "pooled rank transform"
  } else if (is.null(support)) {
    widths <- ranges[, 2] - ranges[, 1]
    support <- cbind(ranges[, 1] - .01 * widths, ranges[, 2] + .01 * widths)
  } else {
    if (!is.matrix(support) || !is.numeric(support) || is.complex(support) ||
        !identical(dim(support), c(ncol(x), 2L)) || any(!is.finite(support)) ||
        any(support[, 1] >= support[, 2]))
      stop("'support' must have finite lower < upper bounds, one row per original variable.", call. = FALSE)
    if (any(ranges[, 1] < support[, 1] | ranges[, 2] > support[, 2]))
      stop("All observations must lie inside 'support'.", call. = FALSE)
  }
  x <- x[, keep, drop = FALSE]; support <- support[keep, , drop = FALSE]
  width <- support[, 2] - support[, 1]
  if (any(!is.finite(support)) || any(!is.finite(width)) || any(width <= 0))
    stop("Data scale is too extreme for safe range scaling; rescale the input.", call. = FALSE)
  unit <- sweep(sweep(x, 2, support[, 1], "-"), 2, width, "/")
  knots <- NULL
  if (transform == "rank") {
    knots <- vector("list", ncol(x))
    for (j in seq_len(ncol(x))) {
      r <- rank(x[, j], ties.method = "average")
      unit[, j] <- (r - min(r)) / (max(r) - min(r))
      ord <- order(x[, j]); ord <- ord[!duplicated(x[ord, j])]
      knots[[j]] <- list(unit = unit[ord, j], original = x[ord, j])
    }
  }
  list(unit = unit, original = x, support = support, group = as.integer(g),
       variables = colnames(x), group_names = labels, dropped = dropped,
       support_source = source, knots = knots)
}

.pt_conditions <- function(tree) {
  conditions <- rep("all observations", nrow(tree))
  for (j in seq_len(nrow(tree))[-1L]) {
    par <- tree$parent_id[j]
    clause <- paste0(tree$variable[par], if (tree$left[par] == j) " <= " else " > ",
                     format(tree$threshold[par], digits = 6, trim = TRUE))
    conditions[j] <- if (par == 1L) clause else paste(conditions[par], clause, sep = " & ")
  }
  conditions
}

# Deterministic low-discrepancy integration: does not consume the user's RNG.
.pt_radical <- function(n, base) {
  k <- seq_len(n); out <- numeric(n); f <- 1 / base
  while (any(k > 0)) { out <- out + f * (k %% base); k <- k %/% base; f <- f / base }
  out
}
.pt_effects <- function(tree, n) {
  ans <- data.frame(effect_abs = rep(NA_real_, nrow(tree)),
                    lower = NA_real_, upper = NA_real_, prob_positive = NA_real_)
  # Log-gamma construction avoids beta draws rounded to 0/1 near boundaries.
  u <- lapply(c(2, 3, 5, 7, 11, 13, 17, 19), function(b) .pt_radical(n, b))
  loggamma <- function(shape, j) {
    if (shape < 1) log(stats::qgamma(u[[j]], shape + 1)) + log(u[[j + 1]]) / shape
    else log(stats::qgamma(u[[j]], shape))
  }
  for (j in which(tree$left > 0)) {
    delta <- loggamma(tree$a1[j], 1) - loggamma(tree$b1[j], 3) -
      loggamma(tree$a2[j], 5) + loggamma(tree$b2[j], 7)
    # Symmetrize quadrature across the two beta distributions for group exchange.
    other <- loggamma(tree$a2[j], 1) - loggamma(tree$b2[j], 3) -
      loggamma(tree$a1[j], 5) + loggamma(tree$b1[j], 7)
    delta <- c(delta, -other)
    p <- tree$prob_difference[j]
    ans$effect_abs[j] <- p * mean(abs(delta))
    ans$prob_positive[j] <- p * mean(delta > 0)
    vals <- c(delta, 0); w <- c(rep(p / length(delta), length(delta)), 1 - p)
    ord <- order(vals); cw <- cumsum(w[ord]); cw[length(cw)] <- 1
    ans$lower[j] <- vals[ord[which(cw >= .025)[1]]]
    ans$upper[j] <- vals[ord[which(cw >= .975)[1]]]
  }
  ans
}
