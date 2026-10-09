#' Settings for sequential Monte Carlo
#' @param particles Number of particles. Increase this to assess Monte Carlo stability.
#' @param ess_threshold Resample when ESS is below this fraction of particles.
#' @param resample_power Power kappa in (0, 1] of the resampling probabilities.
#'   Ancestors are drawn with probabilities proportional to W^kappa and receive
#'   weights proportional to W / W^kappa (Awaya and Ma, 2024, Section 3.1).
#'   The default 0.5 is the paper's choice; 1 gives standard resampling.
#' @param keep_particles Keep every final particle tree (can use substantial memory).
#' @param keep_data Store original data for region plots.
#' @param effect_draws Deterministic quadrature size for absolute effects and
#'   intervals (between 256 and 1048576).
#' @return A validated list of computational settings.
#' @examples
#' pt_control(particles = 512)
#' @export
pt_control <- function(particles = 256L, ess_threshold = 0.5,
                       keep_particles = FALSE, keep_data = TRUE,
                       effect_draws = 4096L, resample_power = 0.5) {
  .pt_integer(particles, "particles", 2, 100000)
  .pt_number(ess_threshold, "ess_threshold", 0, 1)
  .pt_number(resample_power, "resample_power", 0, 1)
  if (resample_power == 0) stop("'resample_power' must be greater than 0.", call. = FALSE)
  .pt_flag(keep_particles, "keep_particles")
  .pt_flag(keep_data, "keep_data")
  .pt_integer(effect_draws, "effect_draws", 256, 1048576)
  structure(list(particles = as.integer(particles), ess_threshold = ess_threshold,
                 keep_particles = keep_particles, keep_data = keep_data,
                 effect_draws = as.integer(effect_draws),
                 resample_power = resample_power), class = "pt_control")
}

#' Finite hidden Markov Polya tree specification
#' @param max_depth Maximum split depth (1 to 20); the root has depth zero.
#' @param split_grid Odd number (1 to 127) of equally spaced relative split
#'   candidates k/(split_grid+1).
#' @param min_node Nodes with fewer pooled observations are uniform leaves.
#' @param gamma Probability parameter for entering the difference state, in (0, 1).
#' @param beta Depth decay of persistence in the difference state (1 to 100);
#'   1 gives the transition matrix of Awaya and Ma (2024, Supplement F.2.1).
#' @param prune Transition probability into absorbing sharing, in (0, 1).
#' @param precision Total beta precision (between 1e-6 and 1e6), centered at the geometric split fraction.
#' @param eta Concentration (0 to 100) of the split prior around the midpoint. Candidate
#'   weights are proportional to exp(-eta * pooled_node_size * abs(fraction-0.5)).
#' @param midpoint_spike Use the spike-and-slab location prior of Awaya and Ma
#'   (2024, Appendix A). The midpoint receives its prior mass from the split
#'   prior above; once a node is divided at the midpoint, all its descendants
#'   are divided at their midpoints. Set to FALSE for the location prior (2)
#'   at every node.
#' @return A validated model specification. These settings change the model,
#'   including the finite resolution and pooled-data stopping rule.
#' @examples
#' pt_prior(max_depth = 6, split_grid = 7)
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
#' @param x Numeric matrix, numeric data frame or numeric vector (one variable)
#'   for group 1, or pooled data.
#' @param y Numeric matrix, data frame or vector for group 2. Use either y or group.
#' @param group Two-level labels for rows of pooled x. Factor levels set group
#'   order; other labels use order of first appearance (so `c(1, 1, 0, 0)` makes
#'   "1" group 1). The sign of `effect` refers to this order.
#' @param support Optional matrix with one row per original variable and two
#'   columns (lower, upper); a length-two vector is accepted for one variable.
#'   For range coordinates, omitted support uses pooled ranges padded by 1 percent.
#' @param transform Common coordinate transformation. 'range' uses linear
#'   scaling. 'rank' maps pooled midranks to the unit interval, preserves ties and stores
#'   an inverse interpolation for explanations. It changes the partition prior
#'   on the original scale and conditions on an estimated marginal transform.
#'   Explicit support is supported only with 'range'.
#'   'range' is the linear rescaling to the unit cube of Awaya and Ma (2024).
#'   The default is 'rank' when support is omitted and 'range' when supplied.
#' @param prior Model settings from [pt_prior()].
#' @param control Computational settings from [pt_control()].
#' @return An object of class `ptcompare_fit`, a list with components
#' \describe{
#'   \item{prob_null, prob_difference}{Posterior probability that the two
#'     distributions are equal, and its complement, under the finite model.}
#'   \item{tree}{Representative tree, one row per node: `node_id`, `parent_id`,
#'     `left`, `right` (child IDs, 0 for leaves), `depth`, `variable_id`,
#'     `variable`, `threshold` (original units), `cut_unit` and `fraction`
#'     (unit-scale cut and relative split location), `n1`, `n2`,
#'     `prob_difference` (local posterior probability of state D), `effect`
#'     (posterior mean signed log odds ratio), `a1`, `b1`, `a2`, `b2` (beta
#'     posterior parameters under D), `stop_reason` and `condition`. Split
#'     columns are NA for leaves. Absolute effects and intervals are computed
#'     by [regions()].}
#'   \item{diagnostics}{See [diagnostics()].}
#'   \item{particles}{With `keep_particles = TRUE`, the final particle trees as
#'     engine tables (unit-scale cuts only, no names or thresholds).}
#'   \item{group_names, variables, sample_sizes, prior, control, support,
#'     dropped_variables, support_source, transform, transform_knots, data,
#'     group, call}{Inputs and preprocessing.}
#' }
#' @details Terminal regions have common uniform conditional distributions.
#'   Global equality concerns this finite model. Data-derived support and node
#'   stopping are conditioned on the pooled data. No group-specific scaling is
#'   applied. Constant columns are removed with a warning. Missing and non-finite
#'   observations are rejected. Ties go to the left child at a split. If the
#'   pooled sample is smaller than `min_node`, no split is possible and the
#'   finite model gives P(equal) = 1; a warning is issued.
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
  tree$threshold[internal] <- vapply(which(internal), function(i) {
    j <- tree$variable_id[i]
    thr <- if (transform == "range") input$support[j, 1] + tree$cut_unit[i] *
      (input$support[j, 2] - input$support[j, 1])
    else stats::approx(input$knots[[j]]$unit, input$knots[[j]]$original,
                       xout = tree$cut_unit[i], ties = "ordered")$y
    .pt_snap(thr, tree$cut_unit[i], input$original[, j], input$unit[, j])
  }, numeric(1))
  tree$condition <- .pt_conditions(tree)
  if (!any(internal))
    warning("No split was possible (pooled sample size below 'min_node'); ",
            "the finite model gives P(equal) = 1.", call. = FALSE)
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
  if (is.null(dim(x)) && is.numeric(x)) x <- matrix(x, ncol = 1L)
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
    if (is.null(dim(support)) && length(support) == 2L) support <- matrix(support, 1L)
    if (!is.matrix(support) || !is.numeric(support) || is.complex(support) ||
        !identical(dim(support), c(ncol(x), 2L)) || any(!is.finite(support)) ||
        any(support[, 1] >= support[, 2]))
      stop("'support' must have finite lower < upper bounds, one row per original variable.", call. = FALSE)
    if (any(ranges[, 1] < support[, 1] | ranges[, 2] > support[, 2]))
      stop("All observations must lie inside 'support'.", call. = FALSE)
  }
  x <- x[, keep, drop = FALSE]; support <- support[keep, , drop = FALSE]
  knots <- NULL
  if (transform == "range") {
    width <- support[, 2] - support[, 1]
    if (any(!is.finite(support)) || any(!is.finite(width)) || any(width <= 0))
      stop("Data scale is too extreme for safe range scaling; rescale the input.", call. = FALSE)
    unit <- sweep(sweep(x, 2, support[, 1], "-"), 2, width, "/")
  } else {
    unit <- x
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

# The engine splits on unit coordinates (unit <= cut goes left). Move the
# original-scale threshold, only if rounding requires it, so that
# original <= threshold reproduces exactly the same assignment.
.pt_snap <- function(thr, cut, original, unit) {
  left <- unit <= cut
  lmax <- if (any(left)) max(original[left]) else -Inf
  rmin <- if (any(!left)) min(original[!left]) else Inf
  if (thr >= lmax && thr < rmin) return(thr)
  if (is.finite(rmin) && is.finite(lmax)) (lmax + rmin) / 2 else lmax
}
.pt_conditions <- function(tree) {
  conditions <- rep("all observations", nrow(tree))
  for (j in seq_len(nrow(tree))[-1L]) {
    par <- tree$parent_id[j]
    clause <- paste0(tree$variable[par], if (tree$left[par] == j) " <= " else " > ",
                     format(tree$threshold[par], digits = 7, trim = TRUE))
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
  na <- rep(NA_real_, nrow(tree))
  ans <- data.frame(effect_abs = na, lower = na, upper = na, prob_positive = na)
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
