#' Computational settings
#'
#' Settings that control how [compare_pt()] carries out its computation.
#' They affect accuracy and running time. The defaults are suitable for most
#' analyses.
#'
#' @param particles Number of candidate trees computed in parallel (the
#'   particles of sequential Monte Carlo). More particles give more accurate
#'   results but take longer. Try `particles = 1024` or more for a final
#'   analysis.
#' @param ess_threshold Fraction of `particles` (between 0 and 1). When the
#'   effective number of particles falls below this fraction, the particles
#'   are resampled.
#' @param resample_power Number in (0, 1] used when resampling. Particles are
#'   drawn with probabilities proportional to their weights raised to this
#'   power. The default 0.5 follows Awaya and Ma (2024); 1 gives standard
#'   resampling.
#' @param keep_particles If `TRUE`, all final candidate trees are kept in the
#'   result (component `particles`). This can use a lot of memory.
#' @param keep_data If `TRUE`, the data are kept in the result. They are
#'   needed for `plot(fit, type = "regions")` and `plot(fit, type = "levels")`.
#' @param effect_draws Number of points (256 to 1048576) used to compute
#'   `effect_abs`, `lower`, `upper` and `prob_positive` in [regions()].
#'   Larger values are more accurate and slower.
#' @return A `pt_control` object to pass to [compare_pt()].
#' @seealso [pt_prior()] for the statistical model settings.
#' @examples
#' pt_control()
#' pt_control(particles = 1024)
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

#' Model settings
#'
#' Settings of the statistical model used by [compare_pt()]: how finely the
#' data space may be divided into regions, where the divisions may be placed,
#' and how likely the two groups are to differ before seeing the data.
#'
#' @param max_depth Maximum number of successive divisions (1 to 20). Larger
#'   values allow smaller regions and take longer.
#' @param split_grid Number of candidate division points for a variable within
#'   a region (an odd number from 1 to 127). The candidates are equally spaced
#'   and include the middle of the region.
#' @param min_node Regions with fewer observations than this (both groups
#'   together) are not divided further.
#' @param gamma Prior tendency of the two groups to differ, between 0 and 1.
#'   Larger values make differences more likely before seeing the data.
#' @param beta Number from 1 to 100. Values above 1 make differences in smaller
#'   (deeper) regions less likely before seeing the data. The default 1
#'   follows Awaya and Ma (2024).
#' @param prune Prior probability, between 0 and 1, that a region and all its
#'   subregions are the same in the two groups.
#' @param precision Positive number. Larger values make the proportion of
#'   observations falling on each side of a division closer, before seeing the
#'   data, to the proportion of the region's width on that side.
#' @param eta Non-negative number. Larger values favor division points closer
#'   to the middle of a region, more strongly in regions with many
#'   observations. `eta = 0` treats all candidates equally.
#' @param midpoint_spike If `TRUE` (default), once a region is divided at its
#'   middle, all its subregions are also divided at their middles, as in
#'   Awaya and Ma (2024, Appendix A). If `FALSE`, every candidate division
#'   point is allowed in every region.
#' @return A `pt_prior` object to pass to [compare_pt()].
#' @seealso [pt_control()] for the computational settings.
#' @examples
#' pt_prior()
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

#' Compare the distributions of two groups
#'
#' Compares the distributions of two groups of observations on one or more
#' numeric variables, using the hidden Markov Pólya tree method of Awaya and
#' Ma (2024). It gives the posterior probability that the two distributions
#' differ, and finds regions of the data space where they differ. Each region
#' is described by conditions on the variables, such as `age <= 40`.
#'
#' The method divides the data space into regions step by step, building a
#' tree of divisions, and compares how the two groups fall on the two sides of
#' each division. Many candidate trees are explored; the one that best fits
#' the data (the representative tree) is used to describe where the groups
#' differ.
#'
#' @param x Data for group 1: a numeric matrix or data frame with one row per
#'   observation and one column per variable, or a numeric vector for a single
#'   variable. When `group` is used, `x` contains the data of both groups.
#' @param y Data for group 2, with the same columns as `x`. Give either `y` or
#'   `group`.
#' @param group Group labels for the rows of `x`, with exactly two different
#'   values. For a factor, the first level is group 1. Otherwise, the label
#'   that appears first is group 1.
#' @param support Optional known range of the variables: a matrix with one row
#'   per variable and two columns (lower and upper limit), or a vector
#'   `c(lower, upper)` for a single variable. Used with `transform = "range"`.
#' @param transform How the variables are put on a common scale before
#'   dividing.
#'   * `"rank"` (default when `support` is not given): each variable is
#'     replaced by the rank of its values in both groups combined. Results do
#'     not change if a variable is transformed monotonically (for example by
#'     `log()`).
#'   * `"range"` (default when `support` is given): each variable is rescaled
#'     linearly to the interval from 0 to 1, as in Awaya and Ma (2024). Without
#'     `support`, the observed range is extended by 1% on each side.
#'
#'   In both cases, thresholds and regions are reported in the original units.
#' @param prior Model settings from [pt_prior()].
#' @param control Computational settings from [pt_control()].
#' @return A `ptcompare_fit` object. Use [summary()][summary.ptcompare_fit],
#'   [regions()] and [plot()][plot.ptcompare_fit] to look at the results. Its
#'   main components are:
#'   \describe{
#'     \item{`prob_difference`}{Posterior probability that the two
#'       distributions differ.}
#'     \item{`prob_null`}{Posterior probability that they are the same
#'       (`1 - prob_difference`).}
#'     \item{`tree`}{The representative tree, one row per region. [regions()]
#'       gives an easier-to-read table of its divisions.}
#'     \item{`diagnostics`}{Information about the computation; see
#'       [diagnostics()].}
#'     \item{`group_names`, `sample_sizes`, `variables`}{Group names, number of
#'       observations per group and names of the variables used.}
#'     \item{`prior`, `control`, `transform`}{The settings used.}
#'   }
#' @details Variables that take the same value in all observations are
#'   removed with a warning. Missing (`NA`) and infinite values are not
#'   allowed; remove those rows first. If there are fewer than `min_node`
#'   observations in total (see [pt_prior()]), no division is possible and
#'   `prob_difference` is 0.
#' @references Awaya, N. and Ma, L. (2024). Hidden Markov Pólya Trees for
#'   High-Dimensional Distributions. *Journal of the American Statistical
#'   Association*, 119(545), 189--201. \doi{10.1080/01621459.2022.2105223}
#' @seealso [regions()], [summary()][summary.ptcompare_fit],
#'   [plot()][plot.ptcompare_fit], [pt_prior()], [pt_control()]
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- matrix(rnorm(200, mean = 0.8), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32))
#' fit
#' summary(fit)
#'
#' # Data in one data frame with a group column
#' two <- droplevels(subset(iris, Species != "setosa"))
#' fit2 <- compare_pt(two[, 1:4], group = two$Species,
#'                    control = pt_control(particles = 32))
#' regions(fit2, prob = 0.9)
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
    warning("No split was possible: there are fewer observations than 'min_node' ",
            "(see ?pt_prior). prob_difference is 0.", call. = FALSE)
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
