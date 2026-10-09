#' Regions where the two groups differ
#'
#' Lists the divisions of the representative tree and how much the two groups
#' differ at each one. Use it to find where, and in which direction, the
#' distributions differ.
#'
#' Each row is one division of a region into two parts. Inside the region
#' described by `condition`, the observations are divided into those with
#' `variable <= threshold` (the left part) and the rest (the right part). The
#' row tells whether the two groups fall into the left part in different
#' proportions.
#'
#' @param x Result of [compare_pt()].
#' @param prob Show only rows with `prob_difference` at least this value
#'   (between 0 and 1). Use `prob = 0` to show all divisions.
#' @param min_n Show only regions with at least this many observations (both
#'   groups together).
#' @param ... Not used.
#' @return A data frame with one row per division, sorted by
#'   `prob_difference` (largest first). Columns:
#'   \describe{
#'     \item{`node_id`, `parent_id`, `depth`}{Position of the region in the
#'       tree (as in `plot(fit)`); depth 0 is the whole data space.}
#'     \item{`condition`}{The region, as conditions on the variables.}
#'     \item{`variable`, `threshold`}{The division: the left part is
#'       `variable <= threshold`.}
#'     \item{`n1`, `n2`}{Number of observations of group 1 and group 2 in the
#'       region.}
#'     \item{`prob_difference`}{Posterior probability that the two groups
#'       fall into the left part in different proportions.}
#'     \item{`effect`}{Size and direction of the difference: the log odds
#'       ratio of falling into the left part, group 1 versus group 2 (posterior
#'       mean). Positive values mean group 1 is more likely to be in the left
#'       part; negative values mean group 2 is.}
#'     \item{`effect_abs`}{Posterior mean of the absolute log odds ratio: the
#'       size of the difference regardless of direction.}
#'     \item{`lower`, `upper`}{95% credible interval for the log odds ratio.
#'       A limit is exactly 0 when the groups may not differ at this division.}
#'     \item{`prob_positive`}{Posterior probability that group 1 is more
#'       likely than group 2 to be in the left part.}
#'   }
#' @seealso [compare_pt()], `plot(fit, type = "region", node = ...)` to see the
#'   data of one region.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- matrix(rnorm(200, mean = 0.8), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32))
#' regions(fit)
#' regions(fit, prob = 0.9, min_n = 30)
#' @export
regions <- function(x, prob = 0.5, min_n = 0L, ...) {
  .pt_fit(x); .pt_dots(...); .pt_number(prob, "prob", 0, 1)
  .pt_integer(min_n, "min_n", 0, .Machine$integer.max)
  a <- x$tree
  keep <- a$left > 0 & !is.na(a$prob_difference) & a$prob_difference >= prob & a$n1+a$n2 >= min_n
  .pt_region_table(x, which(keep))
}
# Effect quadrature is computed only for the requested nodes.
.pt_region_table <- function(x, rows) {
  a <- x$tree[rows, , drop = FALSE]
  a <- cbind(a, .pt_effects(a, x$control$effect_draws))
  cols <- c("node_id", "parent_id", "depth", "condition", "variable", "threshold",
            "n1", "n2", "prob_difference", "effect", "effect_abs", "lower", "upper", "prob_positive")
  a <- a[order(-a$prob_difference, -a$effect_abs), cols, drop = FALSE]
  rownames(a) <- NULL; a
}

#' Computational diagnostics
#'
#' Information about the computation of [compare_pt()], such as the number of
#' effective particles and the running time. `plot(fit, type = "diagnostics")`
#' shows the main quantities as graphs.
#'
#' @param x Result of [compare_pt()].
#' @param ... Not used.
#' @return A list with components:
#'   \describe{
#'     \item{`final_ess`}{Effective number of particles at the end (at most
#'       `particles`).}
#'     \item{`resamples`}{Number of resampling steps.}
#'     \item{`elapsed_seconds`}{Running time of [compare_pt()] in seconds.}
#'     \item{`trace`}{Data frame with one row per step of the computation:
#'       effective number of particles before resampling (`ess`), number of
#'       different trees (`unique_trees`), number of distinct initial
#'       particles that the current particles descend from
#'       (`root_ancestors`), and whether resampling took place (`resampled`).}
#'     \item{`weights`, `particle_null`, `particle_scores`}{For each final
#'       particle: its weight, its probability that the groups are the same,
#'       and the fit score used to choose the representative tree.}
#'     \item{`log_evidence`}{Estimated log marginal likelihood of the data.}
#'     \item{`depth_limited_leaves`, `depth_limited_fraction`}{Number of
#'       regions of the representative tree that reached `max_depth` with at
#'       least `min_node` observations, and the fraction of observations in
#'       them.}
#'     \item{`numerical_stop_leaves`}{Number of regions too narrow to divide
#'       further.}
#'     \item{`unique_proposals`, `peak_particle_payload_bytes`}{Number of
#'       distinct division steps computed, and the approximate memory used by
#'       the particles (bytes).}
#'   }
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- matrix(rnorm(200, mean = 0.8), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32))
#' d <- diagnostics(fit)
#' d$final_ess
#' d$elapsed_seconds
#' @export
diagnostics <- function(x, ...) { .pt_fit(x); .pt_dots(...); x$diagnostics }

#' Summarize a comparison
#'
#' `print(fit)` shows the posterior probability that the two distributions
#' differ. `summary(fit)` also shows the settings, a short report of the
#' computation and the divisions with the largest `prob_difference` (see
#' [regions()] for the meaning of the columns).
#'
#' @param x,object Result of [compare_pt()] (or of `summary()`).
#' @param top Number of divisions to show.
#' @param ... Not used.
#' @return `summary()` returns a `summary.ptcompare_fit` object, which is
#'   printed. Its component `regions` holds the divisions shown. `print()`
#'   returns its argument invisibly.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- matrix(rnorm(200, mean = 0.8), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32))
#' fit
#' summary(fit)
#' summary(fit, top = 3)
#' @name summary.ptcompare_fit
NULL

#' @rdname summary.ptcompare_fit
#' @export
print.ptcompare_fit <- function(x, ...) {
  cat("Bayesian two-sample comparison (ptcompare)\n")
  cat(sprintf("%s: %d; %s: %d; variables: %d\n", x$group_names[1], x$sample_sizes[1],
              x$group_names[2], x$sample_sizes[2], length(x$variables)))
  cat(sprintf("Posterior P(distributions differ): %.6f\n", x$prob_difference))
  cat(sprintf("Common coordinate transformation: %s\n", x$transform))
  if (!any(x$tree$left > 0)) cat("No division was possible (fewer observations than min_node).\n")
  cat(sprintf("Particles: %d; elapsed: %.2f s; divisions in the representative tree: %d\n",
              x$control$particles, x$diagnostics$elapsed_seconds, sum(x$tree$left > 0)))
  cat("Use summary(), regions() and plot() for details.\n")
  invisible(x)
}

#' @rdname summary.ptcompare_fit
#' @export
summary.ptcompare_fit <- function(object, ..., top = 6L) {
  .pt_integer(top, "top", 0, 10000)
  t <- object$tree; rows <- which(t$left > 0)
  rows <- rows[order(-t$prob_difference[rows])]
  # Include every node tied with the last shown probability before ordering by effect.
  if (length(rows) > top)
    rows <- rows[t$prob_difference[rows] >= if (top > 0) t$prob_difference[rows[top]] else Inf]
  structure(list(prob_null = object$prob_null, prob_difference = object$prob_difference,
                 sample_sizes = object$sample_sizes, group_names = object$group_names,
                 n_variables = length(object$variables), prior = object$prior,
                 transform = object$transform,
                 diagnostics = object$diagnostics,
                 regions = utils::head(.pt_region_table(object, rows), top)),
            class = "summary.ptcompare_fit")
}

#' @rdname summary.ptcompare_fit
#' @export
print.summary.ptcompare_fit <- function(x, ...) {
  cat("Bayesian two-sample comparison\n")
  cat(sprintf("N = %d / %d; variables = %d\n", x$sample_sizes[1], x$sample_sizes[2], x$n_variables))
  cat(sprintf("P(equal) = %.6g; P(different) = %.6g\n", x$prob_null, x$prob_difference))
  cat(sprintf("Common coordinate transformation: %s\n", x$transform))
  tr <- x$diagnostics$trace
  cat(sprintf("Final ESS = %.1f; resampling steps = %d; elapsed = %.2f s\n",
              x$diagnostics$final_ess, x$diagnostics$resamples, x$diagnostics$elapsed_seconds))
  if (nrow(tr)) cat(sprintf("Final root ancestors = %d; unique particle trees = %d\n",
                            tail(tr$root_ancestors, 1), tail(tr$unique_trees, 1)))
  if (x$diagnostics$depth_limited_leaves > 0)
    cat(sprintf("Regions that reached max_depth: %d (%.1f%% of observations)\n",
                x$diagnostics$depth_limited_leaves, 100*x$diagnostics$depth_limited_fraction))
  cat("Divisions with the largest prob_difference (see ?regions):\n")
  a <- x$regions[, c("node_id", "variable", "threshold", "n1", "n2", "prob_difference", "effect", "lower", "upper"), drop = FALSE]
  print(a, row.names = FALSE, digits = 4)
  cat("effect > 0: group 1 is more likely to satisfy variable <= threshold within the region.\n")
  invisible(x)
}

.pt_fit <- function(x) {
  if (!inherits(x, "ptcompare_fit")) stop("'x' must be a ptcompare_fit object.", call. = FALSE)
}
.pt_dots <- function(...) {
  if (length(list(...))) stop("Unused arguments in '...'.", call. = FALSE)
}

#' Plot the results of a comparison
#'
#' Three kinds of plots:
#' * `type = "tree"` (default): the representative tree. Each node is a
#'   region, labeled with its division (`variable <= threshold`) and
#'   `prob_difference`. Larger nodes contain more observations. The color
#'   shows the direction of the difference: blue when group 1 is more likely
#'   to be in the left part, brown when group 2 is.
#' * `type = "region"`: the distributions of the two groups inside one region,
#'   along the variable used to divide it. The dashed line is the threshold.
#' * `type = "diagnostics"`: the effective number of particles and the number
#'   of different trees at each step of the computation.
#'
#' @param x Result of [compare_pt()].
#' @param type `"tree"`, `"region"` or `"diagnostics"`.
#' @param node For `type = "region"`: the `node_id` of the region to show (see
#'   [regions()]). By default, the region with the largest `prob_difference`.
#' @param max_depth For `type = "tree"`: show the tree down to this depth (0 to
#'   6). A `+` marks nodes whose subregions are not shown.
#' @param ... Graphical parameters such as `main` (not used for
#'   `type = "diagnostics"`).
#' @return The data shown, invisibly.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- matrix(rnorm(200, mean = 0.8), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32))
#' plot(fit)
#' plot(fit, type = "region")
#' plot(fit, type = "diagnostics")
#' @export
plot.ptcompare_fit <- function(x, type = c("tree", "region", "diagnostics"),
                              node = NULL, max_depth = 3L, ...) {
  type <- match.arg(type)
  if (type == "diagnostics") {
    tr <- x$diagnostics$trace
    if (!nrow(tr)) { graphics::plot.new(); graphics::title("No splits: uniform terminal model"); return(invisible(tr)) }
    old <- graphics::par(mfrow = c(2, 1), mar = c(3.5, 4, 2, 1)); on.exit(graphics::par(old))
    graphics::plot(tr$step, tr$ess, type = "l", ylim = c(0, x$control$particles),
                   xlab = "SMC step", ylab = "ESS", main = "ESS before resampling")
    graphics::abline(h = x$control$particles*x$control$ess_threshold, lty = 2, col = "grey50")
    graphics::points(tr$step[tr$resampled == 1], tr$ess[tr$resampled == 1], col = "firebrick", pch = 16)
    graphics::matplot(tr$step, cbind(tr$root_ancestors, tr$unique_trees), type = "l", lty = 1:2,
                      col = c("#176B87", "#B65D35"), xlab = "SMC step", ylab = "Count",
                      main = "Particle ancestry and diversity")
    graphics::legend("topright", c("Root ancestors", "Unique trees"), lty = 1:2,
                     col = c("#176B87", "#B65D35"), bty = "n")
    return(invisible(tr))
  }
  if (type == "region") return(.pt_plot_region(x, node, ...))
  .pt_integer(max_depth, "max_depth", 0, 6)
  t <- x$tree; shown <- which(t$depth <= max_depth); pos <- numeric(nrow(t)); pos[1] <- .5
  for (j in shown[-1L]) {
    pa <- t$parent_id[j]; offset <- 2^(-t$depth[j]-1)
    pos[j] <- pos[pa] + if (t$left[pa] == j) -offset else offset
  }
  pal <- grDevices::colorRampPalette(c("#B65D35", "#F5F4EF", "#176B87"))(101)
  eff <- t$effect[shown]; eff[is.na(eff)] <- 0; limit <- max(abs(eff), .1)
  color <- pal[1 + round(100*(eff/limit + 1)/2)]
  old <- graphics::par(mar = c(3, 1, 4, 1)); on.exit(graphics::par(old))
  .pt_call(graphics::plot, list(NA, xlim = c(0, 1), ylim = c(-max(t$depth[shown])-.65, .5),
                 xlab = "", ylab = "", axes = FALSE,
                 main = "Representative partition tree"), list(...))
  for (j in shown[-1L]) {
    pa <- t$parent_id[j]; graphics::segments(pos[pa], -t$depth[pa], pos[j], -t$depth[j], col = "grey70")
  }
  graphics::points(pos[shown], -t$depth[shown], pch = 21, bg = color, col = "grey35",
                   cex = 1.5 + 3*sqrt((t$n1+t$n2)[shown]/sum(x$sample_sizes)))
  labs <- vapply(shown, function(j) {
    if (t$left[j] == 0) return(sprintf("#%d leaf\nn=%d/%d", j, t$n1[j], t$n2[j]))
    paste0("#", j, " ", t$variable[j], "\n<= ", formatC(t$threshold[j], digits = 4, format = "g"),
           if (t$depth[j] == max_depth) " +" else "", "\nP(diff)=",
           formatC(t$prob_difference[j], digits = 2, format = "f"))
  }, character(1))
  graphics::text(pos[shown], -t$depth[shown]-.27, labels = labs,
                 cex = if (max_depth <= 3) .7 else .48)
  graphics::mtext("Blue: group 1 favors left; brown: group 2 favors left. + = hidden descendants",
                  side = 1, line = 1, cex = .7)
  invisible(t[shown, ])
}
.pt_plot_region <- function(x, node, ...) {
  if (is.null(x$data)) stop("Region plots require keep_data = TRUE.", call. = FALSE)
  t <- x$tree
  if (is.null(node)) {
    r <- regions(x, prob = 0)
    if (!nrow(r)) stop("The fitted tree has no internal regions.", call. = FALSE)
    node <- r$node_id[1]
  }
  .pt_integer(node, "node", 1, nrow(t))
  if (t$left[node] == 0) stop("Choose an internal node for a split comparison.", call. = FALSE)
  selected <- rep(TRUE, nrow(x$data)); j <- node
  while (t$parent_id[j] > 0) {
    pa <- t$parent_id[j]; v <- x$data[, t$variable_id[pa]]
    selected <- selected & if (t$left[pa] == j) v <= t$threshold[pa] else v > t$threshold[pa]
    j <- pa
  }
  values <- x$data[selected, t$variable_id[node]]; groups <- x$group[selected]
  if (!length(values)) stop("This region has no observations.", call. = FALSE)
  breaks <- pretty(range(c(values, t$threshold[node])), n = 20)
  if (length(breaks) < 2) breaks <- c(min(values)-.5, max(values)+.5)
  counts <- sapply(0:1, function(g) {
    if (!any(groups == g)) return(rep(0, length(breaks)-1))
    graphics::hist(values[groups == g], breaks = breaks, plot = FALSE)$counts / sum(groups == g)
  })
  old <- graphics::par(mar = c(4, 4, 5, 1)); on.exit(graphics::par(old))
  # Step function over bin edges: each proportion spans its own bin.
  .pt_call(graphics::matplot, list(breaks, rbind(counts, counts[nrow(counts), ]), type = "s",
                    lty = 1, lwd = 2, col = c("#176B87", "#B65D35"), xlab = t$variable[node],
                    ylab = "Within-region proportion per bin", ylim = c(0, max(counts)*1.15),
                    main = sprintf("Region #%d: P(local difference) = %.3f", node, t$prob_difference[node])),
           list(...))
  graphics::abline(v = t$threshold[node], lty = 2)
  graphics::mtext(paste(strwrap(t$condition[node], width = 95), collapse = "\n"), side = 3, line = .3, cex = .7)
  graphics::legend("topright", paste0(x$group_names, " (n=", tabulate(groups+1L, 2), ")"),
                   col = c("#176B87", "#B65D35"), lwd = 2, bty = "n")
  invisible(data.frame(value = values, group = x$group_names[groups+1L]))
}

# Call a plotting function with defaults overridden by user arguments.
.pt_call <- function(f, defaults, user) {
  drop <- names(defaults) != "" & names(defaults) %in% names(user)
  do.call(f, c(defaults[!drop], user))
}

#' @export
print.pt_prior <- function(x, ...) {
  cat("ptcompare model settings (pt_prior)\n")
  cat(sprintf("  %-15s %s\n", names(x), vapply(x, format, character(1))), sep = "")
  invisible(x)
}

#' @export
print.pt_control <- function(x, ...) {
  cat("ptcompare computational settings (pt_control)\n")
  cat(sprintf("  %-15s %s\n", names(x), vapply(x, format, character(1))), sep = "")
  invisible(x)
}
