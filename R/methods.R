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
#' @param x Result of [compare_pt()], or the table returned by `regions()`
#'   when printing it.
#' @param prob Show only rows with `prob_difference` at least this value
#'   (between 0 and 1). Use `prob = 0` to show all divisions.
#' @param min_n Show only regions with at least this many observations of
#'   each group.
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
#' @seealso [compare_pt()]; `plot(fit, type = "regions")` draws these regions
#'   with the data.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), ncol = 2)
#' y <- matrix(rnorm(200, mean = 0.8), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32))
#' regions(fit)
#' regions(fit, prob = 0.9, min_n = 20)
#' @export
regions <- function(x, prob = 0.5, min_n = 0L, ...) {
  .pt_fit(x); .pt_dots(...); .pt_number(prob, "prob", 0, 1)
  .pt_integer(min_n, "min_n", 0, .Machine$integer.max)
  a <- x$tree
  keep <- a$left > 0 & !is.na(a$prob_difference) & a$prob_difference >= prob &
    a$n1 >= min_n & a$n2 >= min_n
  .pt_region_table(x, which(keep))
}

#' @rdname regions
#' @export
print.ptcompare_regions <- function(x, ...) {
  a <- as.data.frame(x)
  # Rounded for reading only; the returned table keeps its full precision. The
  # thresholds in `condition` are already the shortest that stay exact.
  for (nm in intersect(c("prob_difference", "prob_positive"), names(a)))
    a[[nm]] <- round(a[[nm]], 4)
  for (nm in intersect(c("effect", "effect_abs", "lower", "upper"), names(a)))
    a[[nm]] <- round(a[[nm]], 3)
  if (!is.null(a$threshold)) a$threshold <- signif(a$threshold, 6)
  print(a, ...)
  invisible(x)
}
# Effect quadrature is computed only for the requested nodes.
.pt_region_table <- function(x, rows) {
  a <- x$tree[rows, , drop = FALSE]
  a <- cbind(a, .pt_effects(a, x$control$effect_draws))
  cols <- c("node_id", "parent_id", "depth", "condition", "variable", "threshold",
            "n1", "n2", "prob_difference", "effect", "effect_abs", "lower", "upper", "prob_positive")
  a <- a[order(-a$prob_difference, -a$effect_abs), cols, drop = FALSE]
  rownames(a) <- NULL; class(a) <- c("ptcompare_regions", "data.frame"); a
}

#' Computational diagnostics
#'
#' Information about the computation of [compare_pt()], such as the number of
#' effective particles and the running time. `plot(fit, type = "diagnostics")`
#' shows the effective number of particles as a graph.
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
  cat(sprintf("Final ESS = %.1f; resampling steps = %d; elapsed = %.2f s\n",
              x$diagnostics$final_ess, x$diagnostics$resamples, x$diagnostics$elapsed_seconds))
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
#' Four kinds of plots:
#' * `type = "tree"` (default): the representative tree. Each node is a
#'   region, labeled with its division (`variable <= threshold`) and
#'   `prob_difference`. Larger nodes contain more observations. The color
#'   shows the direction of the difference: blue when group 1 is more likely
#'   to be in the left part, brown when group 2 is.
#' * `type = "regions"`: the regions where the two groups differ most, one
#'   panel per region (similar to Figure 14 of the supplement of Awaya and Ma,
#'   2024). The regions are those of [regions()] with `prob_difference` at
#'   least `prob` and at least `min_n` observations of each group, ordered by
#'   `rank_by`; at most `top` are shown. With two or more variables, each
#'   panel is a scatter plot: the horizontal axis is the variable that divides
#'   the region (dashed line: threshold), the vertical axis another variable
#'   that bounds the region. Observations inside the region are shown as blue
#'   squares (group 1) and brown triangles (group 2), the others in gray.
#'   The box is the region seen in the two variables of the panel only. With
#'   three or more variables the region may also be bounded by variables that
#'   the panel does not show, so an observation can lie inside the box and
#'   still be gray; the symbols, not the box, show which observations are in
#'   the region. With two variables the box is the region itself. With one
#'   variable, each panel shows the distributions of the two groups inside
#'   the region.
#' * `type = "levels"`: for two variables only. One scatter plot per depth of
#'   the tree; each region at that depth is shaded by its `prob_difference`
#'   (as in Figure 7 of Soriano and Ma, 2017). Dotted lines show how each
#'   region is divided at the next depth.
#' * `type = "diagnostics"`: the effective number of particles at each step of
#'   the computation, with the resampling threshold as a dashed line and the
#'   steps where resampling took place marked. The other recorded quantities
#'   are in [diagnostics()].
#'
#' `"regions"` and `"levels"` need the data, which are kept unless
#' `pt_control(keep_data = FALSE)` was used.
#'
#' @param x Result of [compare_pt()].
#' @param type `"tree"`, `"regions"`, `"levels"` or `"diagnostics"`.
#' @param node For `type = "regions"`: `node_id`s of regions to show (see
#'   [regions()]), instead of choosing them by `prob`, `min_n` and `rank_by`.
#'   A single region is drawn as one large panel.
#' @param max_depth For `type = "tree"`: show the tree down to this depth (0 to
#'   6; default 3). A `+` marks nodes whose subregions are not shown. For
#'   `type = "levels"`: the deepest level shown (default: all levels, at most 8).
#' @param top For `type = "regions"`: maximum number of regions to show.
#' @param prob For `type = "regions"`: show only regions with
#'   `prob_difference` at least this value.
#' @param rank_by For `type = "regions"`: which regions to show first.
#'   `"probability"` (default): largest `prob_difference`, the posterior
#'   probability of a difference. `"effect"`: largest `effect_abs`, the size of
#'   the difference (see [regions()]).
#' @param min_n For `type = "regions"`: show only regions with at least this
#'   many observations of each group.
#' @param ... Graphical parameters such as `main` (for `"tree"`, and for
#'   `"regions"` with a single region).
#' @return The data shown, invisibly.
#' @references Soriano, J. and Ma, L. (2017). Probabilistic multi-resolution
#'   scanning for two-sample differences. *Journal of the Royal Statistical
#'   Society, Series B*, 79(2), 547--572. \doi{10.1111/rssb.12180}
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(400), ncol = 2)
#' y <- matrix(rnorm(400, mean = 0.8), ncol = 2)
#' fit <- compare_pt(x, y, control = pt_control(particles = 32))
#' plot(fit)
#' plot(fit, type = "regions")
#' plot(fit, type = "regions", top = 4, rank_by = "effect")
#' plot(fit, type = "regions", node = 2)
#' plot(fit, type = "levels", max_depth = 3)
#' plot(fit, type = "diagnostics")
#' @export
plot.ptcompare_fit <- function(x, type = c("tree", "regions", "levels", "diagnostics"),
                              node = NULL, max_depth = 3L, top = 9L, prob = 0.5,
                              rank_by = c("probability", "effect"), min_n = 50L, ...) {
  type <- match.arg(type)
  if (type == "diagnostics") {
    tr <- x$diagnostics$trace
    if (!nrow(tr)) { graphics::plot.new(); graphics::title("No splits: uniform terminal model"); return(invisible(tr)) }
    graphics::plot(tr$step, tr$ess, type = "l", ylim = c(0, x$control$particles),
                   xlab = "SMC step", ylab = "ESS", main = "ESS before resampling")
    graphics::abline(h = x$control$particles*x$control$ess_threshold, lty = 2, col = "grey50")
    graphics::points(tr$step[tr$resampled == 1], tr$ess[tr$resampled == 1], col = "firebrick", pch = 16)
    return(invisible(tr))
  }
  if (type == "regions")
    return(.pt_plot_regions(x, node, top, prob, match.arg(rank_by), min_n, list(...)))
  if (type == "levels") return(.pt_plot_levels(x, if (missing(max_depth)) NULL else max_depth))
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
# One variable: distributions of the two groups inside the region, as step
# histograms over the bins (proportion of each group's observations per bin).
.pt_hist_panel <- function(x, node, b, main, user = list(), label_cex = 1) {
  t <- x$tree; inside <- .pt_inside(x, b, node)
  values <- x$data[inside, 1]; groups <- x$group[inside]
  breaks <- pretty(range(c(values, t$threshold[node])), n = 20)
  if (length(breaks) < 2) breaks <- c(min(values) - .5, max(values) + .5)
  counts <- sapply(0:1, function(g) {
    if (!any(groups == g)) return(rep(0, length(breaks) - 1))
    graphics::hist(values[groups == g], breaks = breaks, plot = FALSE)$counts / sum(groups == g)
  })
  .pt_call(graphics::matplot, list(breaks, rbind(counts, counts[nrow(counts), ]), type = "s",
                    lty = 1, lwd = 2, col = .pt_cols, xlab = colnames(x$data)[1],
                    ylab = "Proportion per bin", ylim = c(0, max(counts) * 1.1),
                    main = main, cex.lab = label_cex), user)
  graphics::abline(v = t$threshold[node], lty = 2, lwd = 1.5)
  data.frame(node_id = node, x_variable = colnames(x$data)[1], y_variable = NA_character_,
             n1 = t$n1[node], n2 = t$n2[node], prob_difference = t$prob_difference[node])
}

.pt_need_data <- function(x) {
  if (is.null(x$data)) stop("This plot needs the data; use pt_control(keep_data = TRUE).", call. = FALSE)
}
# Scatter-plot symbols: group 1 dark blue squares, group 2 brown triangles.
.pt_cols <- c("#123E7C", "#B65D35")
.pt_pch <- c(0, 2)

# Lower and upper limits (original units) of every node of the representative
# tree, one row per node and one column per variable; -Inf/Inf when unbounded.
.pt_bounds <- function(x) {
  t <- x$tree; d <- ncol(x$data)
  lo <- matrix(-Inf, nrow(t), d); hi <- matrix(Inf, nrow(t), d)
  for (j in seq_len(nrow(t))[-1L]) {
    pa <- t$parent_id[j]; v <- t$variable_id[pa]
    lo[j, ] <- lo[pa, ]; hi[j, ] <- hi[pa, ]
    if (t$left[pa] == j) hi[j, v] <- t$threshold[pa] else lo[j, v] <- t$threshold[pa]
  }
  list(lo = lo, hi = hi)
}
.pt_inside <- function(x, b, node) {
  inside <- rep(TRUE, nrow(x$data))
  for (v in seq_len(ncol(x$data)))
    inside <- inside & x$data[, v] > b$lo[node, v] & x$data[, v] <= b$hi[node, v]
  inside
}
# Vertical variable: the most recent division above the node on another variable.
.pt_partner <- function(t, node, xv, d) {
  j <- node
  while (t$parent_id[j] > 0) {
    j <- t$parent_id[j]
    if (t$variable_id[j] != xv) return(t$variable_id[j])
  }
  setdiff(seq_len(d), xv)[1L]
}
# Draw points of both groups, interleaved so that neither group hides the other.
.pt_points <- function(px, py, g, cex) {
  ord <- order(stats::ave(seq_along(g), g, FUN = seq_along) / tabulate(g + 1L, 2L)[g + 1L])
  graphics::points(px[ord], py[ord], pch = .pt_pch[g[ord] + 1L],
                   col = .pt_cols[g[ord] + 1L], cex = cex, lwd = 1.1)
}
.pt_panel <- function(x, node, b, main, user = list(), cex = .85, label_cex = 1) {
  t <- x$tree; data <- x$data; d <- ncol(data)
  xv <- t$variable_id[node]; yv <- .pt_partner(t, node, xv, d)
  inside <- .pt_inside(x, b, node)
  .pt_call(graphics::plot, list(data[, xv], data[, yv], type = "n", main = main,
                                xlab = colnames(data)[xv], ylab = colnames(data)[yv],
                                cex.lab = label_cex), user)
  graphics::points(data[!inside, xv], data[!inside, yv], pch = 16, cex = cex * .7, col = "grey82")
  .pt_points(data[inside, xv], data[inside, yv], x$group[inside], cex)
  u <- graphics::par("usr")
  xl <- max(b$lo[node, xv], u[1]); xh <- min(b$hi[node, xv], u[2])
  yl <- max(b$lo[node, yv], u[3]); yh <- min(b$hi[node, yv], u[4])
  graphics::rect(xl, yl, xh, yh, lwd = 2.2)
  graphics::segments(t$threshold[node], yl, t$threshold[node], yh, lwd = 2.2, lty = 2)
  data.frame(node_id = node, x_variable = colnames(data)[xv], y_variable = colnames(data)[yv],
             n1 = t$n1[node], n2 = t$n2[node], prob_difference = t$prob_difference[node])
}
.pt_group_legend <- function(x, where, ..., pch = .pt_pch) {
  graphics::legend(where, x$group_names, pch = pch, col = .pt_cols, bty = "n", ...)
}

.pt_plot_regions <- function(x, node, top, prob, rank_by, min_n, user) {
  .pt_need_data(x); t <- x$tree
  if (is.null(node)) {
    .pt_integer(top, "top", 1, 100); .pt_number(prob, "prob", 0, 1)
    .pt_integer(min_n, "min_n", 0, .Machine$integer.max)
    cand <- which(t$left > 0 & t$prob_difference >= prob & t$n1 >= min_n & t$n2 >= min_n)
    if (!length(cand)) {
      message("No region has prob_difference >= ", prob, " and at least ", min_n,
              " observations of each group; try smaller 'prob' or 'min_n'.")
      return(invisible(NULL))
    }
    r <- .pt_region_table(x, cand)
    ord <- if (rank_by == "effect") order(-r$effect_abs, -r$prob_difference) else
      order(-r$prob_difference, -r$effect_abs)
    r <- utils::head(r[ord, , drop = FALSE], top)
  } else {
    if (!is.numeric(node) || anyNA(node) || any(node != floor(node)) ||
        any(node < 1 | node > nrow(t)) || any(t$left[node] == 0))
      stop("'node' must contain node_id values of divided regions (see regions()).", call. = FALSE)
    r <- .pt_region_table(x, node)
    r <- r[match(node, r$node_id), , drop = FALSE]
  }
  k <- nrow(r); one <- k == 1L; d <- ncol(x$data)
  nc <- ceiling(sqrt(k)); nr <- ceiling(k / nc)
  old <- if (one) graphics::par(mar = c(4.2, 4.2, 3.2, 1)) else
    graphics::par(mfrow = c(nr, nc), mar = c(3.2, 3.2, 2.2, .6), mgp = c(2, .6, 0), oma = c(3, 0, 0, 0))
  on.exit(graphics::par(old))
  b <- .pt_bounds(x); out <- vector("list", k)
  for (i in seq_len(k)) {
    j <- r$node_id[i]
    main <- sprintf("#%d  P %.2f  |effect| %.2f", j, r$prob_difference[i], r$effect_abs[i])
    out[[i]] <- if (d == 1L)
      .pt_hist_panel(x, j, b, main, if (one) user else list(), label_cex = if (one) 1 else .95)
    else .pt_panel(x, j, b, main, if (one) user else list(),
                   cex = if (one) .85 else .65, label_cex = if (one) 1 else .95)
  }
  legend_args <- if (d == 1L) list(lty = 1, lwd = 2, pch = NA) else list(pt.lwd = 1.3)
  if (one) {
    do.call(.pt_group_legend, c(list(x, "topright"), legend_args))
  } else {
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot.new()
    do.call(.pt_group_legend, c(list(x, "bottom", horiz = TRUE, xpd = NA, inset = c(0, .003),
                                     cex = 1.4, pt.cex = 1.4), legend_args))
  }
  invisible(cbind(do.call(rbind, out), effect_abs = r$effect_abs))
}

.pt_plot_levels <- function(x, max_depth) {
  .pt_need_data(x)
  if (ncol(x$data) != 2L) stop("type = 'levels' is available for two variables only.", call. = FALSE)
  t <- x$tree
  deepest <- max(t$depth[t$left > 0], 0L)
  if (is.null(max_depth)) max_depth <- min(deepest, 7L)
  .pt_integer(max_depth, "max_depth", 0, 20)
  levels <- 0:min(max_depth, deepest)
  b <- .pt_bounds(x); data <- x$data
  pal <- grDevices::colorRampPalette(c("#FFFFFF", "#D4D4D4", "#A3A3A3"))(101)
  k <- length(levels) + 1L; nc <- ceiling(sqrt(k)); nr <- ceiling(k / nc)
  old <- graphics::par(mfrow = c(nr, nc), mar = c(3.2, 3.2, 2, .6), mgp = c(2, .6, 0))
  on.exit(graphics::par(old))
  for (k in levels) {
    graphics::plot(data[, 1], data[, 2], type = "n", main = sprintf("Level %d", k),
                   xlab = colnames(data)[1], ylab = colnames(data)[2])
    u <- graphics::par("usr")
    for (j in which(t$depth == k)) {
      shade <- if (t$left[j] > 0) pal[1 + round(100 * t$prob_difference[j])] else "#FFFFFF"
      graphics::rect(max(b$lo[j, 1], u[1]), max(b$lo[j, 2], u[3]),
                     min(b$hi[j, 1], u[2]), min(b$hi[j, 2], u[4]),
                     col = shade, border = "grey70", lwd = .6)
    }
    .pt_points(data[, 1], data[, 2], x$group, .6)
    # Dotted lines: how each region at this level is divided at the next level.
    for (j in which(t$depth == k & t$left > 0)) {
      xl <- max(b$lo[j, 1], u[1]); xh <- min(b$hi[j, 1], u[2])
      yl <- max(b$lo[j, 2], u[3]); yh <- min(b$hi[j, 2], u[4]); thr <- t$threshold[j]
      if (t$variable_id[j] == 1L) graphics::segments(thr, yl, thr, yh, lty = 3, lwd = 1.3)
      else graphics::segments(xl, thr, xh, thr, lty = 3, lwd = 1.3)
    }
    graphics::box()
  }
  # Legend panel: shading scale and group symbols.
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  y <- seq(.15, .85, length.out = 102)
  graphics::rect(.1, y[-102], .25, y[-1], col = pal, border = NA)
  graphics::rect(.1, .15, .25, .85)
  graphics::text(.28, c(.15, .5, .85), c("0", "0.5", "1"), adj = 0, cex = 1.4)
  graphics::text(.05, .95, "P(difference)", adj = 0, cex = 1.5)
  graphics::legend(.42, .62, x$group_names, pch = .pt_pch, col = .pt_cols, bty = "n",
                   cex = 1.4, pt.cex = 1.4, pt.lwd = 1.5)
  invisible(t[t$depth %in% levels, c("node_id", "depth", "variable", "threshold", "prob_difference")])
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
