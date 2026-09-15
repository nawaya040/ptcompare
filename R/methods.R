#' Extract local differences on the representative tree
#' @param x Fitted comparison.
#' @param prob Minimum local posterior difference probability. This is an
#'   exploratory filter, not a p-value or a frequentist error guarantee.
#' @param min_n Minimum pooled region size.
#' @param ... Reserved; must be empty.
#' @return A data frame of internal nodes, ordered by probability and effect.
#'   Regions may overlap through ancestor/descendant nesting. Effects and 95%
#'   intervals condition on the representative tree. A positive effect means
#'   group 1 has greater odds of entering the left child, conditional on the node.
#' @export
regions <- function(x, prob = 0.5, min_n = 0L, ...) {
  .pt_fit(x); .pt_dots(...); .pt_number(prob, "prob", 0, 1)
  .pt_integer(min_n, "min_n", 0, .Machine$integer.max)
  a <- x$tree
  a <- a[a$left > 0 & !is.na(a$prob_difference) & a$prob_difference >= prob & a$n1+a$n2 >= min_n, ]
  cols <- c("node_id", "parent_id", "depth", "condition", "variable", "threshold",
            "n1", "n2", "prob_difference", "effect", "effect_abs", "lower", "upper", "prob_positive")
  a <- a[order(-a$prob_difference, -a$effect_abs), cols, drop = FALSE]
  rownames(a) <- NULL; a
}

#' Computational diagnostics
#' @param x Fitted comparison.
#' @param ... Reserved; must be empty.
#' @return A list with the ESS trajectory, resampling counts, particle weights,
#'   ancestry, stopping diagnostics and timing. Payload bytes measure retained
#'   particle structures at iteration boundaries, excluding allocator overhead,
#'   temporary allocations, input data and R output; this is not peak process RSS.
#' @export
diagnostics <- function(x, ...) { .pt_fit(x); .pt_dots(...); x$diagnostics }

#' @export
print.ptcompare_fit <- function(x, ...) {
  .pt_dots(...)
  cat("Bayesian two-sample comparison (ptcompare)\n")
  cat(sprintf("%s: %d; %s: %d; variables: %d\n", x$group_names[1], x$sample_sizes[1],
              x$group_names[2], x$sample_sizes[2], length(x$variables)))
  cat(sprintf("Posterior P(distributions differ): %.6f\n", x$prob_difference))
  cat("Common coordinate transformation:", x$transform, "\n")
  cat(sprintf("Particles: %d; elapsed: %.2f s; representative internal nodes: %d\n",
              x$control$particles, x$diagnostics$elapsed_seconds, sum(x$tree$left > 0)))
  cat("Local summaries condition on the representative tree. Use summary(), regions(), plot().\n")
  invisible(x)
}

#' @export
summary.ptcompare_fit <- function(object, ..., top = 6L) {
  .pt_dots(...); .pt_integer(top, "top", 0, 10000)
  structure(list(prob_null = object$prob_null, prob_difference = object$prob_difference,
                 sample_sizes = object$sample_sizes, group_names = object$group_names,
                 n_variables = length(object$variables), prior = object$prior,
                 transform = object$transform,
                 diagnostics = object$diagnostics,
                 regions = utils::head(regions(object, prob = 0), top)),
            class = "summary.ptcompare_fit")
}

#' @export
print.summary.ptcompare_fit <- function(x, ...) {
  .pt_dots(...)
  cat("Bayesian two-sample comparison\n")
  cat(sprintf("N = %d / %d; variables = %d\n", x$sample_sizes[1], x$sample_sizes[2], x$n_variables))
  cat(sprintf("P(equal) = %.6g; P(different) = %.6g (finite model)\n", x$prob_null, x$prob_difference))
  cat("Common coordinate transformation:", x$transform, "\n")
  tr <- x$diagnostics$trace
  cat(sprintf("Final ESS = %.1f; resampling steps = %d; elapsed = %.2f s\n",
              x$diagnostics$final_ess, x$diagnostics$resamples, x$diagnostics$elapsed_seconds))
  if (nrow(tr)) cat(sprintf("Final root ancestors = %d; unique particle trees = %d\n",
                            tail(tr$root_ancestors, 1), tail(tr$unique_trees, 1)))
  if (x$diagnostics$depth_limited_leaves > 0)
    cat(sprintf("Depth-limited leaves with >= min_node observations: %d (%.1f%% of observations)\n",
                x$diagnostics$depth_limited_leaves, 100*x$diagnostics$depth_limited_fraction))
  cat("Local split differences, conditional on the representative tree:\n")
  a <- x$regions[, c("node_id", "variable", "threshold", "n1", "n2", "prob_difference", "effect", "lower", "upper"), drop = FALSE]
  print(a, row.names = FALSE, digits = 4)
  cat("Positive effect: group 1 has higher conditional odds of the left child.\n")
  cat("Intervals include the zero point mass. Probabilities are not p-values.\n")
  invisible(x)
}

.pt_fit <- function(x) {
  if (!inherits(x, "ptcompare_fit")) stop("'x' must be a ptcompare_fit object.", call. = FALSE)
}
.pt_dots <- function(...) {
  if (length(list(...))) stop("Unused arguments in '...'.", call. = FALSE)
}

#' Plot a comparison tree, local region or computational diagnostics
#' @param x Fitted comparison.
#' @param type One of 'tree', 'region', or 'diagnostics'.
#' @param node Node ID for a local region plot. Defaults to the node with largest
#'   local difference probability.
#' @param max_depth Display tree through this depth. Hidden descendants are
#'   indicated with a plus sign; the fitted model is unchanged.
#' @param ... Reserved; must be empty.
#' @return Invisibly returns the plotted tree rows or region data, or diagnostics.
#' @export
plot.ptcompare_fit <- function(x, type = c("tree", "region", "diagnostics"),
                              node = NULL, max_depth = 3L, ...) {
  .pt_dots(...); type <- match.arg(type)
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
  if (type == "region") return(.pt_plot_region(x, node))
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
  graphics::plot(NA, xlim = c(0, 1), ylim = c(-max(t$depth[shown])-.65, .5),
                 xlab = "", ylab = "", axes = FALSE,
                 main = "Representative partition tree",
                 sub = "Blue: group 1 favors left; brown: group 2 favors left. + = hidden descendants")
  for (j in shown[-1L]) {
    pa <- t$parent_id[j]; graphics::segments(pos[pa], -t$depth[pa], pos[j], -t$depth[j], col = "grey70")
  }
  graphics::points(pos[shown], -t$depth[shown], pch = 21, bg = color, col = "grey35",
                   cex = 1.5 + 3*sqrt((t$n1+t$n2)[shown]/sum(x$sample_sizes)))
  labs <- vapply(shown, function(j) {
    if (t$left[j] == 0) return(sprintf("#%d leaf\nn=%d/%d", j, t$n1[j], t$n2[j]))
    paste0("#", j, " ", t$variable[j], "\n<= ", formatC(t$threshold[j], digits = 3, format = "g"),
           if (t$depth[j] == max_depth) " +" else "", "\nP(diff)=",
           formatC(t$prob_difference[j], digits = 2, format = "f"))
  }, character(1))
  graphics::text(pos[shown], -t$depth[shown]-.27, labels = labs,
                 cex = if (max_depth <= 3) .7 else .48)
  graphics::mtext("Blue: group 1 favors left; brown: group 2 favors left. + = hidden descendants",
                  side = 1, line = 1, cex = .7)
  invisible(t[shown, ])
}
.pt_plot_region <- function(x, node) {
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
  mids <- head(breaks, -1) + diff(breaks)/2
  old <- graphics::par(mar = c(4, 4, 5, 1)); on.exit(graphics::par(old))
  graphics::matplot(mids, counts, type = "s", lty = 1, lwd = 2,
                    col = c("#176B87", "#B65D35"), xlab = t$variable[node],
                    ylab = "Within-region proportion per bin", ylim = c(0, max(counts)*1.15),
                    main = sprintf("Region #%d: P(local difference) = %.3f", node, t$prob_difference[node]))
  graphics::abline(v = t$threshold[node], lty = 2)
  graphics::mtext(paste(strwrap(t$condition[node], width = 95), collapse = "\n"), side = 3, line = .3, cex = .7)
  graphics::legend("topright", paste0(x$group_names, " (n=", tabulate(groups+1L, 2), ")"),
                   col = c("#176B87", "#B65D35"), lwd = 2, bty = "n")
  invisible(data.frame(value = values, group = x$group_names[groups+1L]))
}
