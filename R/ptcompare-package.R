#' ptcompare: compare the distributions of two groups
#'
#' Compares the distributions of two groups of observations on numeric
#' variables and finds the regions where they differ, using the hidden Markov
#' Pólya tree method of Awaya and Ma (2024).
#'
#' A typical analysis:
#' 1. `fit <- compare_pt(x, y)` compares group 1 (`x`) with group 2 (`y`).
#' 2. `fit` or `summary(fit)` shows the posterior probability that the
#'    distributions differ.
#' 3. `regions(fit, prob = 0.9)` lists the regions where the groups differ.
#' 4. `plot(fit)` and `plot(fit, type = "regions")` show the results.
#'
#' Settings are changed with [pt_prior()] (statistical model) and
#' [pt_control()] (computation).
#' @references Awaya, N. and Ma, L. (2024). Hidden Markov Pólya Trees for
#'   High-Dimensional Distributions. *Journal of the American Statistical
#'   Association*, 119(545), 189--201. \doi{10.1080/01621459.2022.2105223}
#' @useDynLib ptcompare, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils head tail
"_PACKAGE"
