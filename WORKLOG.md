# Implementation progress

Authorized scope: implement PLAN stages through explanation and statistical
evaluation, with runnable RStudio examples. No CRAN submission.

- Environment: R 4.5.2 at C:/Program Files/R/R-4.5.2/bin/Rscript.exe;
  Rtools45, Rcpp, testthat, roxygen2 available. MRS not installed.
- Design: finite adaptive HMPT, uniform leaves, three MRS states. Exact partial
  tree likelihood increments via up/down messages, including sibling evidence.
  This replaces the old path-only updates. Finite regular split grid; no legacy
  absorbing midpoint restriction. Pooled-data stopping/support conditioning is
  documented. Global null constrains internal splits only.
- Storage: shared immutable observation blocks at unsplit nodes, compact nodes,
  shared particle trees after resampling; clone only when a particle diverges.
- Implemented: C++ engine and R API, S3 summaries/plots/regions, deterministic
  log-gamma effect quadrature, original-unit explanations, input validation.
- Exact enumeration tests passed for fixed state trees and all small adaptive
  trees with/without resampling, plus exchange/order/ties/scaling/API tests.
- Range study: 1,220 data sets completed. 2/10 dimensional null n=200 and each
  alternative n=50; 50/100 dimensional null n=40, local/correlation n=25; heavy
  tails and unbalanced sizes n=80 null / 40 alternatives. Comparison with MMD,
  held-out rpart classification, and low-dimensional MRS 1.3.2.
- Finding: range-coordinate finite trees can spend their depth on unrelated
  marginal shapes; particle count alone does not resolve this. Added OPTIONAL
  transform='rank', with pooled midranks, tie-preserving inverse interpolation,
  explicit model-conditioning documentation. After paired validation, the final
  default is rank without support, range with explicitly supplied support.
  Studies and saved examples now specify range explicitly where appropriate,
  so their results are unaffected by this default selection.
- Paired rank study complete: 1,220 paired data sets, 2,440 PT fits. Tables,
  confidence intervals, localization measures and figures saved; conclusions
  and limitations are in evaluation/REPORT.md.
- Final R CMD check --no-manual --no-vignettes passed Status: OK on
  Windows/R4.5.2. Tests: 86 expectations passed, zero failures/warnings/skips.
  Logs and source SHA256 hashes are in evaluation/results.
- All four RStudio examples executed, images visually inspected. Outputs in
  evaluation/figures and results/examples.log. No further example rerun needed
  unless behavior changes.
- Completed isolated Rterm benchmarks through n=200,000, d=20, M=32:
  rank 10.56 seconds, approximately 369 MiB peak working set. Fixed shallow
  budget is a feasibility benchmark, not a guarantee of posterior accuracy.
- Authorized scope complete. Final package installed in .local-library;
  ptcompare.Rproj and four runnable examples are ready. CRAN submission,
  Linux/macOS/R-devel checks and further methodological work remain outside
  this completed development/evaluation stage.
- Build calls need require_escalated: sandbox loses Rtools PATH in subprocesses.
  Normal build env PATH prepend C:/rtools45/usr/bin and
  C:/rtools45/x86_64-w64-mingw32.static.posix/bin, LC_ALL=C. R_MAKEVARS_USER is
  unnecessary outside sandbox. Final .local-library lives inside project.

## 2026-10-09

- Plots: `type = "tree"` (default), `"regions"` (largest differences, scatter
  panels for d >= 2, histograms for d = 1; filters prob >= 0.5 and min_n = 50
  per group, rank_by = "probability"), `"levels"` (d = 2 only), and
  `"diagnostics"`. `type = "region"` was merged into `"regions"` (`node =`).
- Pending (user decision): figures that potential users see first, such as in
  a vignette or the README, should show data with two or more variables
  (scatter-type `regions`/`levels` plots), not one-variable histograms.

## 2026-10-27

- `plot(fit, type = "diagnostics")` now draws the ESS trace only. The second
  panel (root ancestors and unique trees) was removed on user instruction:
  particle degeneracy is worth showing as a weak point of the method in a
  paper, but not in the package, where the plot is read as a check on the
  user's own run. The same two numbers were dropped from the `summary()`
  printout for the same reason, the figure and the text having carried the
  same information; `summary()` still reports final ESS, resampling steps and
  elapsed time. Both series stay in `diagnostics(fit)$trace` and their Rd
  entry is unchanged, so nothing is lost programmatically, and the scripts in
  evaluation/ that read `diagnostics(fit)$trace$root_ancestors` are
  unaffected. The `par(mfrow = c(2, 1), mar = ...)` override went with the
  second panel; the single panel uses the device defaults, which also fixes
  the "SMC step" label being clipped on tall PNG devices.
- `?plot.ptcompare_fit` now states that with three or more variables the box
  of a `"regions"` panel is the region projected onto the two variables drawn,
  so an observation can lie inside the box and still be gray; the symbols, not
  the box, mark the region. This is a description of the existing plot, not a
  change to it.
- Added `examples/05-one-and-three-variables.R` (one variable: equal mean and
  variance, different shape, plus a null; three variables: a difference
  confined to a box in X and Y, and a dependence-only difference). PROVISIONAL
  on user instruction: written in this session to inspect d = 1 and d = 3
  behavior, not reviewed to the standard of examples 01-04, and subject to
  revision, renaming or removal. The file header and the README row say so.
  Decide its status before any release step; if it stays, it needs the same
  check the other examples had.
- Observed while writing it, recorded here rather than in user documentation
  (AGENTS.md): with one variable and `midpoint_spike = TRUE` the particles can
  collapse to a single tree (`unique_trees = 1`, no resampling), because there
  is no choice of splitting variable left. The ESS plot is then uninformative
  by construction. Not acted on.
- Note for the 2026-10-09 pending item (first figures users see should have two
  or more variables): example 05 adds one-variable histograms, but it is the
  last example and is not referenced from the README intro, so the first
  figures a reader meets are still the two-variable ones from example 01.
- Display precision, on user comment that the printed numbers are too long.
  `.pt_conditions()` had `format(digits = 7)` hard-coded and the regions table
  was a plain data frame, so it printed at `getOption("digits")`. Conditions
  are now written with the fewest digits that still reproduce the same
  division: `.pt_snap()` leaves slack between the largest value on the left
  and the smallest on the right, so a shortened threshold selects exactly the
  same rows. The digit count is chosen per threshold, not fixed; a fixed
  `signif(, 3)` already failed on 2 of 4 nodes in a d = 3, n = 900 example.
  Redundant clauses are also folded, since a region is an interval in each
  variable and only the tightest bound per side matters, and clauses now
  follow column order rather than the order of the divisions. Verified that
  all conditions reproduce `.pt_inside()` membership for d = 1, 2, 3, 20 and
  for n = 5000 (750 rows in total). `print.ptcompare_regions()` rounds
  probabilities to 4 and effects to 3 decimals for reading only; the returned
  data frame keeps full precision. The `threshold` column is shown at
  `signif(, 6)`, which still pads trailing zeros because a data frame column
  shares one decimal count; left as is.
- Access to the regions that differ most, on user question. Ranking was
  already easy (`regions()` is sorted, `summary()` shows the top 6,
  `plot(type = "regions")` the top 9), but there was no way to get the
  observations: only five exported functions, none returning membership, and
  `fit$data` / `fit$group` were not in the `@return` of `compare_pt`. Retyping
  the `condition` string was the only route. Added `region_rows()` (logical
  vector over the data in the order given) and `region_data()` (those
  observations with a group column, row names = positions), documented
  `data` / `group` in `compare_pt`, and added `rank_by` to `regions()` to
  match the argument `plot()` already had. `.pt_need_data()` said "This plot
  needs the data"; it is now worded for both callers.
- Added `vignettes/ptcompare.Rmd`, the package's introduction, after a review
  of what data could carry a convincing example. Written in English to match
  DESCRIPTION and the Rd files; the README stays Japanese. Three parts: a
  two-variable simulation where 15% of group 2 comes from a planted cluster
  (`region_rows()` then shows 14% of group 2 overall against 68% inside the
  reported region, so the recovery is checked, not asserted); `mclust::wdbc`,
  30 nucleus measurements, where the root split is on size but a depth-6
  region of large nuclei still separates on shape (`Nconcave_mean` below the
  threshold for 76% of benign against 9% of malignant); and a random split of
  the benign cases, which gives P = 0.07 and no region.
- Data selection for the vignette, recorded because the obvious choices are
  traps: `datasets::penguins` and `survival::flchain` are separated by a
  single variable, so the tree degenerates and the method's point is invisible;
  `kernlab::musk` is 80% ties; `ISLR::Caravan` and `kernlab::spam` are counts
  or zero-inflated. `caret::segmentationData` (57 channels, n = 2019, with its
  own train/test split as a built-in null) was the closest analogue of the
  paper's flow cytometry example but drags in ggplot2, recipes and pROC;
  `mclust::wdbc` was chosen instead because mclust imports only base packages.
  Conclusions on wdbc were unchanged for max_depth 6/8/10/12 and for 256/1024
  particles across three seeds.
- Vignette infrastructure: `VignetteBuilder: knitr`, and knitr, rmarkdown and
  mclust added to Suggests. The mclust chunks are guarded with
  `eval = requireNamespace("mclust")` so the vignette still builds without it.
  `R CMD build` produces a 192 KB tarball and `R CMD check --no-manual` is
  Status: OK, no ERRORs, WARNINGs or NOTEs.
- Vignette prose rewritten on user instruction: the first draft read as
  promotional copy (rhetorical openings, claims about what the method is for,
  emphasis on what competing tests cannot do). Rewritten in the register of an
  analysis report, stating what is computed and what the output means. Code
  and figures unchanged, so every reported number is identical. A scatter plot
  of the simulated two-variable data was added before the first fit, on user
  request, so the toy example is visible before it is analysed.
- Vignette reviewed by three independent agents given the rendered document
  only, the same prompt, and no background on the method. All three returned
  "not self-contained" on both questions asked, and converged on the same
  defects. Three were factual errors, each confirmed by computation here:
  (1) exp(effect) was read as a probability ratio ("seven times as likely");
  it is an odds ratio, the proportions being 60.8% against 17.7%. (2) "size
  has largely been used up as a way of telling the diagnoses apart" is false
  inside node 23: Area_mean medians are 674.5 against 1003.5, KS p = 9e-14,
  and the figure printed beneath the sentence shows it. (3) Nconcave_mean was
  described as an average count of concave sections, but its values lie in
  [0, 0.2012], like the other dimensionless shape indices.
- Other fixes from the same review: the stated row ordering of `regions()` was
  wrong (it is by prob_difference); `plot(type = "regions", top = 4)` dropped
  the second-ranked region because its own `min_n` default is 50 against the
  40 passed to `regions()`, so the plot call now passes `min_n = 40` and the
  panels match the table; the method section defined the effect as a gap
  between probabilities and the rest of the document as a log odds ratio; the
  group-to-level mapping (group 1 = B) was never stated; the single-node plot
  of region 23 duplicated a panel already shown; the 30 variable names were
  promised and never printed; the rank transform was never reconciled with
  "reported in the units of the original data"; `max_depth` was used for the
  display limit and the model limit without distinction; `effect_abs` and
  `prob_positive` were printed but never defined; the unsupported claim that a
  per-variable test "would report the same thing" was removed; the max_depth
  refit now prints depth_limited_fraction (0.861 against 0.726) so the check
  can actually fail; the null fit is printed, which supplies its group sizes;
  and every figure has a caption instead of knitr's chunk placeholder.
- Two claims in the reviews were wrong and were not acted on. All three
  inferred heap numbering for `node_id` from the shallow part of the first
  tree and concluded node 23 could not be at depth 6; ids are assigned in
  creation order and `left == 2 * node_id` fails in general. One review called
  `effect_abs` differing from `|effect|` a rounding artefact; they are
  different posterior functionals. Both misreadings were caused by missing
  definitions, which is why the definitions were added. A third complaint, YAML
  appearing as body text, was an artefact of knitting to Markdown for review.
  `R CMD check --no-manual` remains Status: OK; the tarball is 183 KB.
