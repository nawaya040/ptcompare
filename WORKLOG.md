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
