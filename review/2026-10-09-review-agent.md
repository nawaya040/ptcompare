# ptcompare independent review — consolidated report (for AI agents)

- Date: 2026-10-09. Target: `main` @ d47daa4 (R package `ptcompare` 0.0.1).
- Method: three independent reviewers (R1, R2, R3) received identical instructions. Areas: (A) bugs, (B) fidelity to Awaya and Ma (2024, arXiv:2011.03121) apart from the agreed deviations, (C) usability and CRAN readiness. All three were read-only on the repository and ran small experiments in scratch directories.
- The coordinator independently re-verified the items marked **[V]** on R 4.5.2 / Windows.
- Scope rules (see AGENTS.md):
  - Do not pursue or document statistical performance limitations.
  - Keep changes minimal and faithful to the paper.
  - Agreed deviations are not defects: finite model defaults (depth 8 / grid 15 / min_node 20), rank as the default transform with range available, Appendix A midpoint spike on by default, the legacy-code fixes, beta^(-k) on the D row, the representative-tree definition, deterministic effect quadrature, and the H0 definition.

## Consensus verdict

| Area | R1 | R2 | R3 | Consolidated |
|---|---|---|---|---|
| A. Bugs | No correctness bugs; 2 performance defects, 1 validation bug | 1 reporting bug, 1 plot bug | 2 minor bugs, a few nits | **Core statistics correct.** No defect changes posterior quantities. Defects affect reported thresholds and plots, speed, and edge-case validation. |
| B. Fidelity | Faithful | Faithful except κ-resampling | Faithful except κ-resampling | **Faithful except one undocumented item**: the paper's κ-tempered resampling (B-1). |
| C. Usability | Works, check OK, doc gaps | Check OK, doc gaps | Check OK (1 NOTE: new submission), doc gaps | **Near CRAN-ready**; Rd quality is the main gap. |

`R CMD check --as-cran --no-manual`, with remote incoming checks off, was run by all three. R1 and R2 got Status OK. R3 got 1 NOTE, which was only "New submission". The test suite has 100 expectations, all passing, and takes 40–46 s.

## Findings

Fields: ID | severity | reported by | status.
- Severity scale: major / minor / nit. "Major" means user-visible wrong output or large cost; nothing here is critical.
- Each finding gives the location, the problem with evidence, and a minimal fix.

### A — Bugs and performance

**A-1 | major (reporting) | R2, R3 | confirmed [V]: range-transform thresholds can contradict the engine's split membership.**
- **Location:** `R/api.R:120-121` computes `threshold = support[j,1] + cut_unit*(support[j,2]-support[j,1])`. This value is used by `.pt_conditions` (`R/api.R:236-247`), by the region plot filter (`R/methods.R:162`), and by users who subset data with it. The engine assigns rows on the unit scale: `src/engine.cpp:146-150` and `split()` compare `(x-lo)/w <= cut`.
- **Problem:** the two roundings differ. An observation lying exactly on a cut can fall on the other side in original units.
- **Repro [V]:**
  ```r
  x<-matrix(rep(c(5,6,7),each=10)); y<-matrix(rep(c(8,9),each=15)); set.seed(1)
  f<-compare_pt(x,y,support=matrix(c(-0.2,9.4),1),
    prior=pt_prior(max_depth=1,split_grid=3,midpoint_spike=FALSE,eta=0),
    control=pt_control(particles=64,effect_draws=256))
  ```
  This gives `threshold = 6.9999999999999991`. The engine's left child has n1 = 30, but `sum(x <= threshold)` is 20.
- **Frequency:** R3 found 8957 mismatching (support, k) combinations in a search. R2 found 71 node mismatches in 60 fits on discrete data with a user-supplied support. Rank and padded-range fits: 0 mismatches (R1 179/569 nodes, R2, R3 18,000 splits).
- **Minimal fix:** make the reported threshold reproduce the engine's rule. After computing `thr`, snap it into `[max{x_orig : unit(x) <= cut}, min{x_orig : unit(x) > cut})`, for example by taking the largest original value going left when it is above `thr`. Alternatively, store the unit coordinates and do the plot's membership on the unit scale. Add a regression test with the repro.

**A-2 | minor | R1, R2, R3 | confirmed (R2 and R3 rendered PNGs): region-plot histogram is shifted half a bin and drops the last bin.**
- **Location:** `R/methods.R:175`, `graphics::matplot(mids, counts, type = "s", ...)`.
- **Problem:** with `type = "s"`, each level runs from `mids[i]` to `mids[i+1]`, so the bars sit half a bin to the right of the dashed threshold line and the last bin has no horizontal segment.
- **Fix:** `matplot(breaks, rbind(counts, counts[nrow(counts), ]), type = "s", ...)`. This also needs the `x` range of `breaks`.

**A-3 | major (performance) | R1 | confirmed [V]: effect quadrature dominates runtime.**
- **Location:** `R/api.R:126` → `.pt_effects` (`R/api.R:253-285`). `loggamma` calls `stats::qgamma` on `effect_draws` (default 4096) points, about 12 vectors per internal node of the representative tree.
- **Evidence:**
  - [V] default settings, n = 2000/group, d = 10: `.pt_effects` is 81% of 10.0 s; `pt_smc_cpp` is 1.9 s.
  - R1, defaults, n = 1000/group, d = 10: 85%.
  - R1, deep tree with 1255 internal nodes: 97% (50 s of 51.6 s).
- **Minimal fix:** compute effect summaries only for the nodes that are returned or displayed, or cache the uniform-to-gamma quantiles per node shape. Moving the `qgamma` work to C++ (R's `Rf_qgamma`) is also acceptable. Keep the result deterministic and RNG-free. This is a speed change, not a model change.

**A-4 | minor→major, condition-dependent (performance) | R1; R2 and R3 measured negligible | disputed / plausible: per-round payload accounting.**
- **Location:** `src/engine.cpp:298`, `peak=std::max(peak,payload(ps));`. `payload()` at `src/engine.cpp:220` hashes every node's obs-block pointer of every distinct tree, every round. It exists only for the `peak_particle_payload_bytes` diagnostic.
- **Evidence:**
  - R1, deep trees (min_node = 5, max_depth = 12, n = 2000/group): 37.1 s → 15.7 s without it. R1, n = 5000/group: 308 s → 194 s.
  - R2 and R3, default-depth settings: 16.9 s vs 16.7 s and 14.3 s vs 14.1 s.
- **Fix:** sample it every k rounds, maintain it incrementally, or compute it only on request. The fit is unchanged (R1 confirmed identical P(diff)).

**A-5 | minor | R1, R3 | confirmed: the rank transform runs an irrelevant range-width check.**
- **Location:** `R/api.R:211` sets `support <- ranges` for rank. The width check at `R/api.R:224-226` then rejects finite extreme-scale data with "Data scale is too extreme for safe range scaling".
- **Repro:** `compare_pt(matrix(c(-1e308,1e308,0,1,2)), matrix(3:7))`
- **Fix:** skip the width check and the scaling `sweep` when `transform == "rank"`.

**A-6 | minor (silent degenerate output) | R1, R2 | confirmed [V]: no split possible gives P(diff) = 0 silently.**
- **Location:** `src/engine.cpp:116` (stop when pooled n < `min_node`) and `R/api.R:128-153`.
- **Repro:** `compare_pt(matrix(rnorm(8)), matrix(rnorm(8,10)))$prob_difference` returns `0`.
- **Problem:** the value is a property of the finite model, not evidence for equality.
- **Fix:** `warning()` when the root is a leaf ("no split possible: pooled n < min_node; P(H0) = 1 by construction"), and note it in `print`/`summary`.

### B — Fidelity to the paper

**B-1 | minor (algorithmic default; target unchanged) | R2, R3 (R1 did not flag) | confirmed [V]: κ-tempered resampling not implemented and not documented.**
- **Paper:** §3.1, the paragraph before the algorithm box, resamples with `a_t^m ∝ (W_t^m)^κ`, κ ∈ (0,1], and sets new weights `∝ W/a` (following Lu et al. 2013). It recommends κ = 0.5 as "our default choice in all of our later examples". The algorithm box suggests an ESS trigger such as M/10. Legacy `research/HMPT` uses mixing 0.5 and threshold 0.1.
- **Code:** `src/engine.cpp:283-293` uses systematic resampling ∝ W with weights reset to 1/M (κ = 1). The default is `pt_control(ess_threshold = 0.5)` (`R/api.R:9`).
- **Status:** κ = 1 is within the paper's family and the target is unchanged, but this is not among the agreed deviations, and neither the spec nor `?pt_control` mentions it.
- **Minimal faithful fix:**
  - Add `pt_control(resample_power = 0.5)`. Draw ancestors systematically from `a ∝ W^κ` and set the new log-weights to `log W − log a` (then normalize).
  - Decide whether the default `ess_threshold` should follow the paper (0.1).
  - Extend `test-exact.R` with κ = 0.5 under frequent resampling.
  - Document it in `inst/STATISTICAL-SPECIFICATION.md`.
- **Decision needed from the user:** implement (paper default κ = 0.5) or document κ = 1 as a deviation.

Items checked by all three and found faithful:
- Emissions, eq. (4) with the μ-ratio of eq. (6).
- The F.2.1 transition matrix with k = node depth, root prior = D row at k = 0, A absorbing.
- The eq. (2) location prior with pooled n(A), uniform variable choice, −log d.
- The Prop. 3.1 / Cor. 3.1 proposal and incremental weight, with exact outside messages.
- The Appendix A lock: R1 found 0 violations in 569 nodes, and R1's independent enumeration (d = 1, ties on cuts, 7 vs 4, depth 3) agrees with SMC both with and without the spike.
- The H0 recursion, PMAP, the representative tree, and E|logOR|.

### C — Usability and CRAN readiness

| ID | Sev. | By | Location | Issue | Minimal fix |
|---|---|---|---|---|---|
| C-1 | minor | R1 R2 R3 | DESCRIPTION; `man/compare_pt.Rd:29,31,41` | Roxygen markdown not enabled, so the help shows literal `[pt_prior()]` and backticks | Add `Roxygen: list(markdown = TRUE)`; re-roxygenise |
| C-2 | minor | R1 R2 R3 | `man/*.Rd` | Only `compare_pt` has `\examples`; `summary.ptcompare_fit` (`top`) and `print` are undocumented; `\value` of `compare_pt` and `regions` does not itemize components and columns (`prob_null`, `tree`, `cut_unit`, `fraction`, `effect`, `effect_abs`, `lower`, `upper`, `prob_positive`, `a1..b2`, `stop_reason`, `particles`) | Add short examples, a summary Rd, and itemized `\value` |
| C-3 | nit | R1 R2 R3 | `R/api.R:175` | A plain numeric vector (d = 1) is rejected; d = 1 support must be a matrix | Coerce vector → 1-column matrix; length-2 support → `matrix(support, 1)` |
| C-4 | nit | R1 R2 R3 | `inst/STATISTICAL-SPECIFICATION.md` "Scope and conditioning" | "Root support is supplied or derived from the pooled range with 1% padding" is false for the rank default | Reword: padding applies to range only |
| C-5 | minor | R3 | `.pt_dots` `R/methods.R:85`, used by print/plot | `print(fit, digits = 3)` and `plot(fit, main = "x")` error | Ignore `...` in print methods; pass graphics parameters in `plot` |
| C-6 | nit | R1 R3 | `R/methods.R:139` (3 digits), `R/api.R` conditions (6 digits), summary (4 digits) | Rounded thresholds can display a boundary value on the wrong side, e.g. 0.5625 → "0.562" | More digits, or print exact thresholds in `regions()` |
| C-7 | nit | R3 | `src/engine.cpp:46-48` and `table()` | Leaf rows show `cut_unit = 0`, `fraction = 0.5`, `variable_id = 0` | NA for leaves |
| C-8 | nit | R3 | group handling | Numeric labels are ordered by first appearance (`c(1,1,0,0)` → group 1 = "1"), which flips the sign of `effect`; this is documented | Optionally sort numeric/logical labels; at least emphasize in Rd |
| C-9 | nit | R2 | `R/api.R:130-150` | `fit$particles` (`keep_particles = TRUE`) are raw engine tables without names or thresholds; undocumented | Document, or post-process like `fit$tree` |
| C-10 | nit | R2 | `pt_prior`, `plot` Rd | Enforced upper limits (`max_depth` ≤ 20, `split_grid` ≤ 127, `eta`/`beta` ≤ 100, plot depth ≤ 6) are undocumented | State in `@param` |
| C-11 | nit | R1 | `pt_prior()`, `pt_control()` | No print methods; printing shows a raw list with `attr(,"class")` | Add small print methods |
| C-12 | nit | R1 R2 R3 | `R/methods.R` tree plot `sub=`; `R/api.R:117`; `cat(..., "\n")` at `R/methods.R:40,65` | Legend duplicated in clipped `sub=`; dead `j <-`; trailing space | Remove |
| C-13 | nit | R3 | tests | Suite takes 40–46 s (12,000-particle enumeration run 4×) | `skip_on_cran()` for the heavy SMC-vs-enumeration test, or fewer particles |

## Not verified by any reviewer
- Linux, macOS, R-devel, sanitizers and valgrind. Memory safety was assessed by code reading only; all three judged the copy-on-write sharing safe.
- Remote CRAN incoming checks (spelling, URLs).
- Supplement B text and a line-by-line Appendix A.2 derivation. Correctness there rests on code reading and the exact-enumeration tests.
- Example 03 (R1, R3) and the `evaluation/` scripts (by instruction).

## Suggested order of work (minimal changes)
1. A-1. This is the only defect that makes user-visible output (region conditions) wrong. Fix it with a regression test.
2. B-1. This needs a user decision first.
3. A-2, A-5, A-6. Small and local.
4. A-3, plus A-4 if cheap. Speed only; verify identical outputs.
5. C-1 to C-5 before CRAN; the rest are optional.

Reviewer scratch artifacts, outside the repository: `scratchpad/review-r7q2` (R1), `review-k7q2` (R2), `review-q7k2` (R3).

## Resolution (2026-10-09, same day)

The user chose to implement κ = 0.5 as in the paper and to fix the other findings.

| ID | Outcome |
|---|---|
| B-1 | Fixed. `pt_control(resample_power = 0.5)`: ancestors are drawn with probability ∝ W^κ and weighted W/(M·a); the sum is added to log-evidence and renormalized (`src/engine.cpp`, resampling block). `ess_threshold` default unchanged (0.5). Exact-enumeration tests cover no resampling, κ = 0.5 and κ = 1. |
| A-1 | Fixed. `.pt_snap()` in `R/api.R` moves a threshold only when rounding would cross an observation. Regression tests cover the repro and 20 random discrete-data fits. |
| A-2 | Fixed. The step histogram is now drawn over bin edges (checked visually: the jump sits at the true gate 0.35). |
| A-3 | Fixed. Effect quadrature runs lazily in `regions()`/`summary()` for returned nodes only. Default case (n = 2000/group, d = 10): 10.0 s → 2.6 s, plus `summary()` 0.3 s. `fit$tree` keeps the signed `effect`; `effect_abs`/`lower`/`upper`/`prob_positive` moved to `regions()`. |
| A-4 | Fixed. Payload is measured at steps 1, 2, 4, … and at the end. The deep-tree case now takes 8.5 s. |
| A-5 | Fixed. Width check and scaling only for range. |
| A-6 | Fixed. Warning when no split is possible, plus a note in `print()`. |
| C-1, C-2 | Fixed. Roxygen markdown is enabled; examples were added to all exported topics; a new `summary.ptcompare_fit` Rd was added; `\value` of `compare_pt` and `regions` is itemized. |
| C-3, C-4, C-5, C-7, C-9, C-10, C-11, C-12 | Fixed. |
| C-6 | Partly fixed: more digits (conditions 7, tree labels 4). Exact values are in `regions()$threshold`. |
| C-8 | Documentation only; the behavior is unchanged. |
| C-13 | Fixed. `skip_on_cran()` on the SMC-vs-enumeration test. |

Verification:
- Tests: 362 expectations, 0 failures, with `NOT_CRAN=true`.
- `R CMD check --as-cran --no-manual`: 1 NOTE, "unable to verify current time" (network only).
- All four examples run (9 s, 3.6 s, 52 s, 6.2 s).
