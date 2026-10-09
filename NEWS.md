# ptcompare 0.0.1 (development)

- New two-sample R API and a C++17 SMC engine with exact partial-tree messages.
- Compact shared observation blocks, linear-time weight normalization and
  systematic resampling, and histogram-based split counts.
- Representative-tree region summaries, signed and absolute log-odds effects,
  zero-atom credible intervals, plots and particle diagnostics.
- Common pooled-rank coordinates with original-scale explanations; the default
  without supplied support is rank. Supplied support defaults to linear scaling.
- Exact-enumeration tests, statistical evaluation and RStudio examples.
- `pt_prior(midpoint_spike = TRUE)`: the spike-and-slab location prior of
  Awaya and Ma (2024, Appendix A) is now implemented and used by default; once
  a node is divided at its midpoint, all its descendants are too. Results
  saved in `evaluation/` were obtained before this change, i.e. with the
  equivalent of `midpoint_spike = FALSE`. Exact-enumeration tests cover both.
- `pt_control(resample_power = 0.5)`: resampling with probabilities
  proportional to W^kappa and weights W/a (Awaya and Ma, 2024, Section 3.1),
  with the paper's default kappa = 0.5; `resample_power = 1` gives standard
  resampling. Saved evaluation results predate this change (kappa = 1).
- Reported thresholds now reproduce the engine's region membership exactly
  (previously rounding could move an observation lying on a cut).
- Effect summaries (`effect_abs`, intervals, `prob_positive`) are computed by
  `regions()`/`summary()` for the returned nodes only; `fit$tree` keeps the
  signed `effect`. Split columns of leaves are NA.
- A warning is issued when no split is possible. Vector inputs for one
  variable, `print(fit, ...)` and graphical parameters in `plot()` are
  accepted; the region-plot histogram is aligned with its bins.
- New plots: `plot(fit, type = "regions")` shows the regions with the largest
  effects on scatter plots, as in Figure 14 of the supplement of Awaya and Ma
  (2024); `plot(fit, type = "levels")` shades the regions of each tree level by
  their probability of difference for two variables, as in Figure 7 of
  Soriano and Ma (2017). `type = "region"` now draws the region on a scatter
  plot of two variables.
- This is a finite HMPT variant; it does not reproduce legacy HMPT output.
