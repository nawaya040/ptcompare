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
- This is a finite HMPT variant; it does not reproduce legacy HMPT output.
