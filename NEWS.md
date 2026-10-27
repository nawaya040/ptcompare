# ptcompare 0.1.0

First release.

* `compare_pt()` compares two samples of continuous multivariate observations
  with the hidden Markov Pólya tree of Awaya and Ma (2024)
  <doi:10.1080/01621459.2022.2105223>. It returns the posterior probability
  that the two distributions differ, a representative partition tree, and local
  effect sizes for the divisions of that tree. Posterior computation combines
  sequential Monte Carlo over partitions with recursive message passing. The
  two samples are supplied either as two matrices or as one data set with a
  grouping variable.

* `pt_prior()` sets the model: the depth of the tree, the number of candidate
  split points, the state transition probabilities, and the spike-and-slab
  prior on split locations of Awaya and Ma (2024, Appendix A), which is used by
  default. `pt_control()` sets the computational budget: the number of
  particles, the resampling threshold and power, and the number of posterior
  draws used for the effect sizes.

* Variables are put on a common scale before the space is divided.
  `transform = "rank"` uses pooled midranks and is the default;
  `transform = "range"` uses the linear scaling of the paper and is the default
  when `support` is supplied. Thresholds and region conditions are reported in
  the original units.

* `regions()` lists the divisions of the representative tree with their
  probability of a difference and their signed and absolute log odds ratios,
  ordered by either quantity. `region_rows()` and `region_data()` return the
  observations of a region.

* `plot()` draws the representative tree, the regions with the largest
  differences on scatter plots, the regions of each tree depth shaded by their
  probability of a difference (two variables only), and the effective sample
  size of the sequential Monte Carlo. `summary()` and `diagnostics()` report
  the settings and the computation.

* `vignette("ptcompare")` works through a simulated comparison in two variables
  and a comparison of 30 measurements of cell nuclei.
