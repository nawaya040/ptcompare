# Statistical and computational specification (development 0.0.1)

## Scope and conditioning

Two independent iid samples of continuous vectors are compared. The finite model
uses common uniform distributions within terminal rectangles. Root support is
supplied or derived from the pooled range with 1% padding on each side. Pooled
constant columns are excluded. Support selection and the sample-size stopping
rule are conditioned on the pooled observations; this is an empirical specification,
not a fully generative model for those preprocessing decisions. No group-wise
standardization or screening is used. Repeated subjects and confounding need a
different model.

The transform='rank' option rescales pooled midranks to [0,1] using their
minimum and maximum. Ties remain tied. A stored monotone piecewise-linear
inverse returns cutpoints to the original measurement scale; observed region
membership is preserved. This alternative changes the geometric tree prior and
the original-scale terminal distributions. Posterior probabilities condition on
the estimated transform, and are empirically assessed under simulated nulls.
It is not an assertion of a full Bayesian posterior for marginal transformation
uncertainty. Following empirical evaluation, the default is 'rank' when support
is omitted, and 'range' when support is explicitly supplied. The choice is
stored and printed; transform='range' can always be requested explicitly.

## Tree and state model

The root is depth zero. Nodes at max_depth or with fewer than min_node pooled
observations are terminal. Other nodes split along a uniformly selected variable.
For numerical safety, a node is also terminal if a side of its unit-support
rectangle is at most 4*(L+1)*machine_epsilon; such stops are explicitly reported.
Relative fractions are k/(L+1), k=1,...,L, with probabilities proportional to
exp(-eta*n_A*abs(k/(L+1)-1/2)). Ties at the actual floating-point cut go left.
The first release uses the regular discrete prior without the legacy absorbing
midpoint restriction. It is a finite variant in the Awaya--Ma framework; legacy
numerical outputs are not its specification.

The internal state is D (different), S (shared), or A (absorbing shared). At depth k,
transition rows for parent D and S respectively are

    ((1-rho)*gamma*beta^(-k), (1-rho)*(1-gamma*beta^(-k)), rho)
    ((1-rho)*gamma*2^(-k),    (1-rho)*(1-gamma*2^(-k)),    rho).

The A row is (0,0,1). The root prior is the D row at k=0. rho is `prune`.
In D, theta_1 and theta_2 have independent Beta(c*r,c*(1-r)) priors. In S and A,
they share a single beta variable. c is `precision`; r is the split fraction.
Geometric partitioning continues through an absorbing sharing state: that state
enforces equal conditional distributions but does not require a uniform density
at that node. Common terminal uniforms occur at the finite stopping boundary.

## Likelihood and exact messages

Each internal emission is its beta integral multiplied by
r^(-n_left)*(1-r)^(-n_right). This expresses the likelihood relative to a common
uniform on root support. Terminal emissions equal one. Thus no binomial
coefficients occur: the data are observed labeled points, not an unordered count
observation model.

For own-state likelihood U_A(v), child-to-parent message M_B(v), and transition T,

    U_A(v) = emission_A(v) * M_left(v) * M_right(v)
    M_A(u) = sum_v T_depth(A)(u,v) U_A(v).

The root likelihood is sum_v prior_root(v) U_root(v). All calculations use log
space; structural zeros remain zero. Downward outside messages include both the
parent's emission and the sibling's upward message. They give exact own-state
posteriors conditional on the whole currently built tree.

Global H0 requires state != D at every INTERNAL node. Uniform terminal regions
introduce no further possible difference. Its numerator is computed by the same
recursion with each internal D emission set to zero, divided by tree likelihood.
The final H0 probability averages this ratio using normalized SMC weights.

## SMC target and proposal

Each particle extends its next eligible leaf in breadth-first order. A completed
tree becomes absorbing for subsequent iterations. Partial targets are tree prior
times exact partial-tree likelihood, with unexpanded leaves uniform. Given leaf
A and current partial tree, let pi_A(v) be its state probabilities using all
already observed branches. For split j,

    h(j) = sum_v pi_A(v) emission_j(v)
    q(j) = prior(j) h(j) / sum_k prior(k) h(k).

The incremental importance weight is the denominator. This includes sibling
evidence and telescopes to the declared complete-tree likelihood. Ancestors are
updated along the changed path. Systematic resampling from normalized weights
sets child weights to 1/M; it is triggered by ESS. Final log_evidence is the sum
of log normalizing increments (a stochastic estimate, not exact evidence).

The representative tree maximizes log(tree prior)+log(tree likelihood) among
sampled final particles. Its score includes variable and location priors, and is
distinct from an SMC importance weight. It is not a guaranteed global MAP tree.

## Effects and explanation

At node A, delta=logit(theta_1)-logit(theta_2). Positive values favor the left child
for group 1 conditional on A. P(D|data,T) refers to this conditional split, not
P(Q1(A)!=Q2(A)) or all differences anywhere inside A.

The mean signed effect uses digamma differences times P(D). Absolute mean and
95% equal-tailed intervals use deterministic low-discrepancy integration of
independent log-gamma variables. For shape a<1, log Gamma(a) is generated as
log Gamma(a+1)+log(U)/a, avoiding beta values rounded to endpoints. Group-swapped
quadrature is symmetrized; these summaries consume no R RNG. The zero atom of
weight 1-P(D) is included in the quantiles. Integration is approximate; the number
of points is reported in control$effect_draws and can be increased independently.

All region effects and intervals condition on the selected tree. Filtering nodes
by probability is exploratory and does not provide a frequentist FDR guarantee.
Nested regions overlap. Posterior probabilities are not p-values. Independent
SMC reruns assess algorithmic stability; they do not integrate tree selection
uncertainty into conditional local intervals.

## Computation and verification

Histogram plus cumulative sums replace suffix updates for candidate counts.
Stable linear-time normalization and systematic resampling replace quadratic
loops. Observation blocks are shared across particles and split only when their
node changes. Already-identical particle trees share candidate calculations.
Dense boundary arrays and all-particle R output are avoided by default.

tests/testthat/helper-exact.R independently enumerates internal states and all
small adaptive trees. Tests compare evidence, local states, H0 probabilities,
and full tree weights, with and without resampling. Additional tests cover group
exchange, row reordering, ties, endpoints, inverse transforms and public APIs.

Diagnostics include ESS BEFORE resampling, final weights, root ancestor count,
distinct tree count, finite-depth stopping and runtime. Particle payload counts
retained C++ nodes and observation blocks at iteration boundaries; it excludes
temporary objects and R memory and is not a process memory peak. High ESS does
not guarantee diverse ancestry, convergence, or stable interpretation.

Reference: Awaya and Ma (2024), JASA 119(545), 189--201,
doi:10.1080/01621459.2022.2105223. The MRS state model builds on Soriano and Ma
(2017), doi:10.1111/rssb.12180.
