test_that("fixed-tree messages agree with complete state enumeration", {
  t <- data.frame(node_id=1:7,parent_id=c(0,1,1,2,2,3,3),
    left=c(2,4,6,0,0,0,0),right=c(3,5,7,0,0,0,0),depth=c(0,1,1,2,2,2,2),
    fraction=c(.25,.5,.75,.5,.5,.5,.5),n1=c(17,6,11,2,4,8,3),n2=c(21,12,9,10,2,4,5))
  for (p in list(pt_prior(max_depth=2,beta=1.8),pt_prior(max_depth=2,gamma=.7,prune=.1,precision=2))) {
    ref <- ref_fixed(t,p)
    cpp <- ptcompare:::pt_fixed_cpp(as.matrix(t[,c("left","right")]),as.matrix(t[,c("n1","n2")]),
                                     t$fraction,as.integer(t$depth),unclass(p))
    expect_equal(cpp$log_evidence,ref$log_evidence,tolerance=1e-11)
    expect_equal(cpp$prob_null,ref$prob_null,tolerance=1e-11)
    expect_equal(cpp$tree$prob_difference,ref$prob_diff,tolerance=1e-11)
  }
})

test_that("SMC targets the fully enumerated adaptive tree posterior", {
  skip_on_cran()
  x <- rbind(c(.1,.1),c(.2,.8),c(.3,.4),c(.4,.9),c(.8,.3),c(.9,.6),
             c(.1,.7),c(.2,.6),c(.4,.2),c(.6,.5),c(.7,.8),c(.9,.1))
  g <- rep(0:1,each=6)
  for (spike in c(TRUE,FALSE)) {
  p <- pt_prior(max_depth=2,split_grid=3,min_node=2,beta=1.5,eta=.1,midpoint_spike=spike)
  ref <- ref_enumerate(x,g,p)
  expect_equal(sum(ref$weights),1)
  expect_length(ref$weights,if (spike) 152 else 216)
  keys <- vapply(ref$tables,tree_key,character(1))
  # No resampling; frequent resampling with the paper's kappa = 0.5 (default)
  # and, for the spike model, with standard resampling (kappa = 1).
  runs <- list(c(0,.5),c(1,.5)); if (spike) runs <- c(runs,list(c(1,1)))
  for (run in runs) {
    threshold <- run[1]
    set.seed(713)
    fit <- compare_pt(x,group=g,support=cbind(c(0,0),c(1,1)),prior=p,
                     control=pt_control(particles=12000,ess_threshold=threshold,keep_particles=TRUE,
                                        effect_draws=256,resample_power=run[2]))
    expect_equal(fit$prob_null,ref$prob_null,tolerance=.012)
    expect_equal(fit$diagnostics$log_evidence,ref$log_evidence,tolerance=.025)
    empirical <- numeric(length(keys))
    pk <- vapply(fit$particles,tree_key,character(1))
    expect_false(anyNA(match(pk,keys)))
    agg <- tapply(fit$diagnostics$weights,match(pk,keys),sum)
    empirical[as.integer(names(agg))] <- agg
    expect_lt(sum(abs(empirical-ref$weights))/2,.07)
    if (threshold == 1) expect_gt(fit$diagnostics$resamples,0)
  }
  }
})

test_that("unsplit model has unit null probability and likelihood ratio", {
  x <- matrix(c(.1,.2,.8,.9),ncol=1)
  expect_warning(fit <- compare_pt(x,group=c(0,0,1,1),prior=pt_prior(min_node=10),
                                   control=pt_control(particles=8,effect_draws=256)),"No split")
  expect_equal(fit$prob_null,1,tolerance=1e-12)
  expect_equal(fit$diagnostics$log_evidence,0)
  expect_equal(nrow(regions(fit,prob=0)),0)
})
