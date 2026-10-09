test_that("group exchange, row order and interfaces preserve the target", {
  set.seed(11); x <- matrix(rnorm(100),ncol=2); y <- matrix(rnorm(120,.4),ncol=2)
  p <- pt_prior(max_depth=4,split_grid=3,min_node=12); c <- pt_control(particles=64,effect_draws=512)
  set.seed(8); a <- compare_pt(x,y,prior=p,control=c)
  set.seed(8); b <- compare_pt(y,x,prior=p,control=c)
  expect_equal(a$prob_null,b$prob_null,tolerance=1e-10)
  expect_equal(a$tree$prob_difference,b$tree$prob_difference,tolerance=1e-10)
  expect_equal(a$tree$effect,-b$tree$effect,tolerance=1e-10)
  ra <- regions(a,prob=0); rb <- regions(b,prob=0); rb <- rb[match(ra$node_id,rb$node_id),]
  expect_equal(ra$effect_abs,rb$effect_abs,tolerance=1e-10)
  expect_equal(ra$lower,-rb$upper,tolerance=1e-10)
  set.seed(8); z <- compare_pt(x[nrow(x):1,],y[nrow(y):1,],prior=p,control=c)
  expect_equal(a$prob_null,z$prob_null,tolerance=1e-12)
  set.seed(8); z <- compare_pt(rbind(x,y),group=factor(rep(c("A","B"),c(nrow(x),nrow(y)))),prior=p,control=c)
  expect_equal(a$tree,z$tree)
  state <- .Random.seed; regions(a); summary(a); diagnostics(a)
  expect_identical(.Random.seed,state)
})

test_that("ties, endpoints, constant columns and inverse scaling are handled", {
  x <- cbind(a=c(0,.5,.5,1),b=1); y <- cbind(a=c(0,.25,.75,1),b=1)
  p <- pt_prior(max_depth=3,split_grid=3,min_node=2)
  expect_warning(fit <- compare_pt(x,y,support=cbind(c(0,0),c(1,2)),prior=p,
                                   control=pt_control(particles=24,effect_draws=256)),"constant")
  t <- fit$tree
  for(j in which(t$left>0)) {
    expect_equal(t$n1[t$left[j]]+t$n1[t$right[j]],t$n1[j])
    expect_equal(t$n2[t$left[j]]+t$n2[t$right[j]],t$n2[j])
  }
  set.seed(30); a <- compare_pt(x[,1,drop=FALSE],y[,1,drop=FALSE],prior=p,control=pt_control(particles=32,effect_draws=256))
  set.seed(30); b <- compare_pt(100+20*x[,1,drop=FALSE],100+20*y[,1,drop=FALSE],prior=p,control=pt_control(particles=32,effect_draws=256))
  expect_equal(a$prob_null,b$prob_null,tolerance=1e-9)
  expect_equal(100+20*a$tree$threshold,b$tree$threshold,tolerance=1e-9)
})

test_that("inputs fail clearly", {
  x <- matrix(1:8,4,2)
  expect_error(compare_pt(x),"exactly one")
  expect_error(compare_pt(x,x,group=1:4),"exactly one")
  expect_error(compare_pt(x,group=c(1,1,1,2)),"two observations")
  expect_error(compare_pt(x,group=1:4),"two groups")
  expect_error(compare_pt(x,group=c(1,1,NA,2)),"nonmissing")
  expect_error(compare_pt(x*NA,x),"non-finite")
  expect_error(compare_pt(x+1i,x),"numeric")
  expect_error(compare_pt(x*0,x*0),"constant")
  expect_error(compare_pt(data.frame(a=letters[1:4]),x),"numeric")
  expect_error(compare_pt(x,x,support=matrix(0,2,2)),"support")
  expect_error(pt_control(particles=2.5),"integer")
  expect_error(pt_prior(split_grid=4),"odd")
  expect_error(pt_prior(prune=0),"strictly")
  expect_error(pt_prior(beta=.5),"between")
  expect_error(pt_control(resample_power=0),"greater than 0")
  expect_error(pt_control(resample_power=1.5),"between")
  c <- pt_control(); c$particles <- -1
  expect_error(compare_pt(x,x,control=c),"particles")
})

test_that("log-gamma effect quadrature matches known moments", {
  t <- data.frame(left=2,a1=2,b1=4,a2=4,b2=2,prob_difference=.8)
  ans <- ptcompare:::.pt_effects(t,16384)
  set.seed(5); d <- qlogis(rbeta(250000,2,4))-qlogis(rbeta(250000,4,2))
  expect_equal(ans$effect_abs,.8*mean(abs(d)),tolerance=.012)
  expect_lte(ans$lower,0); expect_gte(ans$upper,0)
  expect_true(is.finite(ans$effect_abs))
  t$prob_difference <- 0
  z <- ptcompare:::.pt_effects(t,256)
  expect_equal(z$effect_abs,0); expect_equal(z$lower,0); expect_equal(z$upper,0)
})

test_that("plots work with and without stored data", {
  set.seed(3); x <- matrix(rnorm(80),ncol=2); y <- x+.6
  a <- compare_pt(x,y,control=pt_control(particles=16,effect_draws=256),prior=pt_prior(max_depth=3,split_grid=3))
  b <- compare_pt(x[,1],y[,1],control=pt_control(particles=16,effect_draws=256),prior=pt_prior(max_depth=3,split_grid=3))
  path <- tempfile(fileext=".pdf"); grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)})
  state <- .Random.seed
  expect_silent(plot(a)); expect_silent(plot(a,type="regions",prob=0,min_n=0)); expect_silent(plot(a,type="diagnostics"))
  expect_silent(plot(a,main="custom")); expect_silent(plot(a,type="regions",node=1,main="custom"))
  expect_silent(plot(a,type="regions",node=c(1,2)))
  expect_message(plot(a,type="regions",prob=1,min_n=1000),"No region")
  expect_error(plot(a,type="regions",node=nrow(a$tree)),"node_id")
  expect_silent(plot(a,type="regions",top=4,prob=0,min_n=2)); expect_silent(plot(a,type="levels"))
  expect_silent(plot(b,type="regions",prob=0,min_n=0)); expect_silent(plot(b,type="regions",node=1))
  expect_error(plot(b,type="levels"),"two variables")
  expect_output(print(a,digits=3),"P\\(distributions differ\\)")
  expect_identical(.Random.seed,state)
  a$data <- NULL; expect_error(plot(a,type="regions"),"keep_data")
})

test_that("extreme nested boundary splits terminate safely", {
  set.seed(45)
  f <- compare_pt(matrix(rep(0,10),ncol=1),matrix(rep(1,10),ncol=1),
    support=matrix(c(0,1),1,2),prior=pt_prior(max_depth=20,split_grid=127,eta=0,min_node=2),
    control=pt_control(particles=4,effect_draws=256))
  expect_true(is.finite(f$prob_null))
  expect_true(all(is.finite(f$tree$effect[!is.na(f$tree$effect)])))
  expect_true(any(f$tree$stop_reason %in% c("max_depth","numerical precision")))
})

test_that("pooled ranks preserve labels, ties and original-scale membership", {
  x <- cbind(A=c(0,0,1,2,3,9),B=c(1,2,3,5,7,8))
  y <- cbind(A=c(0,1,1,2,4,10),B=c(1,3,4,6,8,9))
  set.seed(721)
  f <- compare_pt(x,y,transform="rank",prior=pt_prior(max_depth=3,split_grid=7,min_node=2),
                  control=pt_control(particles=32,effect_draws=256))
  data<-rbind(x,y); group<-rep(1:2,each=6); t<-f$tree
  membership<-list(seq_len(nrow(data)))
  for (j in which(t$left>0)) {
    rows<-membership[[j]]; values<-data[rows,t$variable_id[j]]
    membership[[t$left[j]]]<-rows[values<=t$threshold[j]]
    membership[[t$right[j]]]<-rows[values>t$threshold[j]]
  }
  for(j in seq_len(nrow(t))) {
    expect_equal(tabulate(group[membership[[j]]],2L),c(t$n1[j],t$n2[j]))
  }
  set.seed(721)
  h <- compare_pt(y,x,transform="rank",prior=pt_prior(max_depth=3,split_grid=7,min_node=2),
                  control=pt_control(particles=32,effect_draws=256))
  expect_equal(f$prob_null,h$prob_null,tolerance=1e-10)
  expect_error(compare_pt(x,y,transform="rank",support=cbind(c(-1,-1),c(11,11))),"support")
})

test_that("default coordinates are explicit and support-aware", {
  x<-matrix(c(.1,.2,.3,.4),ncol=1); y<-matrix(c(.3,.5,.7,.9),ncol=1)
  c<-pt_control(particles=8,effect_draws=256); p<-pt_prior(min_node=4)
  set.seed(4); a<-compare_pt(x,y,control=c,prior=p)
  set.seed(4); b<-compare_pt(x,y,control=c,prior=p,transform="rank")
  expect_identical(a$transform,"rank")
  expect_equal(a$tree,b$tree)
  z<-compare_pt(x,y,control=c,prior=p,support=matrix(c(0,1),1,2))
  expect_identical(z$transform,"range")
})

test_that("range thresholds reproduce the engine's region membership", {
  # A data value lies exactly on a cut; rounding of lo + cut*width must not
  # move it to the other side in original units.
  x <- matrix(rep(c(5,6,7),each=10)); y <- matrix(rep(c(8,9),each=15))
  set.seed(1)
  f <- compare_pt(x,y,support=matrix(c(-0.2,9.4),1),
                  prior=pt_prior(max_depth=1,split_grid=3,midpoint_spike=FALSE,eta=0),
                  control=pt_control(particles=64,effect_draws=256))
  t <- f$tree
  expect_equal(c(sum(x<=t$threshold[1]),sum(y<=t$threshold[1])),c(t$n1[2],t$n2[2]))
  set.seed(2)
  for (i in 1:20) {
    lo <- -round(runif(1),1); hi <- 9+round(runif(1),1)
    a <- matrix(sample(0:9,40,TRUE),ncol=2); b <- matrix(sample(0:9,40,TRUE),ncol=2)
    f <- compare_pt(a,b,support=cbind(c(lo,lo),c(hi,hi)),
                    prior=pt_prior(max_depth=3,split_grid=7,min_node=4),
                    control=pt_control(particles=16,effect_draws=256))
    t <- f$tree; data <- rbind(a,b); group <- rep(1:2,each=20)
    membership <- list(seq_len(nrow(data)))
    for (j in which(t$left>0)) {
      rows <- membership[[j]]; v <- data[rows,t$variable_id[j]]
      membership[[t$left[j]]] <- rows[v<=t$threshold[j]]
      membership[[t$right[j]]] <- rows[v>t$threshold[j]]
    }
    for (j in seq_len(nrow(t))) expect_equal(tabulate(group[membership[[j]]],2L),c(t$n1[j],t$n2[j]))
  }
})

test_that("vector inputs, extreme scales and leaf fields are handled", {
  set.seed(9); c <- pt_control(particles=8,effect_draws=256); p <- pt_prior(max_depth=2,min_node=4)
  set.seed(1); a <- compare_pt(c(.1,.2,.3,.4,.5),c(.3,.5,.7,.9,1),prior=p,control=c)
  set.seed(1); b <- compare_pt(matrix(c(.1,.2,.3,.4,.5)),matrix(c(.3,.5,.7,.9,1)),prior=p,control=c)
  expect_equal(a$tree,b$tree)
  set.seed(1); z <- compare_pt(c(.1,.2,.3,.4,.5),c(.3,.5,.7,.9,1),support=c(0,1),prior=p,control=c)
  expect_identical(z$transform,"range")
  expect_silent(compare_pt(matrix(c(-1e308,1e308,0,1,2)),matrix(3:7),prior=p,control=c))
  leaf <- a$tree$left==0
  expect_true(all(is.na(a$tree$variable_id[leaf]) & is.na(a$tree$cut_unit[leaf]) & is.na(a$tree$fraction[leaf])))
  expect_output(print(pt_prior()),"midpoint_spike")
  expect_output(print(pt_control()),"resample_power")
})
