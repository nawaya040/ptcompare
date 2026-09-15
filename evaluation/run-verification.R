source("evaluation/common.R")
source("tests/testthat/helper-exact.R")
dir.create("evaluation/results",showWarnings=FALSE,recursive=TRUE)
x<-rbind(c(.1,.1),c(.2,.8),c(.3,.4),c(.4,.9),c(.8,.3),c(.9,.6),
          c(.1,.7),c(.2,.6),c(.4,.2),c(.6,.5),c(.7,.8),c(.9,.1))
g<-rep(0:1,each=6); p<-pt_prior(max_depth=2,split_grid=3,min_node=2,beta=1.5)
ref<-ref_enumerate(x,g,p); keys<-vapply(ref$tables,tree_key,character(1))
rows<-list()
for(M in c(64L,256L,1024L,4096L)) for(seed in 1:6) for(ess in c(0,1)) {
  set.seed(seed+800)
  f<-compare_pt(x,group=g,support=cbind(c(0,0),c(1,1)),prior=p,
    control=pt_control(particles=M,ess_threshold=ess,keep_particles=TRUE,effect_draws=256))
  pk<-vapply(f$particles,tree_key,character(1)); emp<-numeric(length(keys))
  agg<-tapply(f$diagnostics$weights,match(pk,keys),sum); emp[as.integer(names(agg))]<-agg
  rows[[length(rows)+1]]<-data.frame(particles=M,seed=seed,ess_threshold=ess,
    probability_error=f$prob_null-ref$prob_null,log_evidence_error=f$diagnostics$log_evidence-ref$log_evidence,
    total_variation=sum(abs(emp-ref$weights))/2)
}
write.csv(do.call(rbind,rows),"evaluation/results/exact-convergence.csv",row.names=FALSE)

set.seed(881); dat<-simulate_case("location",n1=400,n2=400,d=10)
rows<-list()
for(M in c(64L,256L,1024L)) for(seed in 1:6) {
  set.seed(seed+900)
  f<-compare_pt(dat$x,dat$y,transform="range",prior=pt_prior(max_depth=6,split_grid=7),
    control=pt_control(particles=M,effect_draws=512))
  r<-regions(f,prob=0); dg<-diagnostics(f)
  rows[[length(rows)+1]]<-data.frame(particles=M,seed=seed,prob_difference=f$prob_difference,
    top_variable=r$variable[1],top_threshold=r$threshold[1],top_condition=r$condition[1],
    effect=r$effect[1],seconds=dg$elapsed_seconds,root_ancestors=tail(dg$trace$root_ancestors,1))
}
write.csv(do.call(rbind,rows),"evaluation/results/stability.csv",row.names=FALSE)

# Separate finite-model sensitivity from Monte Carlo variation.
rows<-list()
for(depth in c(4L,6L,8L)) for(grid in c(3L,7L,15L)) for(seed in 1:3) {
  set.seed(seed+1000)
  f<-compare_pt(dat$x,dat$y,transform="range",prior=pt_prior(max_depth=depth,split_grid=grid),
    control=pt_control(particles=256,effect_draws=256))
  rows[[length(rows)+1]]<-data.frame(depth=depth,grid=grid,seed=seed,
    prob_difference=f$prob_difference,seconds=diagnostics(f)$elapsed_seconds,
    depth_fraction=diagnostics(f)$depth_limited_fraction)
}
write.csv(do.call(rbind,rows),"evaluation/results/model-sensitivity.csv",row.names=FALSE)
cat("Exact convergence, stability and model sensitivity completed.\n")
