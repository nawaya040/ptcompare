source("evaluation/common.R")
a<-commandArgs(trailingOnly=TRUE)
n<-as.integer(a[1]); d<-as.integer(a[2]); M<-as.integer(a[3]); out<-a[4]
transform<-if(length(a)>=5) a[5] else "range"
set.seed(2200+n+d); dat<-simulate_case("location",n1=n%/%2,n2=n-n%/%2,d=d)
set.seed(2300)
f<-compare_pt(dat$x,dat$y,prior=pt_prior(max_depth=6,split_grid=7),
  control=pt_control(particles=M,keep_data=FALSE,effect_draws=256),transform=transform)
dg<-diagnostics(f)
write.csv(data.frame(total_n=n,d=d,particles=M,transform=transform,seconds=dg$elapsed_seconds,
  payload_mb=dg$peak_particle_payload_bytes/1024^2,prob_difference=f$prob_difference,
  nodes=nrow(f$tree),root_ancestors=tail(dg$trace$root_ancestors,1)),out,row.names=FALSE)
