source("evaluation/common.R")
f<-readRDS("evaluation/results/example-local-fit.rds")
t<-f$tree
mass<-function(bounds,g) {
  lo<-bounds$lo; hi<-bounds$hi
  if(g==1) return(prod(hi-lo))
  below<-max(0,min(hi[1],.35)-lo[1])
  above<-max(0,hi[1]-max(lo[1],.35))
  below*(hi[2]^5-lo[2]^5)+above*(hi[2]-lo[2])
}
r<-regions(f,prob=.9,min_n=40)
truth<-vapply(r$node_id,function(j) {
  parent<-node_bounds(f,j); child<-node_bounds(f,t$left[j])
  qlogis(mass(child,1)/mass(parent,1))-qlogis(mass(child,2)/mass(parent,2))
},numeric(1))
r$true_effect<-truth
r$correct_direction<-sign(r$effect)==sign(truth)
r$interval_contains_truth<-r$lower<=truth & truth<=r$upper
write.csv(r,"evaluation/results/example-true-effects.csv",row.names=FALSE)
print(r[,c("node_id","effect","true_effect","lower","upper","correct_direction","interval_contains_truth")])
