source("scripts/load-local.R")
dir.create("evaluation/figures",showWarnings=FALSE)
sink("evaluation/results/examples.log",split=TRUE)
files<-sort(list.files("examples",pattern="[.]R$",full.names=TRUE))
times<-list()
for(i in seq_along(files)) {
  env<-new.env(parent=globalenv()); start<-proc.time()[[3]]
  if(i!=3) png(sprintf("evaluation/figures/example%d-%%02d.png",i),1400,900,res=140)
  cat("\nExecuting",files[i],"\n")
  source(files[i],local=env)
  if(i!=3) dev.off()
  times[[i]]<-data.frame(example=files[i],seconds=proc.time()[[3]]-start)
  if(i==1) saveRDS(env$fit,"evaluation/results/example-local-fit.rds")
  if(i==2) {
    saveRDS(env$dependence_fit,"evaluation/results/example-dependence-fit.rds")
    saveRDS(env$null_fit,"evaluation/results/example-null-fit.rds")
  }
  if(i==3) write.csv(env$stability_table,"evaluation/results/example-stability.csv",row.names=FALSE)
  if(i==4) saveRDS(env$rank_fit,"evaluation/results/example-rank-fit.rds")
}
write.csv(do.call(rbind,times),"evaluation/results/example-times.csv",row.names=FALSE)
sink()
