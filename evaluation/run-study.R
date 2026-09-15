# Run from the project root. Results are checkpointed after each data set.
source("evaluation/common.R")
dir.create("evaluation/results",recursive=TRUE,showWarnings=FALSE)
configs<-rbind(
  expand.grid(d=c(2L,10L),scenario=c("null","location","dispersion","local","correlation")),
  expand.grid(d=c(50L,100L),scenario=c("null","local","correlation")),
  data.frame(d=10L,scenario=c("heavy_null","heavy_location","unbalanced_null","unbalanced_location")))
configs$reps<-ifelse(configs$d<=10,ifelse(grepl("null",configs$scenario),200L,50L),
                     ifelse(configs$scenario=="null",40L,25L))
configs$reps[grepl("heavy|unbalanced",configs$scenario)]<-ifelse(grepl("null",configs$scenario[grepl("heavy|unbalanced",configs$scenario)]),80L,40L)
jobs<-do.call(rbind,lapply(seq_len(nrow(configs)),function(i)
  data.frame(config=i,d=configs$d[i],scenario=as.character(configs$scenario[i]),rep=seq_len(configs$reps[i]))))
path<-"evaluation/results/study.csv"
done<-if(file.exists(path)) read.csv(path) else data.frame()
names(done)<-sub("[.]elapsed$","",names(done))
completed<-if(nrow(done)) done$job else integer()
for(job in setdiff(seq_len(nrow(jobs)),completed)) {
  a<-jobs[job,]; n1<-if(grepl("unbalanced",a$scenario)) 120L else 300L
  n2<-if(grepl("unbalanced",a$scenario)) 600L else 300L
  data_seed<-100000L+job; fit_seed<-200000L+job
  set.seed(data_seed); dat<-simulate_case(a$scenario,n1,n2,a$d)
  set.seed(fit_seed)
  fit<-compare_pt(dat$x,dat$y,transform="range",prior=pt_prior(max_depth=6,split_grid=7,min_node=20),
                   control=pt_control(particles=64,keep_data=FALSE,effect_draws=256))
  set.seed(300000L+job); base<-baseline_scores(dat$x,dat$y,mrs=a$d==2)
  dg<-diagnostics(fit); r<-regions(fit,prob=0)
  loc<-if(a$scenario=="local") localization(fit) else c(top_iou=NA,best5_iou=NA)
  row<-data.frame(job=job,a,n1=n1,n2=n2,data_seed=data_seed,fit_seed=fit_seed,
    pt_score=-log(max(fit$prob_null,.Machine$double.xmin)),prob_difference=fit$prob_difference,
    seconds=dg$elapsed_seconds,final_ess=dg$final_ess,root_ancestors=tail(dg$trace$root_ancestors,1),
    payload_mb=dg$peak_particle_payload_bytes/1024^2,depth_fraction=dg$depth_limited_fraction,
    internal_nodes=sum(fit$tree$left>0),top_variable=if(nrow(r))r$variable[1] else NA,
    top_iou=loc[1],best5_iou=loc[2],as.list(base),check.names=FALSE)
  done<-rbind(done,row); write.csv(done,path,row.names=FALSE)
  if(job%%10L==0L) cat(sprintf("%s %d/%d d=%d %s: %.2fs, P(diff)=%.3f\n",
    format(Sys.time(),"%H:%M:%S"),job,nrow(jobs),a$d,a$scenario,dg$elapsed_seconds,fit$prob_difference))
}
writeLines(capture.output(sessionInfo()),"evaluation/results/session-info.txt")
cat("Study completed:",nrow(done),"data sets\n")
