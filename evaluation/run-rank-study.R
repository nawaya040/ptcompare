source("evaluation/common.R")
jobs<-read.csv("evaluation/results/study.csv")
path<-"evaluation/results/rank-study.csv"
done<-if(file.exists(path)) read.csv(path) else data.frame()
completed<-if(nrow(done)) done$job else integer()
for(i in setdiff(seq_len(nrow(jobs)),completed)) {
  a<-jobs[i,]; set.seed(a$data_seed); dat<-simulate_case(a$scenario,a$n1,a$n2,a$d)
  set.seed(a$fit_seed)
  f<-compare_pt(dat$x,dat$y,transform="rank",prior=pt_prior(max_depth=6,split_grid=7,min_node=20),
    control=pt_control(particles=64,keep_data=FALSE,effect_draws=256))
  dg<-diagnostics(f); r<-regions(f,prob=0)
  loc<-if(a$scenario=="local") localization(f) else c(top_iou=NA,best5_iou=NA)
  row<-data.frame(job=i,rank_score=-log(max(f$prob_null,.Machine$double.xmin)),
    rank_probability=f$prob_difference,rank_seconds=dg$elapsed_seconds,
    rank_ancestors=tail(dg$trace$root_ancestors,1),rank_top_variable=r$variable[1],
    rank_top_iou=loc[1],rank_best5_iou=loc[2],rank_depth_fraction=dg$depth_limited_fraction)
  done<-rbind(done,row); write.csv(done,path,row.names=FALSE)
  if(i%%20==0) cat(sprintf("%s rank %d/%d d=%d %s: %.2fs P(diff)=%.3f\n",
    format(Sys.time(),"%H:%M:%S"),i,nrow(jobs),a$d,a$scenario,dg$elapsed_seconds,f$prob_difference))
}
cat("Rank transform paired study complete.\n")
