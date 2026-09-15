source("scripts/load-local.R")
a<-read.csv("evaluation/results/study.csv")
b<-read.csv("evaluation/results/rank-study.csv")
stopifnot(nrow(a)==1220L,nrow(b)==1220L)
z<-merge(a,b,by="job")
z$family<-ifelse(grepl("heavy",z$scenario),"heavy",ifelse(grepl("unbalanced",z$scenario),"unbalanced","normal"))
methods<-c(pt_score="PT range",rank_score="PT rank",mmd="MMD",classifier="Classifier tree",mrs="MRS")
auc<-function(al,nu) {
  n<-length(al); m<-length(nu)
  (sum(rank(c(al,nu),ties.method="average")[seq_len(n)])-n*(n+1)/2)/(n*m)
}
rate<-function(v) {
  t<-binom.test(sum(v),length(v))
  c(mean(v),unname(t$conf.int))
}
rows<-list(); set.seed(84520)
for(d in unique(z$d)) for(s in unique(z$scenario[z$d==d & !grepl("null",z$scenario)])) {
  alt<-z[z$d==d & z$scenario==s,]
  nu<-z[z$d==d & grepl("null",z$scenario) & z$family==alt$family[1],]
  nu<-nu[order(nu$rep),]; ncal<-floor(nrow(nu)/2)
  cal<-nu[seq_len(ncal),]; val<-nu[-seq_len(ncal),]
  for(key in names(methods)) {
    if(all(is.na(alt[[key]]))) next
    k<-ceiling((ncal+1)*.95)
    threshold<-if(k>ncal) Inf else sort(cal[[key]])[k]
    pow<-rate(alt[[key]]>threshold); fp<-rate(val[[key]]>threshold)
    A<-auc(alt[[key]],val[[key]])
    boots<-replicate(400,auc(sample(alt[[key]],replace=TRUE),sample(val[[key]],replace=TRUE)))
    ci<-unname(quantile(boots,c(.025,.975)))
    bayes<-key%in%c("pt_score","rank_score","mrs")
    rows[[length(rows)+1]]<-data.frame(d=d,scenario=s,method=methods[key],n_alt=nrow(alt),
      n_cal=ncal,n_null_validation=nrow(val),auc=A,auc_lower=ci[1],auc_upper=ci[2],
      calibrated_power=pow[1],power_lower=pow[2],power_upper=pow[3],
      validation_fpr=fp[1],fpr_lower=fp[2],fpr_upper=fp[3],threshold=threshold,
      posterior95_power=if(bayes)mean(alt[[key]] > -log(.05)) else NA,
      posterior95_fpr=if(bayes)mean(nu[[key]] > -log(.05)) else NA)
  }
}
comparison<-do.call(rbind,rows)
write.csv(comparison,"evaluation/results/comparison.csv",row.names=FALSE)
null_rates<-list()
for(d in unique(z$d)) for(s in unique(z$scenario[z$d==d & grepl("null",z$scenario)])) {
  nu<-z[z$d==d & z$scenario==s,]
  for(key in c("pt_score","rank_score")) {
    flags<-nu[[key]] > -log(.05); ci<-rate(flags)
    null_rates[[length(null_rates)+1]]<-data.frame(d=d,scenario=s,method=methods[key],
      n=nrow(nu),exceedances=sum(flags),rate=ci[1],lower=ci[2],upper=ci[3])
  }
}
write.csv(do.call(rbind,null_rates),"evaluation/results/posterior-threshold-null.csv",row.names=FALSE)
local<-z[z$scenario=="local",]
local_summary<-aggregate(cbind(top_iou,best5_iou,rank_top_iou,rank_best5_iou)~d,local,mean)
write.csv(local_summary,"evaluation/results/localization.csv",row.names=FALSE)
timing<-aggregate(cbind(seconds,rank_seconds,payload_mb,root_ancestors,rank_ancestors)~d+scenario,z,median)
write.csv(timing,"evaluation/results/study-diagnostics.csv",row.names=FALSE)

dir.create("evaluation/figures",showWarnings=FALSE)
colors<-c("PT range"="#B65D35","PT rank"="#176B87","MMD"="#497A43","Classifier tree"="#8C6BA6","MRS"="#666666")
png("evaluation/figures/auc.png",1800,1100,res=140)
par(mfrow=c(2,2),mar=c(5,4,3,1),oma=c(1,0,2,0))
for(d in c(2,10,50,100)) {
  zz<-comparison[comparison$d==d & comparison$scenario %in% c("location","dispersion","local","correlation"),]
  scenarios<-unique(zz$scenario); ms<-unique(zz$method)
  h<-sapply(scenarios,function(s) zz$auc[match(paste(ms,s),paste(zz$method,zz$scenario))])
  bp<-barplot(h,beside=TRUE,ylim=c(0,1.25),col=colors[ms],names.arg=scenarios,
    ylab="AUC against held-out null samples",main=paste("Dimension",d),border=NA,cex.names=.85,yaxt="n")
  axis(2,at=seq(0,1,.2))
  low<-sapply(scenarios,function(s) zz$auc_lower[match(paste(ms,s),paste(zz$method,zz$scenario))])
  high<-sapply(scenarios,function(s) zz$auc_upper[match(paste(ms,s),paste(zz$method,zz$scenario))])
  segments(bp,low,bp,high,col="grey25",lwd=1)
  abline(h=.5,lty=2,col="grey50")
  legend("top",ms,fill=colors[ms],bty="n",cex=.68,horiz=TRUE)
}
mtext("Finite PT models: coordinate choice and dimension",outer=TRUE,cex=1.25)
dev.off()

e<-read.csv("evaluation/results/exact-convergence.csv")
png("evaluation/figures/exact-convergence.png",1200,650,res=130)
par(mfrow=c(1,2),mar=c(4,4,3,1))
for(v in c("probability_error","total_variation")) {
  yy<-if(v=="probability_error") abs(e[[v]]) else e[[v]]
  plot(e$particles,yy,log="x",pch=ifelse(e$ess_threshold==0,1,16),
    col=ifelse(e$ess_threshold==0,"#B65D35","#176B87"),xlab="Particles",
    ylab=if(v=="probability_error") "Absolute P(null) error" else "Tree posterior total variation",
    main=if(v=="probability_error") "Global null probability" else "Complete tree posterior")
  legend("topright",c("No resampling","Resampling"),pch=c(1,16),col=c("#B65D35","#176B87"),bty="n",cex=.8)
}
dev.off()
cat("Summaries and figures written.\n")
