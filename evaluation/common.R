source("scripts/load-local.R")

simulate_case <- function(scenario, n1=300L, n2=n1, d=2L) {
  heavy <- grepl("heavy", scenario)
  gen <- if (heavy) function(n) matrix(rt(n*d,3),ncol=d) else function(n) matrix(rnorm(n*d),ncol=d)
  x <- gen(n1); y <- gen(n2)
  if (scenario %in% c("location","unbalanced_location")) y[,1] <- y[,1]+.5
  if (scenario == "heavy_location") y[,1] <- y[,1]+.8
  if (scenario == "dispersion") y[,1] <- 1.7*y[,1]
  if (scenario == "local") y[y[,1]<qnorm(.35),2] <- y[y[,1]<qnorm(.35),2]+1.3
  if (scenario == "correlation") y[,2] <- .8*y[,1]+sqrt(1-.8^2)*y[,2]
  colnames(x) <- colnames(y) <- paste0("X",seq_len(d))
  list(x=x,y=y)
}

baseline_scores <- function(x,y,mrs=FALSE) {
  n<-nrow(x); m<-nrow(y); z<-rbind(x,y); g<-c(rep(0L,n),rep(1L,m))
  # Gaussian MMD U-statistic; pooled scaling and median distance bandwidth.
  start<-proc.time()[3]; zsc<-scale(z); D<-as.matrix(dist(zsc))^2
  bw<-median(D[upper.tri(D)]); K<-exp(-D/(2*bw)); diag(K)<-0
  mmd<-sum(K[seq_len(n),seq_len(n)])/(n*(n-1))+
    sum(K[n+seq_len(m),n+seq_len(m)])/(m*(m-1))-2*sum(K[seq_len(n),n+seq_len(m)])/(n*m)
  mmd_time<-proc.time()[3]-start
  # Held-out classifier comparison; train and test are disjoint and balanced.
  start<-proc.time()[3]; size<-min(n,m)
  a<-sample.int(n,size); b<-n+sample.int(m,size); nt<-floor(.6*size)
  train<-c(a[seq_len(nt)],b[seq_len(nt)]); test<-c(a[-seq_len(nt)],b[-seq_len(nt)])
  frame<-data.frame(label=factor(g),z)
  cf<-rpart::rpart(label~.,data=frame[train,],method="class",
    control=rpart::rpart.control(cp=.001,minsplit=20,maxdepth=8,xval=0))
  acc<-mean(predict(cf,frame[test,],type="class")==frame$label[test])
  classifier_time<-proc.time()[3]-start
  mrs_score<-mrs_time<-NA_real_
  if(mrs) {
    start<-proc.time()[3]; rng<-t(apply(z,2,range)); w<-rng[,2]-rng[,1]
    support<-cbind(rng[,1]-.01*w,rng[,2]+.01*w)
    ans<-MRS::mrs(z,g+1L,Omega=support,K=6,alpha=.5,return_tree=FALSE)
    mrs_score<--log(max(ans$PostGlobNull,.Machine$double.xmin)); mrs_time<-proc.time()[3]-start
  }
  c(mmd=unname(mmd),classifier=unname(acc),mrs=unname(mrs_score),mmd_seconds=unname(mmd_time),
    classifier_seconds=unname(classifier_time),mrs_seconds=unname(mrs_time))
}

node_bounds <- function(fit,node) {
  t<-fit$tree; lo<-fit$support[,1]; hi<-fit$support[,2]
  while(t$parent_id[node]>0) {
    pa<-t$parent_id[node]; d<-t$variable_id[pa]
    if(t$left[pa]==node) hi[d]<-min(hi[d],t$threshold[pa]) else lo[d]<-max(lo[d],t$threshold[pa])
    node<-pa
  }
  list(lo=lo,hi=hi)
}
# Jaccard under a standard normal reference for true gate X1 < qnorm(.35).
# "Best top five" is explicitly an oracle diagnostic within the returned list.
localization <- function(fit) {
  r<-regions(fit,prob=0,min_n=30); r<-r[r$node_id!=1,]
  if(!nrow(r)) return(c(top_iou=NA,best5_iou=NA))
  iou<-vapply(head(r$node_id,5),function(node) {
    b<-node_bounds(fit,node); mass<-prod(pnorm(b$hi)-pnorm(b$lo))
    hi<-b$hi; hi[1]<-min(hi[1],qnorm(.35))
    intersection<-prod(pmax(0,pnorm(hi)-pnorm(b$lo)))
    intersection/(mass+.35-intersection)
  },numeric(1))
  c(top_iou=iou[1],best5_iou=max(iou))
}
