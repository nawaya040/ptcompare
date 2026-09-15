# Independent scalar R implementation. Enumerates all internal-node states.
ref_transition <- function(depth, p) {
  a <- p$gamma * p$beta^(-depth); b <- p$gamma * 2^(-depth)
  rbind(c((1-p$prune)*a, (1-p$prune)*(1-a), p$prune),
        c((1-p$prune)*b, (1-p$prune)*(1-b), p$prune), c(0,0,1))
}
ref_lse <- function(v) { m <- max(v); if (is.infinite(m)) m else m + log(sum(exp(v-m))) }
ref_fixed <- function(t, p) {
  ids <- which(t$left > 0)
  if (!length(ids)) return(list(log_evidence=0, prob_null=1, prob_diff=rep(NA_real_,nrow(t))))
  states <- as.matrix(expand.grid(rep(list(1:3), length(ids))))
  logterms <- numeric(nrow(states))
  for (h in seq_len(nrow(states))) {
    state <- integer(nrow(t)); state[ids] <- states[h, ]
    logterm <- 0
    for (j in ids) {
      s <- state[j]; pa <- t$parent_id[j]
      pr <- if (pa == 0) ref_transition(0, p)[1, s] else ref_transition(t$depth[j], p)[state[pa],s]
      l <- t$left[j]; r <- t$right[j]; f <- t$fraction[j]
      a <- p$precision*f; b <- p$precision*(1-f)
      if (s == 1) {
        ml <- lbeta(a+t$n1[l], b+t$n1[r]) + lbeta(a+t$n2[l], b+t$n2[r]) - 2*lbeta(a,b)
      } else ml <- lbeta(a+t$n1[l]+t$n2[l], b+t$n1[r]+t$n2[r])-lbeta(a,b)
      vol <- -(t$n1[l]+t$n2[l])*log(f) -(t$n1[r]+t$n2[r])*log1p(-f)
      logterm <- logterm+log(pr)+ml+vol
    }
    logterms[h] <- logterm
  }
  z <- ref_lse(logterms); w <- exp(logterms-z); pd <- rep(NA_real_, nrow(t))
  for (k in seq_along(ids)) pd[ids[k]] <- sum(w[states[, k] == 1])
  list(log_evidence=z, prob_null=sum(w[rowSums(states == 1) == 0]), prob_diff=pd)
}
ref_enumerate <- function(x, g, p) {
  grow <- function(rows, lo, hi, depth) {
    base <- list(rows=rows, depth=depth)
    if (depth >= p$max_depth || length(rows) < p$min_node)
      return(list(list(node=base, logprior=0)))
    fractions <- seq_len(p$split_grid)/(p$split_grid+1)
    lp <- -p$eta*length(rows)*abs(fractions-.5)
    lp <- lp-ref_lse(lp)-log(ncol(x))
    out <- list()
    for (j in seq_len(ncol(x))) for (k in seq_along(fractions)) {
      cut <- lo[j]+fractions[k]*(hi[j]-lo[j]); il <- rows[x[rows,j]<=cut]; ir <- setdiff(rows,il)
      hil <- hi; hil[j] <- cut; lor <- lo; lor[j] <- cut
      ls <- grow(il,lo,hil,depth+1); rs <- grow(ir,lor,hi,depth+1)
      for (l in ls) for (r in rs) out[[length(out)+1]] <- list(
        node=c(base,list(dim=j, fraction=fractions[k], cut=cut, left=l$node, right=r$node)),
        logprior=lp[k]+l$logprior+r$logprior)
    }
    out
  }
  flatten <- function(node) {
    queue <- list(node); parents <- 0L; out <- list(); i <- 1L
    while (i <= length(queue)) {
      a <- queue[[i]]; l <- r <- 0L
      if (!is.null(a$left)) {
        l <- length(queue)+1L; r <- l+1L
        queue <- c(queue,list(a$left,a$right)); parents <- c(parents,i,i)
      }
      out[[i]] <- data.frame(node_id=i,parent_id=parents[i],left=l,right=r,depth=a$depth,
        variable_id=if(is.null(a$dim)) 0L else a$dim,
        fraction=if(is.null(a$fraction)) .5 else a$fraction,
        n1=sum(g[a$rows]==0),n2=sum(g[a$rows]==1))
      i <- i+1L
    }
    do.call(rbind,out)
  }
  all <- grow(seq_len(nrow(x)),rep(0,ncol(x)),rep(1,ncol(x)),0)
  tables <- lapply(all, function(a) flatten(a$node))
  fixed <- lapply(tables, ref_fixed, p=p)
  lp <- vapply(all, `[[`, numeric(1), "logprior") + vapply(fixed, `[[`, numeric(1), "log_evidence")
  z <- ref_lse(lp); w <- exp(lp-z)
  list(tables=tables, weights=w, log_evidence=z,
       prob_null=sum(w*vapply(fixed, `[[`, numeric(1), "prob_null")))
}
tree_key <- function(t) paste(paste(t$variable_id, round(t$fraction,8), t$left,t$right,sep=":"),collapse="|")
