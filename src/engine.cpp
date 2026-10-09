#include <Rcpp.h>
#include <array>
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace Rcpp;
namespace pt {
using A3 = std::array<double,3>;
const double neginf = -std::numeric_limits<double>::infinity();
double lse(const A3& a) {
  double m = *std::max_element(a.begin(),a.end());
  if (m == neginf) return m;
  return m + std::log(std::exp(a[0]-m)+std::exp(a[1]-m)+std::exp(a[2]-m));
}
double lsev(const std::vector<double>& a) {
  double m = *std::max_element(a.begin(),a.end()), s=0;
  if (!std::isfinite(m)) Rcpp::stop("Non-finite particle or proposal weights.");
  for (double v:a) s += std::exp(v-m);
  return m+std::log(s);
}
struct Prior {
  int depth, grid, min_node;
  bool spike;
  double gamma, beta, prune, precision, eta;
  A3 root;
  std::vector<std::array<A3,3>> trans;
  explicit Prior(const List& p) {
    depth=as<int>(p["max_depth"]); grid=as<int>(p["split_grid"]);
    min_node=as<int>(p["min_node"]); spike=as<bool>(p["midpoint_spike"]); gamma=as<double>(p["gamma"]);
    beta=as<double>(p["beta"]); prune=as<double>(p["prune"]);
    precision=as<double>(p["precision"]); eta=as<double>(p["eta"]);
    root={std::log((1-prune)*gamma),std::log((1-prune)*(1-gamma)),std::log(prune)};
    for(int k=0;k<=depth;++k) {
      double a=gamma*std::pow(beta,-k), b=gamma*std::pow(2.,-k);
      trans.push_back({A3{std::log((1-prune)*a),std::log((1-prune)*(1-a)),std::log(prune)},
                       A3{std::log((1-prune)*b),std::log((1-prune)*(1-b)),std::log(prune)},
                       A3{neginf,neginf,0.}});
    }
  }
};
struct Node {
  int parent=-1,left=-1,right=-1,depth=0,dim=-1,n1=0,n2=0;
  bool numerical_stop=false, midpoint=false;
  double cut=0, fraction=.5;
  A3 emission{0,0,0}, up{0,0,0}, message{0,0,0};
  std::shared_ptr<std::vector<int>> obs;
};
struct Tree { std::vector<Node> nodes; size_t cursor=0; double logprior=0; };
A3 message(const A3& u, int depth, const Prior& p) {
  A3 out;
  for(int v=0;v<3;++v) {
    A3 tmp;
    for(int w=0;w<3;++w) tmp[w]=p.trans[depth][v][w]+u[w];
    out[v]=lse(tmp);
  }
  return out;
}
double evidence(const Tree& t,const Prior& p) {
  A3 tmp; for(int v=0;v<3;++v) tmp[v]=p.root[v]+t.nodes[0].up[v];
  return lse(tmp);
}
void refresh(Tree& t,int id,const Prior& p) {
  while(id>=0) {
    Node& a=t.nodes[id]; a.up=a.emission;
    if(a.left>=0) for(int v=0;v<3;++v)
      a.up[v]+=t.nodes[a.left].message[v]+t.nodes[a.right].message[v];
    a.message=message(a.up,a.depth,p); id=a.parent;
  }
}
// State distribution at a node with that node's subtree removed.
A3 outside(const Tree& t,int id,const Prior& p) {
  std::vector<int> path;
  for(int j=id;t.nodes[j].parent>=0;j=t.nodes[j].parent) path.push_back(j);
  A3 out=p.root;
  for(auto it=path.rbegin();it!=path.rend();++it) {
    int child=*it; const Node& par=t.nodes[t.nodes[child].parent];
    int sib=par.left==child?par.right:par.left;
    A3 tmp;
    for(int v=0;v<3;++v) tmp[v]=out[v]+par.emission[v]+t.nodes[sib].message[v];
    double z=lse(tmp); for(double& v:tmp) v-=z;
    for(int w=0;w<3;++w) {
      A3 terms;
      for(int v=0;v<3;++v) terms[v]=tmp[v]+p.trans[t.nodes[child].depth][v][w];
      out[w]=lse(terms);
    }
  }
  double z=lse(out); for(double& v:out) v-=z;
  return out;
}
A3 emission(int l1,int l2,int n1,int n2,double r,const Prior& p) {
  double a=p.precision*r,b=p.precision*(1-r), base=R::lbeta(a,b);
  double vol=-(l1+l2)*std::log(r)-(n1+n2-l1-l2)*std::log1p(-r);
  double diff=R::lbeta(a+l1,b+n1-l1)+R::lbeta(a+l2,b+n2-l2)-2*base+vol;
  double shared=R::lbeta(a+l1+l2,b+n1+n2-l1-l2)-base+vol;
  return {diff,shared,shared};
}
void bounds(const Tree& t,int id,int d,std::vector<double>& lo,std::vector<double>& hi) {
  lo.assign(d,0); hi.assign(d,1);
  for(int j=id;t.nodes[j].parent>=0;j=t.nodes[j].parent) {
    const Node& par=t.nodes[t.nodes[j].parent];
    if(par.left==j) hi[par.dim]=std::min(hi[par.dim],par.cut);
    else lo[par.dim]=std::max(lo[par.dim],par.cut);
  }
}
struct Option { int dim,l1,l2; bool midpoint; double r,cut,logprior; A3 em; };
struct Proposal { int id=-1; double increment=0; std::vector<Option> options; std::vector<double> cdf; };
Proposal propose(Tree& t,const NumericMatrix& x,const IntegerVector& group,const Prior& p) {
  Proposal prop;
  std::vector<double> lo,hi;
  while(t.cursor<t.nodes.size()) {
    Node& a=t.nodes[t.cursor];
    if(a.depth<p.depth && a.n1+a.n2>=p.min_node) {
      bounds(t,t.cursor,x.ncol(),lo,hi);
      bool safe=true;
      for(int j=0;j<x.ncol();++j)
        if(hi[j]-lo[j] <= 4*(p.grid+1)*std::numeric_limits<double>::epsilon()) safe=false;
      if(safe) break;
      a.numerical_stop=true;
    }
    ++t.cursor;
  }
  if(t.cursor==t.nodes.size()) return prop;
  prop.id=t.cursor;
  const Node& a=t.nodes[prop.id]; int d=x.ncol(), L=p.grid;
  A3 pi=outside(t,prop.id,p);
  std::vector<double> lp(L),logq;
  for(int k=0;k<L;++k) lp[k]=-p.eta*(a.n1+a.n2)*std::abs((k+1.)/(L+1.)-.5);
  double z=lsev(lp); for(double& v:lp) v-=z+std::log((double)d);
  // Awaya and Ma (2024), Appendix A: R(A)=1 iff L(A)=1/2, and R is absorbing.
  // A child of a midpoint split is therefore divided at the midpoint; otherwise
  // the location prior is exactly (2), whose midpoint mass is r(A).
  const int mid=(L-1)/2;
  bool locked=p.spike && a.parent>=0 && t.nodes[a.parent].midpoint;
  if(locked) for(int k=0;k<L;++k) lp[k]=k==mid?-std::log((double)d):neginf;
  prop.options.reserve(d*L); logq.reserve(d*L);
  for(int j=0;j<d;++j) {
    if(j%8==0) checkUserInterrupt();
    std::vector<int> h1(L+1,0),h2(L+1,0);
    std::vector<double> cuts(L);
    for(int k=0;k<L;++k) cuts[k]=lo[j]+(k+1.)/(L+1.)*(hi[j]-lo[j]);
    for(int row:*a.obs) {
      double value=x(row,j),scaled=(value-lo[j])/(hi[j]-lo[j])*(L+1.);
      int bin=std::max(0,std::min(L,(int)std::floor(scaled)));
      // Correct floating-point boundary ambiguity using actual cut values.
      while(bin>0 && value<=cuts[bin-1]) --bin;
      while(bin<L && value>cuts[bin]) ++bin;
      (group[row]==0?h1:h2)[bin]++;
    }
    int l1=0,l2=0;
    for(int k=0;k<L;++k) {
      l1+=h1[k]; l2+=h2[k]; double r=(k+1.)/(L+1.);
      if(lp[k]==neginf) continue;
      A3 em=emission(l1,l2,a.n1,a.n2,r,p), terms;
      for(int v=0;v<3;++v) terms[v]=pi[v]+em[v];
      prop.options.push_back({j,l1,l2,k==mid,r,cuts[k],lp[k],em});
      logq.push_back(lp[k]+lse(terms));
    }
  }
  prop.increment=lsev(logq); double cum=0;
  for(double v:logq) { cum+=std::exp(v-prop.increment); prop.cdf.push_back(cum); }
  prop.cdf.back()=1; return prop;
}
void split(Tree& t,const Proposal& prop,int choice,const NumericMatrix& x,
           const IntegerVector& group,const Prior& p) {
  const Option& o=prop.options[choice]; int id=prop.id, left=t.nodes.size();
  Node a=t.nodes[id], l,r;
  l.parent=r.parent=id; l.depth=r.depth=a.depth+1;
  l.n1=o.l1; l.n2=o.l2; r.n1=a.n1-o.l1; r.n2=a.n2-o.l2;
  l.obs=std::make_shared<std::vector<int>>(); r.obs=std::make_shared<std::vector<int>>();
  l.obs->reserve(l.n1+l.n2); r.obs->reserve(r.n1+r.n2);
  for(int row:*a.obs) (x(row,o.dim)<=o.cut?*l.obs:*r.obs).push_back(row);
  if((int)l.obs->size()!=l.n1+l.n2) stop("Internal split count mismatch.");
  a.left=left; a.right=left+1; a.dim=o.dim; a.cut=o.cut; a.fraction=o.r; a.midpoint=o.midpoint;
  a.emission=o.em; a.obs.reset(); t.nodes[id]=a;
  t.nodes.push_back(l); t.nodes.push_back(r); ++t.cursor;
  t.logprior+=o.logprior; refresh(t,id,p);
}
double null_probability(const Tree& t,const Prior& p) {
  std::vector<A3> msg(t.nodes.size()); A3 rootup{0,0,0};
  for(int j=(int)t.nodes.size()-1;j>=0;--j) {
    const Node& a=t.nodes[j]; A3 u{0,0,0};
    if(a.left>=0) {
      u[0]=neginf;
      for(int v=1;v<3;++v) u[v]=a.emission[v]+msg[a.left][v]+msg[a.right][v];
    }
    msg[j]=message(u,a.depth,p); if(j==0) rootup=u;
  }
  for(int v=0;v<3;++v) rootup[v]+=p.root[v];
  return std::min(1.,std::exp(lse(rootup)-evidence(t,p)));
}
DataFrame table(const Tree& t,const Prior& p) {
  int n=t.nodes.size(); IntegerVector id(n),parent(n),left(n),right(n),depth(n),dim(n),n1(n),n2(n);
  NumericVector cut(n),frac(n),pd(n),a1(n),b1(n),a2(n),b2(n),signed_eff(n);
  CharacterVector reason(n);
  for(int j=0;j<n;++j) {
    const Node& a=t.nodes[j]; id[j]=j+1; parent[j]=a.parent+1; left[j]=a.left+1;
    right[j]=a.right+1; depth[j]=a.depth; n1[j]=a.n1;n2[j]=a.n2;
    bool internal=a.left>=0;
    dim[j]=internal?a.dim+1:NA_INTEGER; cut[j]=internal?a.cut:NA_REAL; frac[j]=internal?a.fraction:NA_REAL;
    reason[j]=a.left>=0?"split":(a.numerical_stop?"numerical precision":(a.depth>=p.depth?"max_depth":"min_node"));
    pd[j]=a1[j]=b1[j]=a2[j]=b2[j]=signed_eff[j]=NA_REAL;
    if(a.left>=0) {
      A3 post=outside(t,j,p); for(int v=0;v<3;++v) post[v]+=a.up[v];
      pd[j]=std::exp(post[0]-lse(post));
      const Node& l=t.nodes[a.left]; const Node& r=t.nodes[a.right];
      a1[j]=p.precision*a.fraction+l.n1; b1[j]=p.precision*(1-a.fraction)+r.n1;
      a2[j]=p.precision*a.fraction+l.n2; b2[j]=p.precision*(1-a.fraction)+r.n2;
      signed_eff[j]=pd[j]*(R::digamma(a1[j])-R::digamma(b1[j])-R::digamma(a2[j])+R::digamma(b2[j]));
    }
  }
  return DataFrame::create(_["node_id"]=id,_["parent_id"]=parent,_["left"]=left,_["right"]=right,
      _["depth"]=depth,_["variable_id"]=dim,_["cut_unit"]=cut,_["fraction"]=frac,
      _["n1"]=n1,_["n2"]=n2,_["prob_difference"]=pd,_["effect"]=signed_eff,
      _["a1"]=a1,_["b1"]=b1,_["a2"]=a2,_["b2"]=b2,_["stop_reason"]=reason);
}
struct Particle { std::shared_ptr<Tree> tree; double logw; int ancestor; };
size_t payload(const std::vector<Particle>& ps) {
  std::unordered_set<const Tree*> ts; std::unordered_set<const std::vector<int>*> blocks;
  size_t bytes=0;
  for(const auto& q:ps) if(ts.insert(q.tree.get()).second) {
    bytes+=sizeof(Tree)+q.tree->nodes.capacity()*sizeof(Node);
    for(const auto& n:q.tree->nodes) if(n.obs && blocks.insert(n.obs.get()).second)
      bytes+=sizeof(std::vector<int>)+sizeof(int)*n.obs->capacity();
  }
  return bytes;
}
} // namespace pt

// [[Rcpp::export]]
List pt_smc_cpp(NumericMatrix x,IntegerVector group,List prior,List control) {
  using namespace pt;
  Prior p(prior); int M=as<int>(control["particles"]);
  double threshold=as<double>(control["ess_threshold"]);
  double kappa=as<double>(control["resample_power"]);
  bool keep=as<bool>(control["keep_particles"]);
  auto initial=std::make_shared<Tree>(); Node root;
  root.obs=std::make_shared<std::vector<int>>(); root.obs->reserve(x.nrow());
  for(int i=0;i<x.nrow();++i) { root.obs->push_back(i); if(group[i]==0) ++root.n1; else ++root.n2; }
  initial->nodes.push_back(root);
  std::vector<Particle> ps;
  for(int i=0;i<M;++i) ps.push_back({initial,-std::log((double)M),i});
  initial.reset(); root.obs.reset();
  std::vector<double> ess_trace; std::vector<int> unique_trace,ancestor_trace,resampled_trace;
  double logZ=0; size_t peak=payload(ps); int rounds=0,resamples=0,proposals=0;
  while(true) {
    checkUserInterrupt();
    std::unordered_map<Tree*,size_t> lookup;
    std::vector<std::vector<int>> groups;
    for(int i=0;i<M;++i) {
      Tree* key=ps[i].tree.get(); auto hit=lookup.find(key);
      if(hit==lookup.end()) { lookup[key]=groups.size(); groups.push_back({i}); }
      else groups[hit->second].push_back(i);
    }
    bool changed=false;
    for(const auto& ids:groups) {
      Proposal q=propose(*ps[ids[0]].tree,x,group,p);
      if(q.id<0) continue;
      changed=true; ++proposals;
      if(ids.size()==1) {
        int i=ids[0]; int c=std::lower_bound(q.cdf.begin(),q.cdf.end(),R::runif(0,1))-q.cdf.begin();
        split(*ps[i].tree,q,c,x,group,p); ps[i].logw+=q.increment;
      } else {
        auto source=ps[ids[0]].tree;
        std::unordered_map<int,std::shared_ptr<Tree>> children;
        for(int i:ids) {
          int c=std::lower_bound(q.cdf.begin(),q.cdf.end(),R::runif(0,1))-q.cdf.begin();
          auto hit=children.find(c);
          if(hit==children.end()) {
            auto child=std::make_shared<Tree>(*source); split(*child,q,c,x,group,p);
            children[c]=child; ps[i].tree=child;
          } else ps[i].tree=hit->second;
          ps[i].logw+=q.increment;
        }
      }
    }
    if(!changed) break;
    ++rounds; std::vector<double> logw; for(const auto& q:ps) logw.push_back(q.logw);
    double z=lsev(logw),sq=0; logZ+=z;
    for(auto& q:ps) { q.logw-=z; sq+=std::exp(2*q.logw); }
    double ess=1/sq;
    bool resample=ess<threshold*M;
    if(resample) {
      // Section 3.1 (after Lu et al., 2013): ancestors are drawn with
      // probabilities a proportional to W^kappa and receive weights W/a.
      ++resamples; std::vector<Particle> next; next.reserve(M);
      std::vector<double> la(M);
      for(int i=0;i<M;++i) la[i]=kappa*ps[i].logw;
      double za=lsev(la); for(double& v:la) v-=za;
      double u=R::runif(0,1)/M,cum=std::exp(la[0]); int j=0;
      for(int i=0;i<M;++i) {
        double target=u+(double)i/M;
        while(j<M-1 && target>cum) cum+=std::exp(la[++j]);
        next.push_back({ps[j].tree,ps[j].logw-la[j]-std::log((double)M),ps[j].ancestor});
      }
      ps.swap(next);
      // The resampled weights average to one in expectation; keep their sum in
      // the evidence estimate and renormalize. For kappa = 1 the sum is one.
      std::vector<double> lw; for(const auto& q:ps) lw.push_back(q.logw);
      double zr=lsev(lw); logZ+=zr; for(auto& q:ps) q.logw-=zr;
    }
    std::unordered_set<Tree*> unique; std::unordered_set<int> ancestors;
    for(const auto& q:ps) { unique.insert(q.tree.get()); ancestors.insert(q.ancestor); }
    ess_trace.push_back(ess); unique_trace.push_back(unique.size());
    ancestor_trace.push_back(ancestors.size()); resampled_trace.push_back(resample);
    // Measured at rounds 1, 2, 4, 8, ... and at the end, to keep its cost small.
    if((rounds&(rounds-1))==0) peak=std::max(peak,payload(ps));
  }
  peak=std::max(peak,payload(ps));
  std::unordered_map<Tree*,double> nulls;
  double global=0,bestscore=neginf; int best=0;
  NumericVector weights(M),pn(M),scores(M); List particles(keep?M:0);
  for(int i=0;i<M;++i) {
    auto& t=*ps[i].tree;
    if(!nulls.count(&t)) nulls[&t]=null_probability(t,p);
    pn[i]=nulls[&t]; weights[i]=std::exp(ps[i].logw); global+=weights[i]*pn[i];
    scores[i]=t.logprior+evidence(t,p);
    if(scores[i]>bestscore) { bestscore=scores[i]; best=i; }
    if(keep) particles[i]=table(t,p);
  }
  DataFrame trace=DataFrame::create(_["step"]=seq_len(rounds),_["ess"]=ess_trace,
      _["unique_trees"]=unique_trace,_["root_ancestors"]=ancestor_trace,_["resampled"]=resampled_trace);
  return List::create(_["tree"]=table(*ps[best].tree,p),_["prob_null"]=std::min(1.,global),
      _["weights"]=weights,_["particle_null"]=pn,_["particle_scores"]=scores,_["trace"]=trace,
      _["log_evidence"]=logZ,_["resamples"]=resamples,_["proposals"]=proposals,
      _["peak_particle_payload_bytes"]=(double)peak,_["particles"]=particles);
}

// Internal verification entry: a fixed tree, independently enumerated in R tests.
// [[Rcpp::export]]
List pt_fixed_cpp(IntegerMatrix children,IntegerMatrix counts,NumericVector fraction,
                  IntegerVector depth,List prior) {
  using namespace pt;
  Prior p(prior); Tree t; int n=children.nrow(); t.nodes.resize(n);
  if(counts.nrow()!=n || counts.ncol()!=2 || children.ncol()!=2 || fraction.size()!=n || depth.size()!=n)
    stop("Invalid fixed tree dimensions.");
  for(int j=0;j<n;++j) {
    Node& a=t.nodes[j]; a.left=children(j,0)-1;a.right=children(j,1)-1;
    a.depth=depth[j];a.fraction=fraction[j];a.n1=counts(j,0);a.n2=counts(j,1);
    if(a.depth<0 || a.depth>p.depth) stop("Invalid fixed tree depth.");
    if(a.left>=0) {
      if(a.left<=j || a.right<=j || a.left>=n || a.right>=n) stop("Invalid fixed tree child.");
      t.nodes[a.left].parent=t.nodes[a.right].parent=j;
      a.emission=emission(counts(a.left,0),counts(a.left,1),a.n1,a.n2,a.fraction,p);
    }
  }
  for(int j=n-1;j>=0;--j) refresh(t,j,p);
  return List::create(_["log_evidence"]=evidence(t,p),_["prob_null"]=null_probability(t,p),_["tree"]=table(t,p));
}
