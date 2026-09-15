// Standalone reproducible kernel comparison, not a benchmark of the legacy package.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <stdexcept>
using Clock=std::chrono::steady_clock;
volatile double sink=0;
template<class F> double time_ms(F f,int reps) {
  auto start=Clock::now();for(int i=0;i<reps;++i) f();
  return std::chrono::duration<double,std::milli>(Clock::now()-start).count()/reps;
}
int main() {
  std::mt19937 rng(713);std::uniform_real_distribution<double> u(0,1);
  std::cout<<"kernel,size,reference_ms,optimized_ms,speedup\n";
  for(int M:{256,1024,4096}) {
    std::vector<double> lw(M),a(M),b(M);for(double& v:lw)v=4*u(rng)-2;
    auto old=[&](){for(int i=0;i<M;++i){double s=0;for(int j=0;j<M;++j)s+=std::exp(lw[j]-lw[i]);a[i]=1/s;}sink=a[M/2];};
    auto opt=[&](){double mx=*std::max_element(lw.begin(),lw.end()),s=0;for(int i=0;i<M;++i){b[i]=std::exp(lw[i]-mx);s+=b[i];}for(double& v:b)v/=s;sink=b[M/2];};
    old();opt();for(int i=0;i<M;++i)if(std::abs(a[i]-b[i])>1e-13)throw std::runtime_error("normalization mismatch");
    double t0=time_ms(old,M<4096?20:3),t1=time_ms(opt,200);
    std::cout<<"normalization,"<<M<<","<<t0<<","<<t1<<","<<t0/t1<<"\n";
  }
  for(int n:{2000,20000,200000}) {
    int L=31,d=10;std::vector<int> bins(n*d),a(d*L),b(d*L),h(d*(L+1));
    for(int& v:bins)v=(int)(u(rng)*(L+1));
    auto old=[&](){std::fill(a.begin(),a.end(),0);for(int j=0;j<d;++j)for(int i=0;i<n;++i)for(int k=bins[j*n+i];k<L;++k)++a[j*L+k];sink=a[0];};
    auto opt=[&](){std::fill(h.begin(),h.end(),0);for(int j=0;j<d;++j){for(int i=0;i<n;++i)++h[j*(L+1)+bins[j*n+i]];int s=0;for(int k=0;k<L;++k){s+=h[j*(L+1)+k];b[j*L+k]=s;}}sink=b[0];};
    old();opt();if(a!=b)throw std::runtime_error("histogram mismatch");
    double t0=time_ms(old,5),t1=time_ms(opt,5);
    std::cout<<"candidate_counts_d10_L31,"<<n<<","<<t0<<","<<t1<<","<<t0/t1<<"\n";
  }
}
