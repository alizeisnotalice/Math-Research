#include <algorithm>
#include <array>
#include <chrono>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

// Exact node ordering/equality; all evaluated values are binary64 exploration.
struct Ray { long long q; std::array<long long,6> r; };
struct Cone { long long det; std::array<Ray,6> rays; };
struct Node { long long n,d; };
bool eq(const Node&a,const Node&b){return (__int128)a.n*b.d==(__int128)b.n*a.d;}
bool lt(const Node&a,const Node&b){return (__int128)a.n*b.d<(__int128)b.n*a.d;}
long double value(const Node&a){return (long double)a.n/a.d;}
struct Kahan { long double s=0,c=0; void add(long double x){long double y=x-c,t=s+y;c=(t-s)-y;s=t;} };
long double jet(long double x,long double e,int k){
  if(!k)return std::atan(x/e); using C=std::complex<long double>; const C ii(0,1),p(e,x),m(e,-x);
  C z=(std::pow(ii,k)/std::pow(p,k)-std::pow(-ii,k)/std::pow(m,k))/((long double)(2*k)*ii);
  return ((k&1)?1.L:-1.L)*z.real();
}
long double dd5(std::array<Node,6> b,long double e){
  std::sort(b.begin(),b.end(),lt);long double t[6][6]{};for(int i=0;i<6;++i)t[i][0]=jet(value(b[i]),e,0);
  for(int o=1;o<6;++o)for(int i=o;i<6;++i)t[i][o]=eq(b[i],b[i-o])?jet(value(b[i]),e,o):(t[i][o-1]-t[i-1][o-1])/(value(b[i])-value(b[i-o]));
  return t[5][5];
}
int main(int argc,char**argv){
  std::string fan_path=argc>1?argv[1]:"fan_integer_input.txt", point_path=argc>2?argv[2]:"sobol_pilot8.txt";
  std::ifstream fan(fan_path),pts(point_path);int cn;fan>>cn;std::vector<Cone> cones(cn);
  for(auto&c:cones){fan>>c.det;for(auto&r:c.rays){fan>>r.q;for(auto&z:r.r)fan>>z;}}
  long long den;int pn;pts>>den>>pn;std::vector<std::array<long long,6>> xs(pn);for(auto&x:xs)for(auto&z:x)pts>>z;
  constexpr long double eps=0.001L;
  std::cout<<"{\n  \"status\":\"exploratory_binary64_H_1_over_1000\",\n  \"LDBL_MANT_DIG\":"<<LDBL_MANT_DIG<<",\n  \"sizeof_long_double\":"<<sizeof(long double)<<",\n  \"points\":[\n";
  for(int p=0;p<pn;++p){auto start=std::chrono::steady_clock::now();Kahan sum,abss;long long equal_events=0,terms=0;std::array<int,6> perm={0,1,2,3,4,5};
    do{for(int mask=0;mask<64;++mask){std::array<long long,6> y{};for(int i=0;i<6;++i)y[i]=((mask>>i)&1?-1:1)*xs[p][perm[i]];
      for(const auto&c:cones){std::array<Node,6>b{};long double alpha=c.det;for(int j=0;j<6;++j){long long n=0;for(int i=0;i<6;++i)n+=y[i]*c.rays[j].r[i];b[j]={n,den*c.rays[j].q};alpha/=c.rays[j].q;}
        auto z=b;std::sort(z.begin(),z.end(),lt);for(int j=1;j<6;++j)if(eq(z[j-1],z[j]))++equal_events;long double term=alpha*dd5(b,eps);sum.add(term);abss.add(std::fabs(term));++terms;}}
    } while(std::next_permutation(perm.begin(),perm.end()));long double h=sum.s/std::pow(2.L*std::numbers::pi_v<long double>,6),secs=std::chrono::duration<long double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"    {\"index\":"<<p<<",\"numerators\":[";for(int i=0;i<6;++i){if(i)std::cout<<",";std::cout<<xs[p][i];}std::cout<<"],\"denominator\":"<<den<<",\"raw_terms\":"<<terms<<",\"equal_adjacent_events\":"<<equal_events<<",\"signed_sum_binary64\":\""<<std::setprecision(21)<<sum.s<<"\",\"absolute_sum_binary64\":\""<<abss.s<<"\",\"cancellation_ratio\":\""<<(abss.s/std::fabs(sum.s))<<"\",\"H_epsilon_binary64\":\""<<h<<"\",\"elapsed_seconds\":\""<<secs<<"\"}"<<(p+1==pn?"\n":",\n");}
  std::cout<<"  ]\n}\n";
}
