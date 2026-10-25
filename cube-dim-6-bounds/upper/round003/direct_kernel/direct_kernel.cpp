#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

// Exact rational node classification; long double only evaluates accepted DD5.
struct Ray { long long q; std::array<long long,6> r; };
struct Cone { long long det; std::array<Ray,6> rays; };
struct Node { long long n, d; };
bool eq(const Node&a,const Node&b){return (__int128)a.n*b.d==(__int128)b.n*a.d;}
bool lt(const Node&a,const Node&b){return (__int128)a.n*b.d<(__int128)b.n*a.d;}
long double val(const Node&a){return (long double)a.n/a.d;}
struct Kahan{long double s=0,c=0; void add(long double x){long double y=x-c,t=s+y;c=(t-s)-y;s=t;}};
long double dd5(std::array<Node,6> v){
 std::sort(v.begin(),v.end(),lt); long double t[6][6]{};
 for(int i=0;i<6;i++)t[i][0]=v[i].n>0?1.L:-1.L;
 for(int o=1;o<6;o++)for(int i=o;i<6;i++)t[i][o]=(t[i][o-1]-t[i-1][o-1])/(val(v[i])-val(v[i-o]));
 return t[5][5];
}
int main(int argc,char**argv){
 int wanted=argc>1?std::stoi(argv[1]):4; std::string points_path=argc>2?argv[2]:"face_points.txt"; std::ifstream fi("fan_integer_input.txt"),pi(points_path); int cn;fi>>cn;std::vector<Cone> cs(cn);
 for(auto&c:cs){fi>>c.det;for(auto&r:c.rays){fi>>r.q;for(auto&z:r.r)fi>>z;}}
 long long den;int pn;pi>>den>>pn;std::vector<std::array<long long,6>> xs(pn);for(auto&x:xs)for(auto&z:x)pi>>z;wanted=std::min(wanted,pn);
 std::cout<<"{\n  \"status\":\"direct_generic_kernel_diagnostic\",\n  \"points\":[\n";
 for(int p=0;p<wanted;p++){
  auto start=std::chrono::steady_clock::now(); Kahan sum;long double abs_sum=0;long long zeros=0,repeats=0,ok=0;std::array<int,6> perm={0,1,2,3,4,5};
  do{for(int mask=0;mask<64;mask++){std::array<long long,6> tx{};for(int i=0;i<6;i++)tx[i]=((mask>>i)&1?-1:1)*xs[p][perm[i]];
   for(const auto&c:cs){std::array<Node,6> b{};bool bad=false;for(int j=0;j<6;j++){long long n=0;for(int i=0;i<6;i++)n+=tx[i]*c.rays[j].r[i];b[j]={n,den*c.rays[j].q};if(!n){zeros++;bad=true;}}
    auto sb=b;std::sort(sb.begin(),sb.end(),lt);for(int j=1;j<6;j++)if(eq(sb[j-1],sb[j])){repeats++;bad=true;}if(bad)continue;
    long double alpha=c.det;for(const auto&r:c.rays)alpha/=r.q;long double term=alpha*dd5(b);sum.add(term);abs_sum+=std::fabs(term);ok++;
   }} }while(std::next_permutation(perm.begin(),perm.end()));
  long double g=sum.s/(128.L*std::pow(std::numbers::pi_v<long double>,5));long double secs=std::chrono::duration<long double>(std::chrono::steady_clock::now()-start).count();
  std::cout<<"    {\"index\":"<<p<<",\"numerators\":[";for(int i=0;i<6;i++){if(i)std::cout<<",";std::cout<<xs[p][i];}
  std::cout<<"],\"denominator\":"<<den<<",\"raw_terms\":967680,\"accepted_terms\":"<<ok<<",\"zero_node_terms\":"<<zeros<<",\"repeated_node_pairs\":"<<repeats<<",\"status\":\""<<(zeros||repeats?"STRUCTURAL_FAILURE":"GENERIC_OK")<<"\",\"signed_sum_long_double\":\""<<std::setprecision(21)<<sum.s<<"\",\"absolute_sum_long_double\":\""<<abs_sum<<"\",\"cancellation_ratio\":\""<<(abs_sum/std::fabs(sum.s))<<"\",\"G_direct_long_double\":\""<<g<<"\",\"elapsed_seconds\":\""<<secs<<"\"}"<<(p+1==wanted?"\n":",\n");
 }
 std::cout<<"  ]\n}\n";
}
