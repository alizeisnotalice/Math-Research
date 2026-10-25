// Arbitrary-precision spatial kernel; independent round012 engine supplies moments. K3 determinant family.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "round012_graph_moments_gmp.hpp"

using I = round012::Integer;
constexpr int SHIFT=5, MASK=31, MAXM=27, MAXD=7, MAXW=11;
using U = uint64_t;
using X = std::array<int,5>;
using Poly = std::map<U,I>;
I add(const I&a,const I&b){return a+b;}
I mul(const I&a,const I&b){return a*b;}
std::string str(const I&x){return x.get_str();}
int edge[5][5], eu[10],ev[10];
std::unique_ptr<round012::Engine> engine;
I fact[MAXM+1];
I moment(U code){return engine->B(code);}
int parity_of(U code){int p=0;for(int e=0;e<10;e++)if((code>>(SHIFT*e))&1)p^=(1<<eu[e])|(1<<ev[e]);return p;}
int degree(U code){int n=0;for(int e=0;e<10;e++)n+=(code>>(SHIFT*e))&MASK;return n;}
int permutation_sign(const std::vector<int>&p){int inv=0;for(size_t i=0;i<p.size();i++)for(size_t j=i+1;j<p.size();j++)inv+=p[i]>p[j];return inv%2?-1:1;}
Poly minor_poly(const std::vector<int>&rows,const std::vector<int>&cols){
 Poly out;std::vector<int> p(rows.size());for(size_t i=0;i<p.size();i++)p[i]=i;
 do{U code=0;for(size_t i=0;i<p.size();i++){int a=rows[i],b=cols[p[i]];if(a!=b)code+=U(1)<<(SHIFT*edge[a][b]);}out[code]+=permutation_sign(p);}while(std::next_permutation(p.begin(),p.end()));
 for(auto it=out.begin();it!=out.end();)if(it->second==0)it=out.erase(it);else++it;return out;
}
struct BTerm{U code;I coefficient;X spatial;int parity;int n;};
using Coeff = std::array<I,MAXM+1>;
using SymPoly = std::map<X,Coeff>;
std::vector<SymPoly> d1(MAXD+1),d2(MAXD+1),d3(MAXD+1),d4(MAXD+1),bm(MAXW+1);
Poly determinant,determinant_squared,determinant_cubed,determinant_fourth;
static_assert(MAXD+20<=MAXM && MAXW+8<32,"graph code must not carry or exceed degree budget");
std::vector<BTerm> bterms;
std::array<std::vector<BTerm>,32> b_by_parity;

bool descending(const X&x){for(int i=1;i<5;i++)if(x[i]>x[i-1])return false;return true;}
void acc(SymPoly &out,const X&x,int n,I value){auto &row=out[x];row[n]=add(row[n],value);}
void process_W(int n,U code,const X&vertex,int parity,I multinomial){
 if(n<=MAXD&&parity==0){
  X x;for(int i=0;i<5;i++)x[i]=vertex[i]/2;
  if(descending(x)){
   for(const auto&t:determinant){int k=n+degree(t.first);I value=mul(mul(multinomial,t.second),moment(code+t.first));acc(d1[n],x,k,value);}
   for(const auto&t:determinant_squared){int k=n+degree(t.first);I value=mul(mul(multinomial,t.second),moment(code+t.first));acc(d2[n],x,k,value);}
   for(const auto&t:determinant_cubed){int k=n+degree(t.first);assert(k<=MAXM);I value=mul(mul(multinomial,t.second),moment(code+t.first));acc(d3[n],x,k,value);}
   for(const auto&t:determinant_fourth){int k=n+degree(t.first);assert(k<=MAXM);I value=mul(mul(multinomial,t.second),moment(code+t.first));acc(d4[n],x,k,value);}
  }
 }
 for(int i=1;i<5;i++)if(vertex[i]>vertex[i-1]+2)return;
 for(const auto&t:b_by_parity[parity]){
  X x;for(int i=0;i<5;i++){assert((vertex[i]+t.spatial[i])%2==0);x[i]=(vertex[i]+t.spatial[i])/2;}
  if(!descending(x))continue;
  I value=mul(mul(multinomial,t.coefficient),moment(code+t.code));acc(bm[n],x,n+t.n,value);
 }
}
void enumerate_W(int e,int remaining,int n,U code,X&vertex,int parity,I divisor){
 if(e==10){if(remaining==0)process_W(n,code,vertex,parity,fact[n]/divisor);return;}
 for(int k=0;k<=remaining;k++){
  vertex[eu[e]]+=k;vertex[ev[e]]+=k;
  enumerate_W(e+1,remaining-k,n,code+(U(k)<<(SHIFT*e)),vertex,
    parity^((k&1)?((1<<eu[e])|(1<<ev[e])):0),mul(divisor,fact[k]));
  vertex[eu[e]]-=k;vertex[ev[e]]-=k;
 }
}
void write_family(std::ofstream&out,const std::vector<SymPoly>&family){
 out<<"[";for(size_t n=0;n<family.size();n++){if(n)out<<",";out<<"{";bool first=true;
  for(const auto&term:family[n]){bool nz=false;for(const I&q:term.second)if(q!=0)nz=true;if(!nz)continue;if(!first)out<<",";first=false;
   out<<"\"";for(int i=0;i<5;i++){if(i)out<<",";out<<term.first[i];}out<<"\":[";
   for(int k=0;k<=MAXM;k++){if(k)out<<",";out<<"\""<<str(term.second[k])<<"\"";}out<<"]";
  }out<<"}";
 }out<<"]";
}
int main(int argc,char**argv){
 if(argc!=3){std::cerr<<"output path and cache path required\n";return 2;}auto began=std::chrono::steady_clock::now();
 int ei=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++){edge[i][j]=edge[j][i]=ei;eu[ei]=i;ev[ei]=j;ei++;}
 engine=std::make_unique<round012::Engine>(argv[2],"7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8",0.0);
 fact[0]=1;for(int n=1;n<=MAXM;n++)fact[n]=mul(fact[n-1],n);
 assert(moment(U(2))==3);assert(moment(U(4))*55==engine->denominator(4)*3);
 determinant=minor_poly({0,1,2,3,4},{0,1,2,3,4});
 for(const auto&a:determinant)for(const auto&b:determinant)determinant_squared[a.first+b.first]+=a.second*b.second;
 for(auto it=determinant_squared.begin();it!=determinant_squared.end();)if(it->second==0)it=determinant_squared.erase(it);else ++it;
 for(const auto&a:determinant_squared)for(const auto&b:determinant)determinant_cubed[a.first+b.first]+=a.second*b.second;
 for(auto it=determinant_cubed.begin();it!=determinant_cubed.end();)if(it->second==0)it=determinant_cubed.erase(it);else ++it;
 for(const auto&a:determinant_cubed)for(const auto&b:determinant)determinant_fourth[a.first+b.first]+=a.second*b.second;
 for(auto it=determinant_fourth.begin();it!=determinant_fourth.end();)if(it->second==0)it=determinant_fourth.erase(it);else ++it;
 for(const auto&t:determinant_fourth){assert(degree(t.first)<=20);for(int e=0;e<10;e++)assert(((t.first>>(SHIFT*e))&MASK)<=8);}
 for(int i=0;i<5;i++)for(int j=i;j<5;j++){
  std::vector<int> rows,cols;for(int k=0;k<5;k++){if(k!=j)rows.push_back(k);if(k!=i)cols.push_back(k);}
  auto poly=minor_poly(rows,cols);int sign=((i+j)%2?-1:1)*(i==j?1:2);X spatial{};spatial[i]++;spatial[j]++;
  for(const auto&t:poly){int parity=parity_of(t.first);assert(parity==(i==j?0:((1<<i)|(1<<j))));BTerm term{t.first,t.second*sign,spatial,parity,degree(t.first)};bterms.push_back(term);b_by_parity[parity].push_back(term);}
 }
 std::cout<<"det terms "<<determinant.size()<<", det squared terms "<<determinant_squared.size()<<", det cubed terms "<<determinant_cubed.size()<<", det fourth terms "<<determinant_fourth.size()<<", cofactor terms "<<bterms.size()<<std::endl;
 for(int n=0;n<=MAXW;n++){X vertex{};enumerate_W(0,n,n,0,vertex,0,1);double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  engine->sync();std::cout<<"completed W power "<<n<<", canonical moments "<<engine->cached_moments()<<", seconds "<<seconds<<std::endl;
 std::string temp=std::string(argv[1])+".partial";std::ofstream out(temp);out<<"{\"method\":\"GMP arbitrary-precision round012 Eulerian graph backend\",\"backend_source_sha256\":\"7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8\",\"denominators\":[";
 for(int n=0;n<=MAXM;n++){if(n)out<<",";out<<"\""<<str(engine->denominator(n))<<"\"";}out<<"],\"d1\":";write_family(out,d1);out<<",\"d2\":";write_family(out,d2);out<<",\"d3\":";write_family(out,d3);out<<",\"d4\":";write_family(out,d4);out<<",\"b\":";write_family(out,bm);
 out<<",\"completed_W_power\":"<<n<<",\"seconds\":"<<seconds<<",\"canonical_moment_count\":"<<engine->cached_moments()<<",\"max_B_bits\":"<<engine->stats().max_B_bits<<"}\n";
 out.close();if(std::rename(temp.c_str(),argv[1]))throw std::runtime_error("atomic checkpoint rename failed");
 }
 std::cout<<"wrote "<<argv[1]<<std::endl;
}
