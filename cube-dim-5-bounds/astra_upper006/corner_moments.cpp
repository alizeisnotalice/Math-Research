// One necessary-corner pilot. All moment arithmetic is arbitrary precision.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "round012_graph_moments_gmp.hpp"
using I=round012::Integer;
using U=round012::Code;
using Poly=std::map<U,I>;
constexpr int MAXN=7,MAXD=27;
constexpr const char* HEADER_SHA="7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8";
int edge[5][5];
std::unique_ptr<round012::Engine> engine;
I fact[MAXN+1];
std::array<Poly,5> determinant;
std::array<std::array<mpq_class,MAXN+1>,5> result;

U add_code(U a,U b){
 auto x=round012::decode(a),y=round012::decode(b);
 for(unsigned i=0;i<10;++i){x[i]+=y[i];assert(x[i]<=31);}
 assert(round012::degree(x)<=MAXD);return round012::encode(x);
}
int permutation_sign(const std::array<int,5>&p){int inv=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)inv+=p[i]>p[j];return inv%2?-1:1;}
void erase_zeros(Poly&p){for(auto it=p.begin();it!=p.end();)if(it->second==0)it=p.erase(it);else ++it;}
void make_determinants(){
 std::array<int,5> p{{0,1,2,3,4}};
 do{round012::Counts k{};for(int i=0;i<5;i++)if(i!=p[i])++k[edge[i][p[i]]];determinant[1][round012::encode(k)]+=permutation_sign(p);}while(std::next_permutation(p.begin(),p.end()));
 erase_zeros(determinant[1]);
 for(int k=2;k<=4;k++){
  for(const auto&a:determinant[k-1])for(const auto&b:determinant[1])determinant[k][add_code(a.first,b.first)]+=a.second*b.second;
  erase_zeros(determinant[k]);
 }
 for(int k=1;k<=4;k++){
  Poly compressed;I l1=0;
  for(const auto&t:determinant[k]){assert(round012::parity(round012::decode(t.first))==0);compressed[engine->canonical(t.first)]+=t.second;l1+=abs(t.second);}
  std::cout<<"det power "<<k<<", raw terms "<<determinant[k].size()<<", orbit terms "<<compressed.size()<<", L1 "<<l1<<std::endl;
  determinant[k]=std::move(compressed);
 }
}
void enumerate_W(int edge_index,int remaining,int n,round012::Counts&k,I divisor,Poly&out){
 if(edge_index==10){if(remaining==0&&round012::parity(k)==0)out[round012::encode(k)]=fact[n]/divisor;return;}
 for(int j=0;j<=remaining;j++){k[edge_index]=j;enumerate_W(edge_index+1,remaining-j,n,k,I(divisor*fact[j]),out);}k[edge_index]=0;
}
void save(const std::string&path,int completed,double seconds){
 std::ofstream o(path+".partial");o<<"{\"status\":\"necessary_corner_only\",\"backend_source_sha256\":\""<<HEADER_SHA<<"\",\"completed_W_power\":"<<completed<<",\"moments\":{";
 for(int k=1;k<=4;k++){if(k>1)o<<",";o<<"\""<<k<<"\":[";for(int n=0;n<=completed;n++){if(n)o<<",";o<<"\""<<result[k][n].get_str()<<"\"";}o<<"]";}o<<"},\"canonical_moments\":"<<engine->cached_moments()<<",\"max_B_bits\":"<<engine->stats().max_B_bits<<",\"seconds\":"<<seconds<<"}\n";o.close();if(std::rename((path+".partial").c_str(),path.c_str()))throw std::runtime_error("atomic checkpoint rename failed");
}
int main(int argc,char**argv){
 if(argc!=3){std::cerr<<"usage: corner_moments output.json cache.tsv\n";return 2;}
 auto began=std::chrono::steady_clock::now();
 int e=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++){edge[i][j]=edge[j][i]=e++;}
 engine=std::make_unique<round012::Engine>(argv[2],HEADER_SHA,0.0);
 fact[0]=1;for(int n=1;n<=MAXN;n++)fact[n]=fact[n-1]*n;
 make_determinants();
 for(int n=0;n<=MAXN;n++){
  Poly w;round012::Counts counts{};enumerate_W(0,n,n,counts,1,w);
  std::cout<<"W "<<n<<" Eulerian terms "<<w.size()<<std::endl;
  for(int k=1;k<=4;k++){
   std::array<I,MAXD+1> sums{};
   for(const auto&d:determinant[k])for(const auto&t:w){U code=add_code(d.first,t.first);unsigned degree=round012::degree(round012::decode(code));sums[degree]+=d.second*t.second*engine->B(code);}
   mpq_class q=0;for(int degree=0;degree<=MAXD;degree++){mpq_class term(sums[degree],engine->denominator(degree));term.canonicalize();q+=term;}
   q.canonicalize();result[k][n]=q;
   std::cout<<" completed d power "<<k<<", canonical moments "<<engine->cached_moments()<<std::endl;
  }
  engine->sync();double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();save(argv[1],n,elapsed);std::cout<<"saved W "<<n<<", seconds "<<elapsed<<std::endl;
 }
}
