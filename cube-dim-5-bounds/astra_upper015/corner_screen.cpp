// New beta=1/2 moments, necessary corner only. No old numerical cache is read.
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "graph_moments_beta_half.hpp"
using I=astra015::Integer;using U=astra015::Code;using Poly=std::map<U,I>;
constexpr const char* HEADER_SHA="bfdf502cac067d2a1d9465c96ac15076c2e897cee89f8121bfed40977fc3fbf0";
int edge[5][5];
U add_code(U a,U b){auto x=astra015::decode(a),y=astra015::decode(b);
 for(unsigned e=0;e<10;++e)x[e]+=y[e];
 if(astra015::degree(x)>7)throw std::invalid_argument("screen degree exceeds7");return astra015::encode(x);}
void trim(Poly&p){for(auto it=p.begin();it!=p.end();)if(it->second==0)it=p.erase(it);else++it;}
Poly product(const Poly&a,const Poly&b){Poly c;for(const auto&x:a)for(const auto&y:b)c[add_code(x.first,y.first)]+=x.second*y.second;trim(c);return c;}
int sign(const std::vector<int>&p){int n=0;for(size_t i=0;i<p.size();++i)for(size_t j=i+1;j<p.size();++j)n+=p[i]>p[j];return n%2?-1:1;}
Poly minor_poly(const std::vector<int>&rows,const std::vector<int>&cols){Poly out;std::vector<int>p(rows.size());for(size_t i=0;i<p.size();++i)p[i]=int(i);
 do{astra015::Counts k{};for(size_t i=0;i<p.size();++i)if(rows[i]!=cols[p[i]])++k[edge[rows[i]][cols[p[i]]]];
 out[astra015::encode(k)]+=sign(p);}while(std::next_permutation(p.begin(),p.end()));trim(out);return out;}
mpq_class expectation(const Poly&p,astra015::Engine&e){mpq_class ans=0;for(const auto&t:p){auto k=astra015::decode(t.first);mpq_class z(t.second*e.B(t.first),e.denominator(astra015::degree(k)));z.canonicalize();ans+=z;}ans.canonicalize();return ans;}
int main(int argc,char**argv){
 if(argc!=3){std::cerr<<"output.json new_beta_cache.tsv\n";return 2;}
 int k=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)edge[i][j]=edge[j][i]=k++;
 astra015::Engine e(argv[2],HEADER_SHA,120.0);
 try{
  Poly one{{0,I(1)}},W;for(int i=0;i<10;++i)W[U(1)<<(5*i)]=1;
  Poly d=minor_poly({0,1,2,3,4},{0,1,2,3,4}),b;
  for(int i=0;i<5;i++)for(int j=0;j<5;j++){std::vector<int>rows,cols;for(int a=0;a<5;a++){if(a!=j)rows.push_back(a);if(a!=i)cols.push_back(a);}
   auto p=minor_poly(rows,cols);for(const auto&t:p)b[t.first]+=((i+j)%2?-1:1)*t.second;}
  trim(b);
  std::map<std::string,std::vector<mpq_class>> result;
  for(const auto&item:std::vector<std::pair<std::string,Poly>>{{"one",one},{"d",d},{"b",b}}){
   Poly f=item.second;int maxn=item.first=="b"?3:2;
   for(int n=0;n<=maxn;n++){result[item.first].push_back(expectation(f,e));if(n<maxn)f=product(f,W);}}
  std::vector<std::pair<astra015::Counts,mpq_class>> tests;
  for(unsigned n=0;n<=6;n++){astra015::Counts c{};c[0]=n;auto value=e.expectation(c);
   if(value!=(n%2?mpq_class(0):mpq_class(1,n+1)))throw std::runtime_error("uniform single-edge check failed");tests.push_back({c,value});}
  for(const auto&c:std::vector<astra015::Counts>{{{2,2,0,0,2,0,0,0,0,0}},{{1,0,1,0,1,0,0,1,0,0}},{{2,2,0,0,0,0,0,0,0,0}},{{2,0,0,0,0,0,0,2,0,0}}})tests.push_back({c,e.expectation(c)});
  e.sync();std::ofstream o(std::string(argv[1])+".partial");o<<"{\"beta\":\"1/2\",\"status\":\"necessary_corner_exact_moments\",\"backend_sha256\":\""<<HEADER_SHA<<"\",\"moments\":{";
  bool first=true;for(const auto&item:result){if(!first)o<<",";first=false;o<<"\""<<item.first<<"\":[";for(size_t n=0;n<item.second.size();++n){if(n)o<<",";o<<"\""<<item.second[n].get_str()<<"\"";}o<<"]";}o<<"},\"independent_check_targets\":[";
  first=true;for(const auto&t:tests){if(!first)o<<",";first=false;o<<"{\"counts\":[";for(unsigned i=0;i<10;i++){if(i)o<<",";o<<t.first[i];}o<<"],\"moment\":\""<<t.second.get_str()<<"\"}";}o<<"],\"cache_count\":"<<e.cached_moments()<<",\"loaded\":"<<e.stats().loaded<<",\"new\":"<<e.stats().computed<<",\"seconds\":"<<e.elapsed()<<",\"max_B_bits\":"<<e.stats().max_B_bits<<"}\n";o.close();
  if(std::rename((std::string(argv[1])+".partial").c_str(),argv[1]))throw std::runtime_error("atomic save failed");
  std::cout<<"saved "<<argv[1]<<" with "<<e.cached_moments()<<" beta-half moments in "<<e.elapsed()<<" sec\n";
 }catch(const astra015::CheckpointStop&x){e.sync();std::cerr<<x.what()<<"; graph journal saved, restart screen from empty sums\n";return 3;}
}
