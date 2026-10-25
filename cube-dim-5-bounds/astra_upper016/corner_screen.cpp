// Two fixed m values, new beta-half mixed moments. One invariant factor only is orbit compressed.
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "../astra_upper015/graph_moments_beta_half.hpp"
using I=astra015::Integer;using U=astra015::Code;using Poly=std::map<U,I>;
constexpr const char* HEADER_SHA="bfdf502cac067d2a1d9465c96ac15076c2e897cee89f8121bfed40977fc3fbf0";
int edge[5][5];I fact[10];
U add_code(U a,U b){auto x=astra015::decode(a),y=astra015::decode(b);for(unsigned e=0;e<10;++e)x[e]+=y[e];
 if(astra015::degree(x)>13)throw std::invalid_argument("screen degree exceeds13");return astra015::encode(x);}
void trim(Poly&p){for(auto it=p.begin();it!=p.end();)if(it->second==0)it=p.erase(it);else++it;}
int sign(const std::vector<int>&p){int n=0;for(size_t i=0;i<p.size();++i)for(size_t j=i+1;j<p.size();++j)n+=p[i]>p[j];return n%2?-1:1;}
Poly minor_poly(const std::vector<int>&rows,const std::vector<int>&cols){Poly out;std::vector<int>p(rows.size());for(size_t i=0;i<p.size();++i)p[i]=int(i);
 do{astra015::Counts k{};for(size_t i=0;i<p.size();++i)if(rows[i]!=cols[p[i]])++k[edge[rows[i]][cols[p[i]]]];
 out[astra015::encode(k)]+=sign(p);}while(std::next_permutation(p.begin(),p.end()));trim(out);return out;}
Poly compress(const Poly&p,astra015::Engine&e){Poly out;for(const auto&t:p)out[e.canonical(t.first)]+=t.second;trim(out);return out;}
void enumerate_W(unsigned edge_index,unsigned remaining,unsigned n,astra015::Counts&k,I divisor,std::array<Poly,32>&out){
 if(edge_index==10){if(!remaining)out[astra015::parity(k)][astra015::encode(k)]=fact[n]/divisor;return;}
 for(unsigned j=0;j<=remaining;++j){k[edge_index]=j;enumerate_W(edge_index+1,remaining-j,n,k,I(divisor*fact[j]),out);}k[edge_index]=0;}
mpq_class expectation(const Poly&p,const std::array<Poly,32>&w,astra015::Engine&e){std::array<I,14>sums{};
 for(const auto&t:p)for(const auto&v:w[astra015::parity(astra015::decode(t.first))]){U k=add_code(t.first,v.first);sums[astra015::degree(astra015::decode(k))]+=t.second*v.second*e.B(k);}
 mpq_class ans=0;for(unsigned n=0;n<=13;++n){mpq_class z(sums[n],e.denominator(n));z.canonicalize();ans+=z;}ans.canonicalize();return ans;}
void save(const std::string&path,const std::map<std::string,std::vector<mpq_class>>&result,unsigned completed,astra015::Engine&e){
 std::ofstream o(path+".partial");o<<"{\"beta\":\"1/2\",\"status\":\"two_fixed_ridge_power_necessary_corner_moments\",\"backend_sha256\":\""<<HEADER_SHA<<"\",\"completed_W_power\":"<<completed<<",\"moments\":{";
 bool first=true;for(const auto&item:result){if(!first)o<<",";first=false;o<<"\""<<item.first<<"\":[";for(size_t n=0;n<item.second.size();++n){if(n)o<<",";o<<"\""<<item.second[n].get_str()<<"\"";}o<<"]";}o<<"},\"cache_count\":"<<e.cached_moments()<<",\"loaded\":"<<e.stats().loaded<<",\"new\":"<<e.stats().computed<<",\"seconds\":"<<e.elapsed()<<",\"max_B_bits\":"<<e.stats().max_B_bits<<"}\n";o.close();if(std::rename((path+".partial").c_str(),path.c_str()))throw std::runtime_error("atomic save failed");}
int main(int argc,char**argv){
 if(argc!=3){std::cerr<<"output.json beta_half_cache.tsv\n";return 2;}
 int k=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)edge[i][j]=edge[j][i]=k++;
 fact[0]=1;for(int n=1;n<=9;n++)fact[n]=fact[n-1]*n;
 astra015::Engine e(argv[2],HEADER_SHA,120.0);
 try{
  Poly one{{0,I(1)}},d=minor_poly({0,1,2,3,4},{0,1,2,3,4}),b;
  for(int i=0;i<5;i++)for(int j=0;j<5;j++){std::vector<int>rows,cols;for(int a=0;a<5;a++){if(a!=j)rows.push_back(a);if(a!=i)cols.push_back(a);}
   auto p=minor_poly(rows,cols);for(const auto&t:p)b[t.first]+=((i+j)%2?-1:1)*t.second;}
  trim(b);d=compress(d,e);b=compress(b,e);
  std::map<std::string,std::vector<mpq_class>> result;
  for(unsigned n=0;n<=9;++n){std::array<Poly,32>w;astra015::Counts counts{};enumerate_W(0,n,n,counts,I(1),w);
   if(n<=6){result["one"].push_back(expectation(one,w,e));result["d"].push_back(expectation(d,w,e));}
   result["b"].push_back(expectation(b,w,e));e.sync();save(argv[1],result,n,e);
   std::cout<<"completed W "<<n<<", moments "<<e.cached_moments()<<", seconds "<<e.elapsed()<<std::endl;
   if(e.elapsed()>120)throw astra015::CheckpointStop("complete-W checkpoint reached");}
 }catch(const astra015::CheckpointStop&x){e.sync();std::cerr<<x.what()<<"; graph journal and prior complete W powers saved; restart sums empty\n";return 3;}
}
