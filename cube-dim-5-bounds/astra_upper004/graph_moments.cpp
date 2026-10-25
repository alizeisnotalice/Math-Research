// Exact signed 128-bit graph moments, with overflow checks on every large operation.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using I = __int128_t;
using U = uint64_t;
using X = std::array<int,5>;
using Poly = std::map<U,int64_t>;
I add(I a,I b){I c;if(__builtin_add_overflow(a,b,&c))throw std::runtime_error("128-bit add overflow");return c;}
I mul(I a,I b){I c;if(__builtin_mul_overflow(a,b,&c))throw std::runtime_error("128-bit product overflow");return c;}
std::string str(I x){if(!x)return "0";bool neg=x<0;if(neg)x=-x;std::string s;while(x){s.push_back('0'+int(x%10));x/=10;}if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s;}
int edge[5][5], eu[10],ev[10];
std::vector<std::array<int,10>> perm_edges;
std::unordered_map<U,U> canonical_cache;
std::unordered_map<U,I> moment_cache;
I denom[18], cumul[18];
int64_t choose[18][18],fact[18];

U canonical(U code){
 auto found=canonical_cache.find(code);if(found!=canonical_cache.end())return found->second;
 std::array<int,10> k{};for(int e=0;e<10;e++)k[e]=(code>>(4*e))&15;
 U best=code;
 for(const auto &p:perm_edges){U q=0;for(int e=0;e<10;e++)q|=U(k[e])<<(4*p[e]);best=std::min(best,q);}
 canonical_cache.emplace(code,best);return best;
}
I moment(U code);
void subset_sum(int e,const std::array<int,10>&k,U sub,int size,int parity,int64_t binom,U whole,I&total){
 if(e==10){
  if(size&&parity==0){I term=mul(mul(cumul[size],moment(whole-sub)),I(binom));total=add(total,term);}return;
 }
 for(int j=0;j<=k[e];j++)subset_sum(e+1,k,sub+(U(j)<<(4*e)),size+j,
  parity^((j&1)?((1<<eu[e])|(1<<ev[e])):0),binom*choose[k[e]][j],whole,total);
}
I moment(U code){
 code=canonical(code);auto found=moment_cache.find(code);if(found!=moment_cache.end())return found->second;
 std::array<int,10> k{};int n=0,parity=0;for(int e=0;e<10;e++){k[e]=(code>>(4*e))&15;n+=k[e];if(k[e]&1)parity^=(1<<eu[e])|(1<<ev[e]);}
 assert(parity==0&&n<=17);I total=0;subset_sum(0,k,0,0,0,1,code,total);assert(total>=0&&total<=denom[n]);moment_cache.emplace(code,total);return total;
}
int parity_of(U code){int p=0;for(int e=0;e<10;e++)if((code>>(4*e))&1)p^=(1<<eu[e])|(1<<ev[e]);return p;}
int degree(U code){int n=0;for(int e=0;e<10;e++)n+=(code>>(4*e))&15;return n;}
int permutation_sign(const std::vector<int>&p){int inv=0;for(size_t i=0;i<p.size();i++)for(size_t j=i+1;j<p.size();j++)inv+=p[i]>p[j];return inv%2?-1:1;}
Poly minor_poly(const std::vector<int>&rows,const std::vector<int>&cols){
 Poly out;std::vector<int> p(rows.size());for(size_t i=0;i<p.size();i++)p[i]=i;
 do{U code=0;for(size_t i=0;i<p.size();i++){int a=rows[i],b=cols[p[i]];if(a!=b)code+=U(1)<<(4*edge[a][b]);}out[code]+=permutation_sign(p);}while(std::next_permutation(p.begin(),p.end()));
 for(auto it=out.begin();it!=out.end();)if(it->second==0)it=out.erase(it);else++it;return out;
}
struct BTerm{U code;int64_t coefficient;X spatial;int parity;int n;};
using Coeff = std::array<I,18>;
using SymPoly = std::map<X,Coeff>;
std::vector<SymPoly> d1(8),d2(8),bm(12);
Poly determinant,determinant_squared;
std::vector<BTerm> bterms;

bool descending(const X&x){for(int i=1;i<5;i++)if(x[i]>x[i-1])return false;return true;}
void acc(SymPoly &out,const X&x,int n,I value){auto &row=out[x];row[n]=add(row[n],value);}
void process_W(int n,U code,const X&vertex,int parity,int64_t multinomial){
 if(n<=7&&parity==0){
  X x;for(int i=0;i<5;i++)x[i]=vertex[i]/2;
  if(descending(x)){
   for(const auto&t:determinant){int k=n+degree(t.first);I value=mul(I(multinomial*t.second),moment(code+t.first));acc(d1[n],x,k,value);}
   for(const auto&t:determinant_squared){int k=n+degree(t.first);I value=mul(I(multinomial*t.second),moment(code+t.first));acc(d2[n],x,k,value);}
  }
 }
 for(const auto&t:bterms){
  if(parity!=t.parity)continue;
  X x;for(int i=0;i<5;i++){assert((vertex[i]+t.spatial[i])%2==0);x[i]=(vertex[i]+t.spatial[i])/2;}
  if(!descending(x))continue;
  I value=mul(I(multinomial*t.coefficient),moment(code+t.code));acc(bm[n],x,n+t.n,value);
 }
}
void enumerate_W(int e,int remaining,int n,U code,X&vertex,int parity,int64_t divisor){
 if(e==10){if(remaining==0)process_W(n,code,vertex,parity,fact[n]/divisor);return;}
 for(int k=0;k<=remaining;k++){
  vertex[eu[e]]+=k;vertex[ev[e]]+=k;
  enumerate_W(e+1,remaining-k,n,code+(U(k)<<(4*e)),vertex,
    parity^((k&1)?((1<<eu[e])|(1<<ev[e])):0),divisor*fact[k]);
  vertex[eu[e]]-=k;vertex[ev[e]]-=k;
 }
}
void write_family(std::ofstream&out,const std::vector<SymPoly>&family){
 out<<"[";for(size_t n=0;n<family.size();n++){if(n)out<<",";out<<"{";bool first=true;
  for(const auto&term:family[n]){bool nz=false;for(I q:term.second)if(q)nz=true;if(!nz)continue;if(!first)out<<",";first=false;
   out<<"\"";for(int i=0;i<5;i++){if(i)out<<",";out<<term.first[i];}out<<"\":[";
   for(int k=0;k<18;k++){if(k)out<<",";out<<"\""<<str(term.second[k])<<"\"";}out<<"]";
  }out<<"}";
 }out<<"]";
}
int main(int argc,char**argv){
 if(argc!=2){std::cerr<<"output path required\n";return 2;}auto began=std::chrono::steady_clock::now();
 int ei=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++){edge[i][j]=edge[j][i]=ei;eu[ei]=i;ev[ei]=j;ei++;}
 std::array<int,5> p={0,1,2,3,4};do{std::array<int,10> q{};for(int e=0;e<10;e++)q[e]=edge[p[eu[e]]][p[ev[e]]];perm_edges.push_back(q);}while(std::next_permutation(p.begin(),p.end()));
 denom[0]=1;cumul[0]=0;for(int n=1;n<=17;n++){denom[n]=mul(denom[n-1],4*n-1);cumul[n]=(n==1?3:mul(cumul[n-1],4*n-7));}
 fact[0]=1;for(int n=1;n<=17;n++)fact[n]=fact[n-1]*n;
 for(int n=0;n<=17;n++){choose[n][0]=choose[n][n]=1;for(int k=1;k<n;k++)choose[n][k]=choose[n-1][k-1]+choose[n-1][k];}
 moment_cache[0]=1;canonical_cache[0]=0;
 assert(moment(U(2))==3);assert(moment(U(4))*55==denom[4]*3);
 determinant=minor_poly({0,1,2,3,4},{0,1,2,3,4});
 for(const auto&a:determinant)for(const auto&b:determinant)determinant_squared[a.first+b.first]+=a.second*b.second;
 for(auto it=determinant_squared.begin();it!=determinant_squared.end();)if(!it->second)it=determinant_squared.erase(it);else ++it;
 for(int i=0;i<5;i++)for(int j=i;j<5;j++){
  std::vector<int> rows,cols;for(int k=0;k<5;k++){if(k!=j)rows.push_back(k);if(k!=i)cols.push_back(k);}
  auto poly=minor_poly(rows,cols);int sign=((i+j)%2?-1:1)*(i==j?1:2);X spatial{};spatial[i]++;spatial[j]++;
  for(const auto&t:poly){int parity=parity_of(t.first);assert(parity==(i==j?0:((1<<i)|(1<<j))));bterms.push_back({t.first,t.second*sign,spatial,parity,degree(t.first)});}
 }
 std::cout<<"det terms "<<determinant.size()<<", det squared terms "<<determinant_squared.size()<<", cofactor terms "<<bterms.size()<<std::endl;
 for(int n=0;n<=11;n++){X vertex{};enumerate_W(0,n,n,0,vertex,0,1);double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
  std::cout<<"completed W power "<<n<<", canonical moments "<<moment_cache.size()<<", raw canonical cache "<<canonical_cache.size()<<", seconds "<<seconds<<std::endl;
 }
 std::ofstream out(argv[1]);out<<"{\"method\":\"exact 128-bit Eulerian graph recurrence with checked arithmetic\",\"denominators\":[";
 for(int n=0;n<=17;n++){if(n)out<<",";out<<"\""<<str(denom[n])<<"\"";}out<<"],\"d1\":";write_family(out,d1);out<<",\"d2\":";write_family(out,d2);out<<",\"b\":";write_family(out,bm);
 out<<",\"canonical_moment_count\":"<<moment_cache.size()<<",\"canonical_cache_count\":"<<canonical_cache.size()<<"}\n";
 std::cout<<"wrote "<<argv[1]<<std::endl;
}
