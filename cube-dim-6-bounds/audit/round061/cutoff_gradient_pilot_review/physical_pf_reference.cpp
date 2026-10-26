#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "arb.h"
#include "fmpz.h"

constexpr slong PREC=768;
constexpr uint64_t FULL=967680;
struct B { arb_t v; B(){arb_init(v);} ~B(){arb_clear(v);}
 B(const B&o){arb_init(v);arb_set(v,o.v);} B& operator=(const B&o){arb_set(v,o.v);return *this;} };
struct Q { int64_t n,d; bool operator==(const Q&o)const{return n==o.n&&d==o.d;} };
struct Hash { size_t operator()(const Q&q)const{return std::hash<int64_t>{}(q.n)^(std::hash<int64_t>{}(q.d)+0x9e3779b97f4a7c15ULL+(uint64_t(q.n)<<1));} };
Q rat(int64_t n,int64_t d){if(d<=0||n==INT64_MIN)throw std::runtime_error("rational range");auto g=std::gcd(n,d);return{n/g,d/g};}
int64_t checked(__int128 x){if(x<=INT64_MIN||x>INT64_MAX)throw std::runtime_error("int64 range");return int64_t(x);}
void setq(B&b,Q q){arb_set_si(b.v,q.n);arb_div_ui(b.v,b.v,q.d,PREC);}
void output(const B&b){if(!arb_is_finite(b.v))throw std::runtime_error("nonfinite output");fmpz_t lo,hi,ex;fmpz_init(lo);fmpz_init(hi);fmpz_init(ex);arb_get_interval_fmpz_2exp(lo,hi,ex,b.v);char*l=fmpz_get_str(nullptr,10,lo),*h=fmpz_get_str(nullptr,10,hi),*e=fmpz_get_str(nullptr,10,ex);std::cout<<"{\"lower_mantissa\":\""<<l<<"\",\"upper_mantissa\":\""<<h<<"\",\"common_binary_exponent\":\""<<e<<"\"}";flint_free(l);flint_free(h);flint_free(e);fmpz_clear(lo);fmpz_clear(hi);fmpz_clear(ex);}
void output6(const std::array<B,6>&v){std::cout<<"[";for(int k=0;k<6;k++){if(k)std::cout<<",";output(v[k]);}std::cout<<"]";}

struct Node { B b; std::array<B,8>phi;
 Node(Q q,const B&a,const B&T){setq(b,q);B x,y;arb_div(x.v,b.v,a.v,PREC);arb_atan(x.v,x.v,PREC);arb_div(y.v,b.v,T.v,PREC);arb_atan(y.v,y.v,PREC);arb_sub(phi[1].v,x.v,y.v,PREC);
  std::array<std::array<B,7>,2>re,im;
  for(int end=0;end<2;end++){const B&u=end?T:a;B den,rr,ii,nr,ni;arb_mul(den.v,u.v,u.v,PREC);arb_addmul(den.v,b.v,b.v,PREC);arb_div(rr.v,u.v,den.v,PREC);arb_div(ii.v,b.v,den.v,PREC);arb_one(re[end][0].v);
   for(int n=1;n<=6;n++){arb_mul(nr.v,re[end][n-1].v,rr.v,PREC);arb_submul(nr.v,im[end][n-1].v,ii.v,PREC);arb_mul(ni.v,re[end][n-1].v,ii.v,PREC);arb_addmul(ni.v,im[end][n-1].v,rr.v,PREC);re[end][n]=nr;im[end][n]=ni;}}
  for(int r=2;r<=7;r++){int rot=(r-6+8)%4;const auto&part=(rot%2)?im:re;arb_sub(phi[r].v,part[0][r-1].v,part[1][r-1].v,PREC);if(rot==1||rot==2)arb_neg(phi[r].v,phi[r].v);arb_div_ui(phi[r].v,phi[r].v,r-1,PREC);}
 }
};
class PhysicalPF {
 B a,T;std::unordered_map<Q,std::unique_ptr<Node>,Hash>cache;
 public: PhysicalPF(){setq(a,rat(1,250));setq(T,rat(1,2));cache.reserve(50000);}
 Node& node(Q q){auto it=cache.find(q);if(it==cache.end())it=cache.emplace(q,std::make_unique<Node>(q,a,T)).first;return *it->second;}
 size_t node_count()const{return cache.size();}
 void evaluate(const std::array<Q,6>&base,const std::array<std::array<B,6>,6>&direction,std::array<B,6>&out){
  std::array<Q,6>key;std::array<Node*,6>nd;std::array<int,6>mult{};std::array<std::array<B,6>,6>delta;int groups=0;
  for(int j=0;j<6;j++){int g=0;for(;g<groups&&!(key[g]==base[j]);g++);if(g==groups){key[g]=base[j];nd[g]=&node(base[j]);groups++;}mult[g]++;for(int k=0;k<6;k++)arb_add(delta[g][k].v,delta[g][k].v,direction[j][k].v,PREC);}
  std::array<B,6>coefficient;
  for(int g=0;g<groups;g++){int m=mult[g];std::array<B,7>q;std::array<B,6>inv;arb_one(q[0].v);
   for(int h=0;h<groups;h++)if(h!=g){B diff;arb_sub(diff.v,nd[g]->b.v,nd[h]->b.v,PREC);if(arb_contains_zero(diff.v))throw std::runtime_error("distinct exact poles not separated");arb_inv(inv[h].v,diff.v,PREC);
    std::array<B,7>factor,next;arb_pow_ui(factor[0].v,inv[h].v,mult[h],PREC);
    for(int l=1;l<=m;l++){arb_mul(factor[l].v,factor[l-1].v,inv[h].v,PREC);arb_mul_ui(factor[l].v,factor[l].v,mult[h]+l-1,PREC);arb_div_ui(factor[l].v,factor[l].v,l,PREC);arb_neg(factor[l].v,factor[l].v);}
    for(int l=0;l<=m;l++)for(int j=0;j<=l;j++)arb_addmul(next[l].v,q[j].v,factor[l-j].v,PREC);q=next;
   }
   B u;for(int r=1;r<=m+1;r++)arb_addmul(u.v,q[m+1-r].v,nd[g]->phi[r].v,PREC);arb_add(coefficient[g].v,coefficient[g].v,u.v,PREC);
   std::array<B,6>bl;for(int l=0;l<m;l++)for(int r=1;r<=m-l;r++)arb_addmul(bl[l].v,q[m-r-l].v,nd[g]->phi[r].v,PREC);
   for(int h=0;h<groups;h++)if(h!=g){B power,v;power=inv[h];for(int l=0;l<m;l++){arb_addmul(v.v,bl[l].v,power.v,PREC);arb_mul(power.v,power.v,inv[h].v,PREC);arb_neg(power.v,power.v);}arb_add(coefficient[h].v,coefficient[h].v,v.v,PREC);}
  }
  for(int k=0;k<6;k++){arb_zero(out[k].v);for(int g=0;g<groups;g++)arb_addmul(out[k].v,delta[g][k].v,coefficient[g].v,PREC);if(!arb_is_finite(out[k].v))throw std::runtime_error("nonfinite physical gradient");}
 }
};
struct Cone {int64_t det;std::array<int64_t,6>q;std::array<std::array<int64_t,6>,6>ray;std::array<std::array<B,6>,6>d;B alpha;};
std::array<Cone,21> read_fan(const char*path){std::ifstream in(path);int count;in>>count;if(!in||count!=21)throw std::runtime_error("fan count");std::array<Cone,21>fan;
 for(auto&c:fan){in>>c.det;if(!in||c.det<=0||c.det>1000000000)throw std::runtime_error("fan det");arb_set_si(c.alpha.v,c.det);
  for(int j=0;j<6;j++){in>>c.q[j];if(!in||c.q[j]<=0||c.q[j]>1000000)throw std::runtime_error("fan q");arb_div_ui(c.alpha.v,c.alpha.v,c.q[j],PREC);for(int k=0;k<6;k++){in>>c.ray[j][k];if(!in||c.ray[j][k]<-1000000||c.ray[j][k]>1000000)throw std::runtime_error("fan ray");setq(c.d[j][k],rat(c.ray[j][k],c.q[j]));}}}
 std::string extra;if(in>>extra||in.bad())throw std::runtime_error("fan termination");return fan;}

void scalars(){PhysicalPF engine;const std::array<std::array<Q,6>,5>sets={{{{{-3,8},{-1,9},{0,1},{2,7},{3,5},{5,6}}},{{{0,1},{0,1},{1,7},{1,7},{-2,5},{3,4}}},{{{0,1},{0,1},{0,1},{0,1},{0,1},{0,1}}},{{{2,11},{2,11},{2,11},{2,11},{2,11},{2,11}}},{{{-2,13},{-2,13},{-2,13},{0,1},{0,1},{5,17}}}}};
 std::array<std::array<B,6>,6>d;for(int j=0;j<6;j++)for(int k=0;k<6;k++)setq(d[j][k],rat(((j+k)%2)?-1:1,j+k+1));
 std::cout<<"{\"status\":\"FIVE_NEW_PHYSICAL_PF_SCALARS\",\"precision_bits\":768,\"rows\":[";
 for(int c=0;c<5;c++){std::array<B,6>out;engine.evaluate(sets[c],d,out);if(c)std::cout<<",";std::cout<<"{\"configuration\":"<<c+3<<",\"gradient\":";output6(out);std::cout<<"}";}std::cout<<"]}\n";
}

void real_case(const std::string&mode,const std::string&case_id,const char*fan_path,int block,int blocks){
 auto fan=read_fan(fan_path);std::array<int64_t,6>x;int64_t den;
 if(case_id=="low-primary"){den=100000000000000LL;x={den,99999998999999LL,99999997999988LL,99999996999867LL,57255226550170LL,43659674865LL};}
 else if(case_id=="S3-primary"){den=5;x={5,5,5,5,3,0};}else throw std::runtime_error("unregistered case");
 if((blocks!=1&&blocks!=4)||block<0||block>=blocks)throw std::runtime_error("unregistered partition");
 std::array<bool,FULL>sample{};for(uint64_t k=0;k<128;k++)sample[7559*k]=true;
 PhysicalPF engine;std::array<B,6>total;std::array<int,6>perm={0,1,2,3,4,5};uint64_t index=0,done=0,zero=0,repeated=0,per_index=0;auto tic=std::chrono::steady_clock::now();
 do{for(int mask=0;mask<64;mask++)for(auto&c:fan){bool wanted=mode=="benchmark"?sample[index]:!sample[index]&&int(per_index/(720/blocks))==block;index++;if(!wanted)continue;
  std::array<int,6>sign;for(int k=0;k<6;k++)sign[k]=(mask>>k&1)?-1:1;
  std::array<Q,6>base;std::array<std::array<B,6>,6>d;bool haszero=false,hasrepeat=false;
  for(int j=0;j<6;j++){__int128 num=0;for(int k=0;k<6;k++)num+=(__int128)c.ray[j][k]*sign[k]*x[perm[k]];base[j]=rat(checked(num),checked((__int128)den*c.q[j]));haszero|=base[j].n==0;for(int l=0;l<j;l++)hasrepeat|=base[l]==base[j];for(int k=0;k<6;k++){d[j][perm[k]]=c.d[j][k];if(sign[k]<0)arb_neg(d[j][perm[k]].v,d[j][perm[k]].v);}}
  std::array<B,6>gradient;engine.evaluate(base,d,gradient);for(int k=0;k<6;k++)arb_addmul(total[k].v,c.alpha.v,gradient[k].v,PREC);done++;zero+=haszero;repeated+=hasrepeat;
 }per_index++;}while(std::next_permutation(perm.begin(),perm.end()));
 if(index!=FULL||(mode=="benchmark"&&done!=128))throw std::runtime_error("tuple count");
 B norm;arb_const_pi(norm.v,PREC);arb_mul_2exp_si(norm.v,norm.v,1);arb_pow_ui(norm.v,norm.v,6,PREC);for(auto&v:total)arb_div(v.v,v.v,norm.v,PREC);
 std::cout<<"{\"status\":\"COMPLETE_REGISTERED_PARTIAL_GRADIENT\",\"case_id\":\""<<case_id<<"\",\"mode\":\""<<mode<<"\",\"precision_bits\":768,\"block\":"<<block<<",\"blocks\":"<<blocks<<",\"computed_tuples\":"<<done<<",\"physical_fields\":6,\"tuple_field_contributions\":"<<6*done<<",\"zero_base_tuples\":"<<zero<<",\"repeated_base_tuples\":"<<repeated<<",\"cached_nodes\":"<<engine.node_count()<<",\"kernel_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-tic).count()<<",\"gradient_partial_normalized\":";output6(total);std::cout<<",\"full_reference\":false,\"benchmark_indices\":\"7559*k,k=0..127\",\"skip_benchmark_in_remaining_mode\":true}\n";
}
int main(int argc,char**argv){try{static_assert(sizeof(arb_struct)==48&&sizeof(fmpz)==8&&sizeof(ulong)==8);if(std::string(flint_version)!="3.6.0")throw std::runtime_error("FLINT ABI");if(argc==2&&std::string(argv[1])=="scalars")scalars();else if(argc==4&&std::string(argv[1])=="benchmark")real_case("benchmark",argv[2],argv[3],0,1);else if(argc==6&&std::string(argv[1])=="remaining")real_case("remaining",argv[2],argv[3],std::stoi(argv[4]),std::stoi(argv[5]));else throw std::runtime_error("scalars | benchmark CASE FAN | remaining CASE FAN BLOCK BLOCKS");}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}return 0;}
