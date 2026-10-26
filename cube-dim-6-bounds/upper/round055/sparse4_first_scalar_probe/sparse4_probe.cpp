#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include "arb.h"
#include "fmpz.h"
static constexpr slong PREC=512;
struct Ball{arb_t a;Ball(){arb_init(a);}~Ball(){arb_clear(a);}Ball(const Ball&b){arb_init(a);arb_set(a,b.a);}Ball&operator=(const Ball&b){arb_set(a,b.a);return*this;}};
static int64_t integer(const std::string&s){if(s.empty()||s=="-0"||s[0]=='+'||(s[0]=='0'&&s.size()>1)||(s[0]=='-'&&(s.size()<2||s[1]=='0')))throw std::runtime_error("noncanonical integer");int64_t n=0;auto r=std::from_chars(s.data(),s.data()+s.size(),n);if(r.ec!=std::errc()||r.ptr!=s.data()+s.size())throw std::runtime_error("invalid integer");return n;}
static int64_t next(std::istream&in){std::string s;if(!(in>>s))throw std::runtime_error("truncated input");return integer(s);}
static int64_t checked(__int128 v){if(v<=std::numeric_limits<int64_t>::min()||v>std::numeric_limits<int64_t>::max())throw std::runtime_error("int64 exact range");return(int64_t)v;}
struct Rat{int64_t n,d;bool operator==(const Rat&r)const{return n==r.n&&d==r.d;}};
static Rat rational(int64_t n,int64_t d){if(d<=0||n==std::numeric_limits<int64_t>::min())throw std::runtime_error("rational domain");int64_t g=std::gcd(n,d);return{n/g,d/g};}
static bool less(Rat a,Rat b){return(__int128)a.n*b.d<(__int128)b.n*a.d;}
struct RH{size_t operator()(Rat r)const{return std::hash<int64_t>{}(r.n)^(std::hash<int64_t>{}(r.d)+0x9e3779b97f4a7c15ULL+(uint64_t(r.n)<<1));}};
static void interval(std::ostream&os,const Ball&v){if(!arb_is_finite(v.a))throw std::runtime_error("nonfinite output");fmpz_t l,u,e;fmpz_init(l);fmpz_init(u);fmpz_init(e);arb_get_interval_fmpz_2exp(l,u,e,v.a);char*ls=fmpz_get_str(nullptr,10,l),*us=fmpz_get_str(nullptr,10,u),*es=fmpz_get_str(nullptr,10,e);os<<"{\"lower_mantissa\":\""<<ls<<"\",\"upper_mantissa\":\""<<us<<"\",\"common_binary_exponent\":\""<<es<<"\"}";flint_free(ls);flint_free(us);flint_free(es);fmpz_clear(l);fmpz_clear(u);fmpz_clear(e);}
struct Node{Ball x,G0;std::array<Ball,6>H;bool ready=false;
 Node(Rat r,const Ball&T,const Ball&halfpi){Ball y;arb_set_si(x.a,r.n);arb_div_ui(x.a,x.a,r.d,PREC);arb_div(y.a,x.a,T.a,PREC);arb_atan(H[0].a,y.a,PREC);if(r.n>0)arb_set(G0.a,halfpi.a);else if(r.n<0)arb_neg(G0.a,halfpi.a);}
 void jets(const Ball&T){if(ready)return;Ball D,tmp;std::array<Ball,5>g;arb_mul(D.a,x.a,x.a,PREC);arb_addmul(D.a,T.a,T.a,PREC);arb_div(g[0].a,T.a,D.a,PREC);for(int k=1;k<5;k++){arb_mul(tmp.a,x.a,g[k-1].a,PREC);arb_mul_2exp_si(tmp.a,tmp.a,1);if(k>=2)arb_add(tmp.a,tmp.a,g[k-2].a,PREC);arb_neg(tmp.a,tmp.a);arb_div(g[k].a,tmp.a,D.a,PREC);}for(int k=1;k<=5;k++)arb_div_ui(H[k].a,g[k-1].a,k,PREC);ready=true;}
};
class JointDD{const Ball&T;Ball halfpi;std::unordered_map<Rat,std::unique_ptr<Node>,RH>cache;
 Node&node(Rat r){auto i=cache.find(r);if(i==cache.end())i=cache.emplace(r,std::make_unique<Node>(r,T,halfpi)).first;return*i->second;}
public:explicit JointDD(const Ball&t):T(t){arb_const_pi(halfpi.a,PREC);arb_mul_2exp_si(halfpi.a,halfpi.a,-1);cache.reserve(131072);}size_t nodes()const{return cache.size();}
 void evaluate(Ball&G,Ball&H,const std::array<Rat,6>&z){std::array<Node*,6>n;std::array<Ball,6>g,h;Ball d,inv,tmp;
  for(int j=0;j<6;j++){n[j]=&node(z[j]);arb_set(g[j].a,n[j]->G0.a);arb_set(h[j].a,n[j]->H[0].a);}
  for(int k=1;k<6;k++)for(int j=0;j<6-k;j++){
   if(z[j]==z[j+k]){n[j]->jets(T);arb_set(h[j].a,n[j]->H[k].a);arb_zero(g[j].a);}
   else{arb_sub(d.a,n[j+k]->x.a,n[j]->x.a,PREC);if(arb_contains_zero(d.a))throw std::runtime_error("distinct exact nodes inseparable");arb_inv(inv.a,d.a,PREC);arb_sub(tmp.a,h[j+1].a,h[j].a,PREC);arb_mul(h[j].a,tmp.a,inv.a,PREC);arb_sub(tmp.a,g[j+1].a,g[j].a,PREC);arb_mul(g[j].a,tmp.a,inv.a,PREC);}
  }arb_set(G.a,g[0].a);arb_set(H.a,h[0].a);if(!arb_is_finite(G.a)||!arb_is_finite(H.a))throw std::runtime_error("nonfinite DD field");
 }
};
struct Ray{int64_t q;std::array<int64_t,6>r;};
struct Cone{int64_t det;std::array<int,6>ids;Ball alpha24;};
struct Input{std::vector<Ray>rays;std::vector<Cone>cones;std::vector<std::array<int64_t,6>>images;};
static Input load(const char*path){std::ifstream in(path);if(!in.is_open())throw std::runtime_error("input unavailable");std::string header;if(!(in>>header)||header!="C055_SPARSE4_BETA1_80_X050_T32_825")throw std::runtime_error("input model header");if(next(in)!=478||next(in)!=13229||next(in)!=1920)throw std::runtime_error("complete input counts");if(next(in)!=50||next(in)!=24||next(in)!=32||next(in)!=825)throw std::runtime_error("fixed point/time/orbit protocol");
 Input d;d.rays.resize(478);d.cones.resize(13229);d.images.resize(1920);
 for(auto&v:d.rays){v.q=next(in);if(v.q<=0||v.q>20773)throw std::runtime_error("rate range");bool nonzero=false;for(auto&x:v.r){x=next(in);if(x<0||x>520)throw std::runtime_error("ray range");nonzero|=x!=0;}if(!nonzero)throw std::runtime_error("zero ray");}
 for(auto&c:d.cones){c.det=next(in);if(c.det<=0||c.det>2048000000LL)throw std::runtime_error("det range");std::set<int>ids;for(auto&i:c.ids){int64_t k=next(in);if(k<0||k>=478)throw std::runtime_error("ray id");i=(int)k;ids.insert(i);}if(ids.size()!=6)throw std::runtime_error("repeated cone ray");}
 std::set<std::array<int64_t,6>>seen;const std::array<int64_t,6>expected={0,27,45,50,50,50};
 for(auto&y:d.images){std::array<int64_t,6>abs,neg;for(int i=0;i<6;i++){y[i]=next(in);if(y[i]<-50||y[i]>50)throw std::runtime_error("image range");abs[i]=std::abs(y[i]);neg[i]=-y[i];}std::sort(abs.begin(),abs.end());if(abs!=expected||!(y<neg)||!seen.insert(y).second)throw std::runtime_error("orbit uniqueness/parity");}
 std::string extra;if(in>>extra)throw std::runtime_error("extra input tokens");if(in.bad())throw std::runtime_error("input IO failure");return d;
}
static void run_full(const char*path,const char*cpPath,const char*tokenPath){std::ifstream gate(tokenPath);std::string token,extra;if(!(gate>>token)||token!="C055_MAIN_RELEASED_FINAL_DEPENDENCIES_CHECKED"||(gate>>extra))throw std::runtime_error("full not released");std::ifstream prior(cpPath);if(prior.is_open())throw std::runtime_error("checkpoint exists;refuse repeat");
 auto in=load(path);auto tic=std::chrono::steady_clock::now();Ball T;arb_set_ui(T.a,32);arb_div_ui(T.a,T.a,825,PREC);
 for(auto&c:in.cones){arb_set_si(c.alpha24.a,c.det);arb_mul_ui(c.alpha24.a,c.alpha24.a,24,PREC);for(int id:c.ids)arb_div_ui(c.alpha24.a,c.alpha24.a,in.rays[id].q,PREC);}
 JointDD dd(T);Ball g,h,G,H,K,scale;uint64_t terms=0,zeros=0,repeated=0;std::ofstream cp(cpPath,std::ios::out);if(!cp)throw std::runtime_error("checkpoint open failure");
 for(size_t image=0;image<in.images.size();image++){
  auto&y=in.images[image];for(auto&c:in.cones){std::array<Rat,6>z;bool zero=false;for(int j=0;j<6;j++){auto&r=in.rays[c.ids[j]];__int128 n=0;for(int k=0;k<6;k++)n+=(__int128)y[k]*r.r[k];z[j]=rational(checked(n),checked((__int128)50*r.q));zero|=z[j].n==0;}
   std::sort(z.begin(),z.end(),less);bool rep=false;for(int j=1;j<6;j++)rep|=z[j]==z[j-1];dd.evaluate(g,h,z);arb_addmul(G.a,c.alpha24.a,g.a,PREC);arb_addmul(H.a,c.alpha24.a,h.a,PREC);terms++;zeros+=zero;repeated+=rep;
  }
  if((image+1)%32==0||image+1==in.images.size()){auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-tic).count();cp<<"{\"phase\":\"running\",\"complete\":false,\"images_complete\":"<<image+1<<",\"paired_image_cone_terms\":"<<terms<<",\"logical_terms_each_field\":"<<terms*24<<",\"G_DD5_evaluations\":"<<terms<<",\"H_DD5_evaluations\":"<<terms<<",\"total_DD5_evaluations\":"<<2*terms<<",\"unique_nodes\":"<<dd.nodes()<<",\"seconds\":"<<elapsed<<"}\n";cp.flush();if(!cp)throw std::runtime_error("checkpoint write failure");}
 }
 if(terms!=25399680ULL||terms*24!=609592320ULL)throw std::runtime_error("full field count mismatch");arb_const_pi(scale.a,PREC);arb_mul_2exp_si(scale.a,scale.a,1);arb_pow_ui(scale.a,scale.a,6,PREC);arb_div(G.a,G.a,scale.a,PREC);arb_div(H.a,H.a,scale.a,PREC);arb_sub(K.a,G.a,H.a,PREC);if(!arb_is_finite(G.a)||!arb_is_finite(H.a)||!arb_is_finite(K.a))throw std::runtime_error("nonfinite full result");auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-tic).count();
 std::cout<<"{\"status\":\"COMPLETE_ONE_JOINT_SCALAR_PROBE\",\"precision_bits\":512,\"paired_image_cone_terms\":"<<terms<<",\"G_logical_terms\":"<<terms*24<<",\"H_logical_terms\":"<<terms*24<<",\"G_DD5_evaluations\":"<<terms<<",\"H_DD5_evaluations\":"<<terms<<",\"total_DD5_evaluations\":"<<2*terms<<",\"zero_node_paired_terms\":"<<zeros<<",\"confluent_paired_terms\":"<<repeated<<",\"unique_nodes\":"<<dd.nodes()<<",\"kernel_seconds\":"<<elapsed<<",\"G\":";interval(std::cout,G);std::cout<<",\"H_T\":";interval(std::cout,H);std::cout<<",\"K_T\":";interval(std::cout,K);std::cout<<"}\n";std::cout.flush();if(!std::cout)throw std::runtime_error("stdout failure");cp<<"{\"phase\":\"computed\",\"complete\":true,\"paired_image_cone_terms\":"<<terms<<",\"G_DD5_evaluations\":"<<terms<<",\"H_DD5_evaluations\":"<<terms<<",\"total_DD5_evaluations\":"<<2*terms<<",\"seconds\":"<<elapsed<<"}\n";cp.flush();if(!cp)throw std::runtime_error("final checkpoint failure");
}
int main(int argc,char**argv){try{static_assert(sizeof(ulong)==8&&sizeof(fmpz)==8&&sizeof(arb_struct)==48);if(std::string(flint_version)!="3.6.0")throw std::runtime_error("FLINT ABI mismatch");
 if(argc==2&&std::string(argv[1])=="abi")std::cout<<"{\"mode\":\"ABI_ONLY\",\"flint\":\""<<flint_version<<"\",\"ulong\":"<<sizeof(ulong)<<",\"fmpz\":"<<sizeof(fmpz)<<",\"arb_struct\":"<<sizeof(arb_struct)<<",\"positive_time_evaluations\":0,\"full_calls\":0}\n";
 else if(argc==3&&std::string(argv[1])=="validate"){auto in=load(argv[2]);std::cout<<"{\"mode\":\"VALIDATE_ONLY_NO_KERNEL\",\"rays\":"<<in.rays.size()<<",\"cones\":"<<in.cones.size()<<",\"images\":"<<in.images.size()<<",\"paired_terms\":25399680,\"logical_terms_each_field\":609592320,\"total_DD5_evaluations_if_released\":50799360,\"positive_time_evaluations\":0,\"full_calls\":0}\n";}
 else if(argc==5&&std::string(argv[1])=="full")run_full(argv[2],argv[3],argv[4]);else throw std::runtime_error("abi|validate INPUT|full INPUT CHECKPOINT RELEASE_TOKEN");
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}return 0;}
