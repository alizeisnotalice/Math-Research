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
static constexpr slong PREC=512;
struct Ball{arb_t a;Ball(){arb_init(a);}~Ball(){arb_clear(a);}Ball(const Ball&b){arb_init(a);arb_set(a,b.a);}Ball&operator=(const Ball&b){arb_set(a,b.a);return*this;}};
struct Rat{int64_t n,d;bool operator==(const Rat&o)const{return n==o.n&&d==o.d;}};
struct ZeroNode:std::runtime_error{ZeroNode():std::runtime_error("exact zero Fourier node; whole evaluation unsupported"){} };
static int64_t checked(__int128 n){if(n<=std::numeric_limits<int64_t>::min()||n>std::numeric_limits<int64_t>::max())throw std::runtime_error("int64 exact range exceeded");return(int64_t)n;}
static Rat rat(int64_t n,int64_t d){if(d<=0||n==std::numeric_limits<int64_t>::min())throw std::runtime_error("rational range");int64_t g=std::gcd(n,d);return{n/g,d/g};}
static bool lt(Rat a,Rat b){return(__int128)a.n*b.d<(__int128)b.n*a.d;}
struct RH{size_t operator()(Rat r)const{return std::hash<int64_t>{}(r.n)^(std::hash<int64_t>{}(r.d)+0x9e3779b97f4a7c15ULL+(uint64_t(r.n)<<1));}};
static void endpoints(const Ball&b){fmpz_t l,u,e;fmpz_init(l);fmpz_init(u);fmpz_init(e);arb_get_interval_fmpz_2exp(l,u,e,b.a);char*ls=fmpz_get_str(nullptr,10,l),*us=fmpz_get_str(nullptr,10,u),*es=fmpz_get_str(nullptr,10,e);std::cout<<"{\"lower_mantissa\":\""<<ls<<"\",\"upper_mantissa\":\""<<us<<"\",\"common_binary_exponent\":\""<<es<<"\"}";flint_free(ls);flint_free(us);flint_free(es);fmpz_clear(l);fmpz_clear(u);fmpz_clear(e);}
struct Node{
 Ball x,D,re,im;std::array<Ball,7>jet;bool ready=false;
 Node(Rat r,const Ball&T){if(r.n==0)throw ZeroNode();Ball tmp;arb_set_si(x.a,r.n);arb_div_ui(x.a,x.a,r.d,PREC);arb_mul(D.a,x.a,x.a,PREC);arb_addmul(D.a,T.a,T.a,PREC);arb_div(re.a,T.a,D.a,PREC);arb_div(im.a,x.a,D.a,PREC);arb_div(tmp.a,T.a,x.a,PREC);arb_atan(jet[0].a,tmp.a,PREC);}
 void jets(const Ball&T){if(ready)return;std::array<Ball,6>g;Ball tmp;arb_set(g[0].a,re.a);for(int k=1;k<6;k++){arb_mul(tmp.a,x.a,g[k-1].a,PREC);arb_mul_2exp_si(tmp.a,tmp.a,1);if(k>=2)arb_add(tmp.a,tmp.a,g[k-2].a,PREC);arb_neg(tmp.a,tmp.a);arb_div(g[k].a,tmp.a,D.a,PREC);}for(int k=1;k<=6;k++){arb_div_ui(jet[k].a,g[k-1].a,k,PREC);arb_neg(jet[k].a,jet[k].a);}ready=true;}
};
class ShortEngine{
 const Ball&T;std::unordered_map<Rat,std::unique_ptr<Node>,RH>nodes;
 Node&node(Rat r){auto i=nodes.find(r);if(i==nodes.end())i=nodes.emplace(r,std::make_unique<Node>(r,T)).first;return*i->second;}
 void hermite(Ball&out,std::vector<Rat>z){std::sort(z.begin(),z.end(),lt);std::array<Ball,7>b;Ball d,tmp;int n=z.size();for(int j=0;j<n;j++){if(z[j].n==0)throw ZeroNode();arb_set(b[j].a,node(z[j]).jet[0].a);}for(int k=1;k<n;k++)for(int j=0;j<n-k;j++){if(z[j]==z[j+k]){auto&v=node(z[j]);v.jets(T);arb_set(b[j].a,v.jet[k].a);}else{arb_sub(d.a,node(z[j+k]).x.a,node(z[j]).x.a,PREC);if(arb_contains_zero(d.a))throw std::runtime_error("precision cannot separate distinct exact nodes");arb_sub(tmp.a,b[j+1].a,b[j].a,PREC);arb_div(b[j].a,tmp.a,d.a,PREC);}}arb_set(out.a,b[0].a);}
public:
 explicit ShortEngine(const Ball&t):T(t){nodes.reserve(50000);}size_t node_count()const{return nodes.size();}
 bool evaluate(const std::array<Rat,6>&z,Ball&value,std::array<Ball,6>&partial,Ball&density){
  arb_zero(value.a);for(auto&v:partial)arb_zero(v.a);bool repeated=false;
  for(int i=0;i<6;i++){if(z[i].n==0)throw ZeroNode();for(int j=0;j<i;j++)repeated|=z[i]==z[j];}
  std::array<Node*,6>x;for(int i=0;i<6;i++)x[i]=&node(z[i]);
  Ball pr,pi,nr,ni,tmp;arb_one(pr.a);
  for(int i=0;i<6;i++){arb_mul(nr.a,pr.a,x[i]->re.a,PREC);arb_submul(nr.a,pi.a,x[i]->im.a,PREC);arb_mul(ni.a,pr.a,x[i]->im.a,PREC);arb_addmul(ni.a,pi.a,x[i]->re.a,PREC);arb_set(pr.a,nr.a);arb_set(pi.a,ni.a);}arb_set(density.a,pr.a);
  if(repeated){std::vector<Rat>v(z.begin(),z.end());hermite(value,v);for(int i=0;i<6;i++){auto q=v;q.push_back(z[i]);hermite(partial[i],q);}return true;}
  std::array<Ball,6>prod,w,c;std::array<Ball,15>diff;for(auto&v:prod)arb_one(v.a);
  int k=0;for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,k++){arb_sub(diff[k].a,x[i]->x.a,x[j]->x.a,PREC);if(arb_contains_zero(diff[k].a))throw std::runtime_error("precision cannot separate distinct exact nodes");arb_mul(prod[i].a,prod[i].a,diff[k].a,PREC);arb_neg(tmp.a,diff[k].a);arb_mul(prod[j].a,prod[j].a,tmp.a,PREC);}
  for(int i=0;i<6;i++){arb_inv(w[i].a,prod[i].a,PREC);arb_mul(c[i].a,x[i]->jet[0].a,w[i].a,PREC);arb_add(value.a,value.a,c[i].a,PREC);arb_mul(partial[i].a,x[i]->re.a,w[i].a,PREC);arb_neg(partial[i].a,partial[i].a);}
  k=0;for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,k++){arb_add(tmp.a,c[i].a,c[j].a,PREC);arb_div(tmp.a,tmp.a,diff[k].a,PREC);arb_sub(partial[i].a,partial[i].a,tmp.a,PREC);arb_add(partial[j].a,partial[j].a,tmp.a,PREC);}return false;
 }
};
struct Ray{int64_t q;std::array<int64_t,6>v;};struct Cone{int64_t det;std::array<Ray,6>r;Ball alpha2;};
static std::vector<Cone>fan(const char*path){std::ifstream in(path);int n;in>>n;if(n!=21)throw std::runtime_error("21 fan cells required");std::vector<Cone>c(n);for(auto&v:c){in>>v.det;if(v.det<=0||v.det>1000000000)throw std::runtime_error("fan determinant range");arb_set_si(v.alpha2.a,v.det);arb_mul_2exp_si(v.alpha2.a,v.alpha2.a,1);for(auto&r:v.r){in>>r.q;if(r.q<=0||r.q>1000000)throw std::runtime_error("fan rate range");for(auto&a:r.v){in>>a;if(a<-1000000||a>1000000)throw std::runtime_error("fan ray range");}arb_div_ui(v.alpha2.a,v.alpha2.a,r.q,PREC);}}if(!in)throw std::runtime_error("fan parse");return c;}
static void anchors(){Ball T;arb_set_ui(T.a,4);arb_div_ui(T.a,T.a,125,PREC);ShortEngine eng(T);std::array<std::array<Rat,6>,6>sets={{{{{3,1},{-2,1},{-3,1},{1,1},{-1,1},{2,1}}},{{{-2,1},{-1,1},{-1,1},{1,1},{2,1},{3,1}}},{{{1,1},{1,1},{1,1},{1,1},{1,1},{1,1}}},{{{-2,1},{-2,1},{-2,1},{1,1},{1,1},{3,1}}},{{{-1,7},{-1,7},{2,9},{2,9},{2,9},{4,5}}},{{{-2,1},{-1,1},{0,1},{1,1},{2,1},{3,1}}}}};std::cout<<"{\"T\":\"4/125\",\"anchors\":[";for(int i=0;i<6;i++){if(i)std::cout<<",";Ball v,p;std::array<Ball,6>d;try{eng.evaluate(sets[i],v,d,p);std::cout<<"{\"value\":";endpoints(v);std::cout<<",\"partials\":[";for(int j=0;j<6;j++){if(j)std::cout<<",";endpoints(d[j]);}std::cout<<"],\"p_cone\":";endpoints(p);std::cout<<"}";}catch(const ZeroNode&){std::cout<<"{\"unsupported_zero_node\":true}";}}std::cout<<"]}\n";}
static void full(int argc,char**argv){if(argc!=12)throw std::runtime_error("full FAN T_NUM T_DEN DEN X0 X1 X2 X3 X4 X5");auto c=fan(argv[2]);int64_t tn=std::stoll(argv[3]),td=std::stoll(argv[4]),den=std::stoll(argv[5]);if(tn<=0||td<=0||den<=0)throw std::runtime_error("positive T and coordinate denominator required");Ball T;arb_set_si(T.a,tn);arb_div_ui(T.a,T.a,td,PREC);std::array<int64_t,6>x;for(int i=0;i<6;i++){x[i]=std::stoll(argv[6+i]);if(x[i]==std::numeric_limits<int64_t>::min())throw std::runtime_error("point int64 range");}
 auto tic=std::chrono::steady_clock::now();ShortEngine eng(T);Ball value,density,total,pTotal,tmp,scale,euler;std::array<Ball,6>partial,divided,gy,grad;std::array<int,6>perm={0,1,2,3,4,5};uint64_t terms=0,computed=0,repeated=0;
 do{for(unsigned mask=0;mask<32;mask++){std::array<int,6>s;std::array<int64_t,6>y;for(int j=0;j<6;j++){s[j]=(mask>>j&1)?-1:1;y[j]=checked((__int128)s[j]*x[perm[j]]);}for(auto&cone:c){std::array<Rat,6>z;for(int m=0;m<6;m++){__int128 n=0;for(int j=0;j<6;j++)n+=(__int128)y[j]*cone.r[m].v[j];z[m]=rat(checked(n),checked((__int128)den*cone.r[m].q));if(z[m].n==0)throw ZeroNode();}bool rep=eng.evaluate(z,value,partial,density);repeated+=2*rep;arb_addmul(total.a,cone.alpha2.a,value.a,PREC);arb_addmul(pTotal.a,cone.alpha2.a,density.a,PREC);
 for(auto&g:gy)arb_zero(g.a);for(int m=0;m<6;m++){arb_div_ui(divided[m].a,partial[m].a,cone.r[m].q,PREC);for(int j=0;j<6;j++)if(cone.r[m].v[j])arb_addmul_si(gy[j].a,divided[m].a,cone.r[m].v[j],PREC);}for(int j=0;j<6;j++){if(s[j]>0)arb_addmul(grad[perm[j]].a,cone.alpha2.a,gy[j].a,PREC);else arb_submul(grad[perm[j]].a,cone.alpha2.a,gy[j].a,PREC);}terms+=2;computed++;
 }}}while(std::next_permutation(perm.begin(),perm.end()));
 if(terms!=967680)throw std::runtime_error("complete group count");arb_const_pi(scale.a,PREC);arb_mul_2exp_si(scale.a,scale.a,1);arb_pow_ui(scale.a,scale.a,6,PREC);arb_div(total.a,total.a,scale.a,PREC);arb_div(pTotal.a,pTotal.a,scale.a,PREC);for(auto&g:grad)arb_div(g.a,g.a,scale.a,PREC);arb_mul_ui(euler.a,total.a,5,PREC);arb_addmul(euler.a,T.a,pTotal.a,PREC);for(int i=0;i<6;i++){arb_set_si(tmp.a,x[i]);arb_div_ui(tmp.a,tmp.a,den,PREC);arb_addmul(euler.a,tmp.a,grad[i].a,PREC);}if(!arb_is_finite(total.a)||!arb_is_finite(pTotal.a)||!arb_contains_zero(euler.a))throw std::runtime_error("finite/Euler failure");for(auto&g:grad)if(!arb_is_finite(g.a))throw std::runtime_error("nonfinite gradient");auto end=std::chrono::steady_clock::now();
 std::cout<<"{\"status\":\"generic_nonzero_supported\",\"precision_bits\":"<<PREC<<",\"full_group_terms\":"<<terms<<",\"computed_sign_pairs\":"<<computed<<",\"confluent_terms\":"<<repeated<<",\"zero_node_terms\":0,\"unique_node_values\":"<<eng.node_count()<<",\"kernel_seconds\":"<<std::chrono::duration<double>(end-tic).count()<<",\"K\":";endpoints(total);std::cout<<",\"gradient\":[";for(int i=0;i<6;i++){if(i)std::cout<<",";endpoints(grad[i]);}std::cout<<"],\"p_T\":";endpoints(pTotal);std::cout<<",\"Euler_x_grad_plus5K_plusTp\":";endpoints(euler);std::cout<<"}\n";
}
int main(int argc,char**argv){try{static_assert(sizeof(ulong)==8&&sizeof(fmpz)==8&&sizeof(arb_struct)==48);if(std::string(flint_version)!="3.6.0")throw std::runtime_error("FLINT ABI");if(argc==2&&std::string(argv[1])=="anchors")anchors();else if(argc>1&&std::string(argv[1])=="full")full(argc,argv);else throw std::runtime_error("anchors|full required");}catch(const ZeroNode&e){std::cout<<"{\"status\":\"unsupported_zero_node\",\"full_group_completed\":false}\n";std::cerr<<e.what()<<"\n";return 3;}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}return 0;}
