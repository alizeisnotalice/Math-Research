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
struct Ball {
    arb_t a;
    Ball(){arb_init(a);} ~Ball(){arb_clear(a);}
    Ball(const Ball& b){arb_init(a);arb_set(a,b.a);}
    Ball& operator=(const Ball&b){arb_set(a,b.a);return *this;}
};
struct Rat {int64_t n,d;bool operator==(const Rat&o)const{return n==o.n&&d==o.d;}};
static Rat rational(int64_t n,int64_t d){
    if(d<=0||n==std::numeric_limits<int64_t>::min())throw std::runtime_error("invalid rational range");
    int64_t g=std::gcd(n,d);return {n/g,d/g};
}
static bool lt(const Rat&a,const Rat&b){return (__int128)a.n*b.d<(__int128)b.n*a.d;}
struct RH {size_t operator()(const Rat&a)const{return std::hash<int64_t>{}(a.n)^(std::hash<int64_t>{}(a.d)+0x9e3779b97f4a7c15ULL+(uint64_t(a.n)<<1));}};
using Key=std::array<Rat,6>;
struct KH {size_t operator()(const Key&k)const{size_t h=0;for(auto&r:k)h^=RH{}(r)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
struct Node {
    Ball x;std::array<Ball,6> jet;bool ready=false;
    Node(Rat r,const Ball&T){Ball b;arb_set_si(x.a,r.n);arb_div_ui(x.a,x.a,r.d,PREC);arb_div(b.a,x.a,T.a,PREC);arb_atan(jet[0].a,b.a,PREC);}
    void jets(const Ball&T){
        if(ready)return;Ball D,tmp;std::array<Ball,5>g;
        arb_mul(D.a,x.a,x.a,PREC);arb_addmul(D.a,T.a,T.a,PREC);arb_div(g[0].a,T.a,D.a,PREC);
        for(int k=1;k<5;k++){
            arb_mul(tmp.a,x.a,g[k-1].a,PREC);arb_mul_2exp_si(tmp.a,tmp.a,1);
            if(k>=2)arb_add(tmp.a,tmp.a,g[k-2].a,PREC);
            arb_neg(tmp.a,tmp.a);arb_div(g[k].a,tmp.a,D.a,PREC);
        }
        for(int k=1;k<=5;k++)arb_div_ui(jet[k].a,g[k-1].a,k,PREC);
        ready=true;
    }
};
class DD {
    const Ball&T;std::unordered_map<Rat,std::unique_ptr<Node>,RH> cache;
public:
    explicit DD(const Ball&t):T(t){cache.reserve(50000);}
    size_t nodes()const{return cache.size();}
    void eval(Ball&out,const Key&k){
        std::array<Node*,6>n;std::array<Ball,6>b;Ball d,tmp;
        for(int j=0;j<6;j++){
            auto it=cache.find(k[j]);
            if(it==cache.end())it=cache.emplace(k[j],std::make_unique<Node>(k[j],T)).first;
            n[j]=it->second.get();arb_set(b[j].a,n[j]->jet[0].a);
        }
        for(int order=1;order<6;order++)for(int j=0;j<6-order;j++){
            if(k[j]==k[j+order]){n[j]->jets(T);arb_set(b[j].a,n[j]->jet[order].a);}
            else{
                arb_sub(d.a,n[j+order]->x.a,n[j]->x.a,PREC);
                if(arb_contains_zero(d.a))throw std::runtime_error("insufficient precision separates exact nodes");
                arb_sub(tmp.a,b[j+1].a,b[j].a,PREC);arb_div(b[j].a,tmp.a,d.a,PREC);
            }
        }
        arb_set(out.a,b[0].a);
        if(!arb_is_finite(out.a))throw std::runtime_error("nonfinite DD");
    }
};

// This is a FORMAL per-tuple epsilon^0 coefficient. It is not physical G
// until whole-fan, finite-part theorem and off-true-S2 gates are discharged.
void formal_G_constant(Ball&out,Key k){
 std::sort(k.begin(),k.end(),lt);std::array<Ball,6>b,x;Ball pi,d,tmp;
 arb_const_pi(pi.a,PREC);arb_mul_2exp_si(pi.a,pi.a,-1);
 for(int j=0;j<6;j++){
  arb_set_si(x[j].a,k[j].n);arb_div_ui(x[j].a,x[j].a,k[j].d,PREC);
  if(k[j].n>0)arb_set(b[j].a,pi.a);else if(k[j].n<0)arb_neg(b[j].a,pi.a);else arb_zero(b[j].a);
 }
 for(int order=1;order<6;order++)for(int j=0;j<6-order;j++){
  if(k[j]==k[j+order])arb_zero(b[j].a); // constant coefficient of every higher jet
  else{arb_sub(d.a,x[j+order].a,x[j].a,PREC);if(arb_contains_zero(d.a))throw std::runtime_error("distinct exact nodes inseparable");arb_sub(tmp.a,b[j+1].a,b[j].a,PREC);arb_div(b[j].a,tmp.a,d.a,PREC);}
 }
 arb_set(out.a,b[0].a);
}
// Exported scalar-only core ensures the real Arb code is compiled. The ABI
// executable deliberately provides no positive-time scalar or full-frame CLI.
extern "C" void codex_sparse4_FORMAL_constant_and_H(const int64_t*n,const int64_t*d,const arb_t time,arb_t formal_constant,arb_t H){
 if(!arb_is_positive(time))throw std::runtime_error("require positive common time");
 Key z;for(int i=0;i<6;i++)z[i]=rational(n[i],d[i]);std::sort(z.begin(),z.end(),lt);
 Ball T,g,h;arb_set(T.a,time);DD dd(T);dd.eval(h,z);formal_G_constant(g,z);arb_set(formal_constant,g.a);arb_set(H,h.a);
}
int main(int argc,char**argv){
 try{static_assert(sizeof(ulong)==8&&sizeof(fmpz)==8&&sizeof(arb_struct)==48);if(std::string(flint_version)!="3.6.0")throw std::runtime_error("FLINT ABI version mismatch");
 if(argc!=2||std::string(argv[1])!="abi")throw std::runtime_error("PREPARATION ONLY: no full frame or numeric scalar CLI; F050/B052/off-S2/main release gates pending");
 std::cout<<"{\"mode\":\"ABI_ONLY_NO_KERNEL\",\"flint\":\""<<flint_version<<"\",\"ulong_bytes\":"<<sizeof(ulong)<<",\"fmpz_bytes\":"<<sizeof(fmpz)<<",\"arb_struct_bytes\":"<<sizeof(arb_struct)<<",\"precision_bits\":"<<PREC<<",\"full_kernel_calls\":0,\"positive_time_evaluations\":0}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}return 0;
}
