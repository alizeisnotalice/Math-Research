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
static Key canonical(Key k){
    std::sort(k.begin(),k.end(),lt);Key neg;
    for(int i=0;i<6;i++)neg[i]={-k[5-i].n,k[5-i].d};
    if(std::lexicographical_compare(neg.begin(),neg.end(),k.begin(),k.end(),lt))return neg;
    return k;
}
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
static void print_interval(const Ball&b){
    fmpz_t l,u,e;fmpz_init(l);fmpz_init(u);fmpz_init(e);
    arb_get_interval_fmpz_2exp(l,u,e,b.a);
    char*ls=fmpz_get_str(nullptr,10,l),*us=fmpz_get_str(nullptr,10,u),*es=fmpz_get_str(nullptr,10,e);
    std::cout<<"{\"lower_mantissa\":\""<<ls<<"\",\"upper_mantissa\":\""<<us<<"\",\"common_binary_exponent\":\""<<es<<"\"}";
    flint_free(ls);flint_free(us);flint_free(es);fmpz_clear(l);fmpz_clear(u);fmpz_clear(e);
}
struct Ray {int64_t q;std::array<int64_t,6>v;};
struct Cone {int64_t det;std::array<Ray,6>r;Ball alpha;};
static std::vector<Cone> load(const char*path){
    std::ifstream f(path);int count;if(!(f>>count)||count!=21)throw std::runtime_error("expected21fan cones");
    std::vector<Cone>c(count);
    for(auto&v:c){f>>v.det;arb_set_si(v.alpha.a,v.det);for(auto&r:v.r){f>>r.q;if(r.q<=0)throw std::runtime_error("rate");for(auto&x:r.v)f>>x;arb_div_ui(v.alpha.a,v.alpha.a,r.q,PREC);}}
    if(!f)throw std::runtime_error("fan parse");return c;
}
static int64_t checked(__int128 n){if(n<=std::numeric_limits<int64_t>::min()||n>std::numeric_limits<int64_t>::max())throw std::runtime_error("exact integer range exceeded");return(int64_t)n;}
static void anchors(){
    Ball T;arb_one(T.a);arb_div_ui(T.a,T.a,250,PREC);DD dd(T);
    std::array<Key,6> a={{{{{-2,1},{-1,1},{0,1},{1,1},{2,1},{3,1}}},{{{-2,1},{-1,1},{0,1},{0,1},{2,1},{3,1}}},{{{0,1},{0,1},{0,1},{0,1},{0,1},{0,1}}},{{{1,3},{1,3},{1,3},{1,3},{1,3},{1,3}}},{{{-1,1},{-1,1},{0,1},{0,1},{1,1},{1,1}}},{{{-1,7},{-1,7},{2,9},{2,9},{2,9},{4,5}}}}};
    std::cout<<"{\"flint_version\":\""<<flint_version<<"\",\"precision_bits\":"<<PREC<<",\"sizeof_arb\":"<<sizeof(arb_struct)<<",\"anchors\":[";
    for(int i=0;i<6;i++){Ball v;dd.eval(v,canonical(a[i]));if(i)std::cout<<",";print_interval(v);}std::cout<<"]}\n";
}
static void full(int argc,char**argv){
    if(argc!=12)throw std::runtime_error("full FAN T_NUM T_DEN X_DEN X0 X1 X2 X3 X4 X5");
    auto cones=load(argv[2]);Ball T;arb_set_si(T.a,std::stoll(argv[3]));arb_div_ui(T.a,T.a,std::stoll(argv[4]),PREC);
    int64_t den=std::stoll(argv[5]);std::array<int64_t,6>x;for(int i=0;i<6;i++)x[i]=std::stoll(argv[6+i]);
    if(den<=0||arb_is_nonpositive(T.a))throw std::runtime_error("positive den/T required");
    auto start=std::chrono::steady_clock::now();
    DD dd(T);Ball total,value,tmp;std::vector<Ball>alpha2(21);
    for(int c=0;c<21;c++)arb_mul_2exp_si(alpha2[c].a,cones[c].alpha.a,1);
    std::array<int,6>perm={0,1,2,3,4,5};uint64_t terms=0,repeat_terms=0,zero_terms=0;
    uint64_t computed=0;
    // mask^63 is the exact global negative. DD5[atan] is unchanged by
    // negating all six nodes. Both full-group weights are retained in alpha2.
    do{for(unsigned mask=0;mask<32;mask++){
        std::array<int64_t,6>y;for(int j=0;j<6;j++)y[j]=((mask>>j)&1)?checked(-(__int128)x[perm[j]]):x[perm[j]];
        for(int c=0;c<21;c++){
            Key k;bool zero=false;
            for(int j=0;j<6;j++){__int128 n=0;for(int q=0;q<6;q++)n+=(__int128)y[q]*cones[c].r[j].v[q];k[j]=rational(checked(n),checked((__int128)den*cones[c].r[j].q));zero|=k[j].n==0;}
            k=canonical(k);bool repeated=false;for(int j=1;j<6;j++)repeated|=k[j]==k[j-1];
            dd.eval(value,k);arb_addmul(total.a,alpha2[c].a,value.a,PREC);
            terms+=2;computed++;repeat_terms+=2*repeated;zero_terms+=2*zero;
        }
    }}while(std::next_permutation(perm.begin(),perm.end()));
    if(terms!=967680)throw std::runtime_error("logical group count");
    arb_const_pi(tmp.a,PREC);arb_mul_2exp_si(tmp.a,tmp.a,1);arb_pow_ui(tmp.a,tmp.a,6,PREC);arb_div(total.a,total.a,tmp.a,PREC);
    if(!arb_is_finite(total.a))throw std::runtime_error("nonfinite full sum");
    auto finish=std::chrono::steady_clock::now();
    std::cout<<"{\"flint_version\":\""<<flint_version<<"\",\"precision_bits\":"<<PREC<<",\"full_group_terms\":"<<terms<<",\"zero_node_terms\":"<<zero_terms<<",\"confluent_terms\":"<<repeat_terms<<",\"computed_DD_evaluations\":"<<computed<<",\"unique_node_values\":"<<dd.nodes()<<",\"exact_global_sign_pairing\":true,\"total_kernel_seconds\":"<<std::chrono::duration<double>(finish-start).count()<<",\"dyadic_interval\":";
    print_interval(total);std::cout<<"}\n";
}
int main(int argc,char**argv){
    try{static_assert(sizeof(ulong)==8&&sizeof(fmpz)==8&&sizeof(arb_struct)==48);if(std::string(flint_version)!="3.6.0")throw std::runtime_error("FLINT ABI version mismatch");
        if(argc==2&&std::string(argv[1])=="anchors")anchors();else if(argc>1&&std::string(argv[1])=="full")full(argc,argv);else throw std::runtime_error("anchors|full mode required");
    }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}return 0;
}
