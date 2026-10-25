#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <chrono>
#include "arb.h"
#include "fmpz.h"
constexpr slong P=768;
struct B{arb_t a;B(){arb_init(a);}~B(){arb_clear(a);}B(const B&b){arb_init(a);arb_set(a,b.a);}B&operator=(const B&b){arb_set(a,b.a);return *this;}};
using Key=std::pair<int64_t,int64_t>;
struct N{B b,f,fp,fpp;};
struct C{long det;std::array<long,6>q;std::array<std::array<long,6>,6>r;B weight;std::array<std::array<B,6>,6>d;};
void rational(B&b,int64_t n,int64_t d){if(d<=0)throw std::runtime_error("bad denominator");arb_set_si(b.a,n);arb_div_ui(b.a,b.a,d,P);}
void output(const B&b){if(!arb_is_finite(b.a))throw std::runtime_error("nonfinite result");fmpz_t l,u,e;fmpz_init(l);fmpz_init(u);fmpz_init(e);arb_get_interval_fmpz_2exp(l,u,e,b.a);char*x=fmpz_get_str(nullptr,10,l),*y=fmpz_get_str(nullptr,10,u),*z=fmpz_get_str(nullptr,10,e);std::cout<<"{\"lower_mantissa\":\""<<x<<"\",\"upper_mantissa\":\""<<y<<"\",\"common_binary_exponent\":\""<<z<<"\"}";flint_free(x);flint_free(y);flint_free(z);fmpz_clear(l);fmpz_clear(u);fmpz_clear(e);}
int main(int argc,char**argv){try{
 static_assert(sizeof(arb_struct)==48&&sizeof(fmpz)==8&&sizeof(ulong)==8);
 if(std::string(flint_version)!="3.6.0"||argc!=2)throw std::runtime_error("ABI or argv");
 auto start=std::chrono::steady_clock::now();std::ifstream in(argv[1]);int count;in>>count;if(count!=21)throw std::runtime_error("fan count");
 std::array<C,21>cones;for(auto&c:cones){in>>c.det;if(c.det<=0)throw std::runtime_error("det");arb_set_si(c.weight.a,c.det);for(int i=0;i<6;i++){in>>c.q[i];if(c.q[i]<=0||c.q[i]>1000000)throw std::runtime_error("q");arb_div_ui(c.weight.a,c.weight.a,c.q[i],P);for(int k=0;k<6;k++){in>>c.r[i][k];if(std::abs(c.r[i][k])>1000000)throw std::runtime_error("ray range");rational(c.d[i][k],c.r[i][k],c.q[i]);}}}std::string extra;if(!in||in>>extra)throw std::runtime_error("input termination");
 const int64_t den=1000000007;const std::array<int64_t,6>x={1000000007,998836847,892738545,807464540,467291939,96968636};
 B a,T;rational(a,1,250);rational(T,6,125);std::map<Key,N>cache;
 auto node=[&](Key k)->N&{auto it=cache.find(k);if(it!=cache.end())return it->second;auto &n=cache[k];rational(n.b,k.first,k.second);B t,u,Da,DT,v,w;
  // Separate principal arctangents, independent of author's atan quotient.
  arb_div(t.a,n.b.a,a.a,P);arb_atan(t.a,t.a,P);arb_div(u.a,n.b.a,T.a,P);arb_atan(u.a,u.a,P);arb_sub(n.f.a,t.a,u.a,P);
  arb_mul(Da.a,n.b.a,n.b.a,P);arb_addmul(Da.a,a.a,a.a,P);arb_mul(DT.a,n.b.a,n.b.a,P);arb_addmul(DT.a,T.a,T.a,P);
  arb_div(v.a,a.a,Da.a,P);arb_div(w.a,T.a,DT.a,P);arb_sub(n.fp.a,v.a,w.a,P);
  arb_div(v.a,v.a,Da.a,P);arb_div(w.a,w.a,DT.a,P);arb_sub(v.a,w.a,v.a,P);arb_mul(n.fpp.a,v.a,n.b.a,P);arb_mul_2exp_si(n.fpp.a,n.fpp.a,1);
  return n;};
 std::array<int,6>per={0,1,2,3,4,5};std::array<B,2>total;uint64_t terms=0,zeros=0,repeats=0;
 do{std::array<int,3>pos;for(int k=0;k<6;k++){if(per[k]==1)pos[0]=k;if(per[k]==4)pos[1]=k;if(per[k]==5)pos[2]=k;}
 for(int mask=0;mask<64;mask++){std::array<int,6>sg;for(int k=0;k<6;k++)sg[k]=(mask>>k&1)?-1:1;
 for(auto&c:cones){std::array<Key,6>keys;std::array<N*,6>nd;std::array<std::array<B,6>,3>d,A;std::array<B,6>w,Cdiag,Cmix;
  for(int i=0;i<6;i++){__int128 nn=0;for(int k=0;k<6;k++)nn+=(__int128)c.r[i][k]*sg[k]*x[per[k]];__int128 dd=(__int128)den*c.q[i];if(nn<=INT64_MIN||nn>INT64_MAX||dd<=0||dd>INT64_MAX)throw std::runtime_error("rational range");int64_t n=(int64_t)nn,z=(int64_t)dd;auto g=std::gcd(n,z);keys[i]={n/g,z/g};if(!n){zeros++;throw std::runtime_error("unsupported zero base node: full reference stopped");}for(int j=0;j<i;j++)if(keys[i]==keys[j]){repeats++;throw std::runtime_error("unsupported repeated base node: full reference stopped");}nd[i]=&node(keys[i]);arb_one(w[i].a);for(int k=0;k<3;k++){arb_set(d[k][i].a,c.d[i][pos[k]].a);if(sg[pos[k]]<0)arb_neg(d[k][i].a,d[k][i].a);}}
  B diff,inv,minus,qd,qm,tmp;std::array<B,3>g;
  // A_i= sum_j (d_i-d_j)/(b_i-b_j); C_i=sum of products of these ratios.
  for(int i=0;i<6;i++)for(int j=i+1;j<6;j++){arb_sub(diff.a,nd[i]->b.a,nd[j]->b.a,P);if(arb_contains_zero(diff.a))throw std::runtime_error("cannot separate rational nodes");arb_inv(inv.a,diff.a,P);arb_mul(w[i].a,w[i].a,inv.a,P);arb_neg(minus.a,inv.a);arb_mul(w[j].a,w[j].a,minus.a,P);
   for(int k=0;k<3;k++){arb_sub(g[k].a,d[k][i].a,d[k][j].a,P);arb_mul(g[k].a,g[k].a,inv.a,P);arb_add(A[k][i].a,A[k][i].a,g[k].a,P);arb_add(A[k][j].a,A[k][j].a,g[k].a,P);}
   arb_mul(qd.a,g[0].a,g[0].a,P);arb_mul(qm.a,g[1].a,g[2].a,P);arb_add(Cdiag[i].a,Cdiag[i].a,qd.a,P);arb_add(Cdiag[j].a,Cdiag[j].a,qd.a,P);arb_add(Cmix[i].a,Cmix[i].a,qm.a,P);arb_add(Cmix[j].a,Cmix[j].a,qm.a,P);
  }
  std::array<B,2>sum;B v1,v2,v3;
  for(int i=0;i<6;i++){
   arb_mul(v1.a,d[0][i].a,d[0][i].a,P);arb_mul(v1.a,v1.a,nd[i]->fpp.a,P);
   arb_mul(v2.a,d[0][i].a,A[0][i].a,P);arb_mul(v2.a,v2.a,nd[i]->fp.a,P);arb_mul_2exp_si(v2.a,v2.a,1);arb_sub(v1.a,v1.a,v2.a,P);
   arb_mul(v3.a,A[0][i].a,A[0][i].a,P);arb_add(v3.a,v3.a,Cdiag[i].a,P);arb_addmul(v1.a,nd[i]->f.a,v3.a,P);arb_addmul(sum[0].a,w[i].a,v1.a,P);
   arb_mul(v1.a,d[1][i].a,d[2][i].a,P);arb_mul(v1.a,v1.a,nd[i]->fpp.a,P);
   arb_mul(v2.a,d[1][i].a,A[2][i].a,P);arb_addmul(v2.a,d[2][i].a,A[1][i].a,P);arb_submul(v1.a,nd[i]->fp.a,v2.a,P);
   arb_mul(v3.a,A[1][i].a,A[2][i].a,P);arb_add(v3.a,v3.a,Cmix[i].a,P);arb_addmul(v1.a,nd[i]->f.a,v3.a,P);arb_addmul(sum[1].a,w[i].a,v1.a,P);
  }
  for(int k=0;k<2;k++)arb_addmul(total[k].a,c.weight.a,sum[k].a,P);terms++;
 }}}
 while(std::next_permutation(per.begin(),per.end()));
 if(terms!=967680)throw std::runtime_error("incomplete group");B scale;arb_const_pi(scale.a,P);arb_mul_2exp_si(scale.a,scale.a,1);arb_pow_ui(scale.a,scale.a,6,P);for(auto&t:total)arb_div(t.a,t.a,scale.a,P);
 std::cout<<"{\"status\":\"COMPLETE_TWO_COMPONENT_REFERENCES\",\"precision_bits\":768,\"logical_terms_each\":"<<terms<<",\"computed_terms_each\":"<<terms<<",\"zero_base_terms\":"<<zeros<<",\"repeated_base_terms\":"<<repeats<<",\"cached_nodes\":"<<cache.size()<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"x2x2\":";output(total[0]);std::cout<<",\"x5x6\":";output(total[1]);std::cout<<",\"full_21_independent_reference\":false}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}return 0;}
