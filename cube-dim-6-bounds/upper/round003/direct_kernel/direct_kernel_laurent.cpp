#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>
struct Ray{long long q;std::array<long long,6>r;};struct Cone{long long det;std::array<Ray,6>rays;};struct Node{long long n,d;};
bool eq(const Node&a,const Node&b){return (__int128)a.n*b.d==(__int128)b.n*a.d;}bool lt(const Node&a,const Node&b){return (__int128)a.n*b.d<(__int128)b.n*a.d;}long double val(const Node&a){return (long double)a.n/a.d;}
struct V{long double p=0,e1=0,e3=0,e5=0;};V sub(V a,V b){return {a.p-b.p,a.e1-b.e1,a.e3-b.e3,a.e5-b.e5};}V divv(V a,long double d){return {a.p/d,a.e1/d,a.e3/d,a.e5/d};}V jet(const Node&x,int ord){if(x.n)return ord?V{}:V{(x.n>0?0.5L:-0.5L),0,0,0};if(ord==1)return{0,1,0,0};if(ord==3)return{0,0,-1.L/3,0};if(ord==5)return{0,0,0,1.L/5};return{};}
struct K{long double s=0,c=0,a=0;void add(long double x){long double y=x-c,t=s+y;c=(t-s)-y;s=t;a+=std::fabs(x);}};
V dd(std::array<Node,6>x){std::sort(x.begin(),x.end(),lt);V t[6][6]{};for(int i=0;i<6;i++)t[i][0]=jet(x[i],0);for(int o=1;o<6;o++)for(int i=o;i<6;i++)t[i][o]=eq(x[i],x[i-o])?jet(x[i],o):divv(sub(t[i][o-1],t[i-1][o-1]),val(x[i])-val(x[i-o]));return t[5][5];}
int main(int argc,char**argv){int wanted=argc>1?std::stoi(argv[1]):4;std::string path=argc>2?argv[2]:"face_points.txt";std::ifstream fi("fan_integer_input.txt"),pi(path);int cn;fi>>cn;std::vector<Cone>cs(cn);for(auto&c:cs){fi>>c.det;for(auto&r:c.rays){fi>>r.q;for(auto&z:r.r)fi>>z;}}long long D;int pn;pi>>D>>pn;std::vector<std::array<long long,6>>xs(pn);for(auto&x:xs)for(auto&z:x)pi>>z;wanted=std::min(wanted,pn);std::cout<<"{\n\"status\":\"laurent_long_double_diagnostic_only\",\n\"points\":[\n";
  for(int p=0;p<wanted;p++){auto begin=std::chrono::steady_clock::now();K P,E1,E3,E5;long long zeros=0,repeats=0;std::array<int,6>perm={0,1,2,3,4,5};do{for(int mask=0;mask<64;mask++){std::array<long long,6>y{};for(int i=0;i<6;i++)y[i]=((mask>>i)&1?-1:1)*xs[p][perm[i]];for(const auto&c:cs){std::array<Node,6>b{};for(int j=0;j<6;j++){long long n=0;for(int i=0;i<6;i++)n+=y[i]*c.rays[j].r[i];b[j]={n,D*c.rays[j].q};if(!n)zeros++;}auto z=b;std::sort(z.begin(),z.end(),lt);for(int j=1;j<6;j++)if(eq(z[j-1],z[j]))repeats++;long double a=c.det;for(const auto&r:c.rays)a/=r.q;V q=dd(b);P.add(a*q.p);E1.add(a*q.e1);E3.add(a*q.e3);E5.add(a*q.e5);}}}while(std::next_permutation(perm.begin(),perm.end()));auto ratio=[](K z){return z.s==0?0.L:z.a/std::fabs(z.s);};long double elapsed=std::chrono::duration<long double>(std::chrono::steady_clock::now()-begin).count();std::cout<<"{\"index\":"<<p<<",\"zero_nodes\":"<<zeros<<",\"repeated_pairs\":"<<repeats<<",\"pi_coeff\":\""<<std::setprecision(21)<<P.s<<"\",\"eps_m1\":\""<<E1.s<<"\",\"eps_m3\":\""<<E3.s<<"\",\"eps_m5\":\""<<E5.s<<"\",\"cancellation_ratio_pi\":\""<<ratio(P)<<"\",\"cancellation_ratio_m1\":\""<<ratio(E1)<<"\",\"cancellation_ratio_m3\":\""<<ratio(E3)<<"\",\"cancellation_ratio_m5\":\""<<ratio(E5)<<"\",\"G_over_pi5_long_double\":\""<<P.s/64<<"\",\"elapsed_seconds\":\""<<elapsed<<"\"}"<<(p+1==wanted?"\n":",\n");}std::cout<<"]\n}\n";}
