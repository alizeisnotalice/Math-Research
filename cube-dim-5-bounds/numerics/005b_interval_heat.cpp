#pragma STDC FENV_ACCESS ON
// Production-candidate repair of 005.  It makes binary64 endpoint encoding,
// floating-point preconditions, and provenance explicit.  It remains a
// conditional interval implementation rather than a formally verified library.
#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef CUBE5_BUILD_COMMAND
#error "Pass -DCUBE5_BUILD_COMMAND=... to record the exact build command"
#endif
#ifndef CUBE5_SOURCE_SHA256
#error "Pass -DCUBE5_SOURCE_SHA256=... to record the source hash"
#endif

static_assert(std::numeric_limits<double>::is_iec559);
static_assert(std::numeric_limits<double>::radix == 2);
static_assert(std::numeric_limits<double>::digits == 53);
static_assert(sizeof(double) == 8);

struct I { double lo, hi; };
static const double INF=std::numeric_limits<double>::infinity();
static void finite_or_fail(double x, const char* where) {
  if (!std::isfinite(x)) throw std::runtime_error(std::string("nonfinite ")+where);
}
static double down(double x) { finite_or_fail(x,"input"); double y=std::nextafter(x,-INF); finite_or_fail(y,"nextafter-down"); return y; }
static double up(double x) { finite_or_fail(x,"input"); double y=std::nextafter(x,INF); finite_or_fail(y,"nextafter-up"); return y; }
static I exact(double x) { finite_or_fail(x,"exact"); return {x,x}; }
static I input_ratio(int n,int d) { double x=double(n)/double(d); return {down(x),up(x)}; }
static I add(I a,I b) { return {down(a.lo+b.lo),up(a.hi+b.hi)}; }
static I neg(I a) { return {-a.hi,-a.lo}; }
static I sub(I a,I b) { return add(a,neg(b)); }
static I mul(I a,I b) {
  double x[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
  for(double v:x) finite_or_fail(v,"multiply");
  return {down(*std::min_element(x,x+4)),up(*std::max_element(x,x+4))};
}
static I reciprocal_positive(I a) {
  if(!(a.lo>0) || !std::isfinite(a.hi)) throw std::runtime_error("nonpositive denominator");
  return {down(1.0/a.hi),up(1.0/a.lo)};
}
static I divide(I a,I positive) { return mul(a,reciprocal_positive(positive)); }
struct C { I re,im; };
static C cmul(C a,C b) { return {sub(mul(a.re,b.re),mul(a.im,b.im)),add(mul(a.re,b.im),mul(a.im,b.re))}; }
static C cinv(C a) { I den=add(mul(a.re,a.re),mul(a.im,a.im)); I inv=reciprocal_positive(den); return {mul(a.re,inv),mul(neg(a.im),inv)}; }

struct Term { std::array<std::array<int,5>,5> ray; std::array<int,5> rate; int det; };
static const int base[7][5][5] = {
 {{1,1,0,0,0},{1,0,0,0,0},{4,1,1,1,1},{3,1,1,1,0},{2,1,1,0,0}},
 {{1,1,0,0,0},{4,1,1,1,1},{2,1,1,1,1},{3,1,1,1,0},{2,1,1,0,0}},
 {{1,1,1,1,0},{1,1,0,0,0},{2,1,1,1,1},{3,1,1,1,0},{2,1,1,0,0}},
 {{1,1,1,1,0},{2,2,2,1,1},{1,1,0,0,0},{2,1,1,1,1},{2,1,1,0,0}},
 {{1,1,1,1,0},{2,2,2,1,1},{1,1,0,0,0},{3,3,2,2,2},{2,1,1,1,1}},
 {{1,1,1,1,1},{1,1,1,1,0},{2,2,2,1,1},{3,3,2,2,2},{2,1,1,1,1}},
 {{1,1,1,1,0},{1,1,1,0,0},{2,2,2,1,1},{1,1,0,0,0},{2,1,1,0,0}}
};
static const int dets[7]={1,2,2,2,2,1,1};
static const int rates[7][5]={{16,16,64,48,32},{16,64,36,48,32},{24,16,36,48,32},{24,48,16,36,32},{24,48,16,72,36},{30,24,48,72,36},{24,24,48,16,32}};
static std::vector<Term> make_terms() {
  std::vector<Term> out; out.reserve(26880); std::array<int,5> p={0,1,2,3,4};
  do { for(int mask=0;mask<32;++mask) for(int c=0;c<7;++c) {
    Term z; z.det=dets[c];
    for(int j=0;j<5;++j) { z.rate[j]=rates[c][j]; for(int row=0;row<5;++row) z.ray[j][row]=((mask>>row)&1?-1:1)*base[c][j][p[row]]; }
    out.push_back(z);
  }} while(std::next_permutation(p.begin(),p.end()));
  return out;
}

static I heat(const std::array<I,5>& x,double time,const std::vector<Term>& terms) {
  I sum=exact(0);
  for(const auto& term:terms) {
    C product={exact(double(term.det)),exact(0)};
    for(int j=0;j<5;++j) {
      I dot=exact(0); for(int k=0;k<5;++k) dot=add(dot,mul(exact(double(term.ray[j][k])),x[k]));
      C denominator={I{down(time*term.rate[j]),up(time*term.rate[j])},neg(dot)};
      product=cmul(product,cinv(denominator));
    }
    sum=add(sum,product.re);
  }
  I pi={down(0x1.921fb54442d18p+1),up(0x1.921fb54442d18p+1)};
  I norm=exact(1), two_pi=mul(exact(2),pi); for(int k=0;k<5;++k) norm=mul(norm,two_pi);
  return divide(sum,norm);
}

static std::string hex64(double x) {
  std::ostringstream s; s<<"0x"<<std::hex<<std::setw(16)<<std::setfill('0')<<std::bit_cast<std::uint64_t>(x); return s.str();
}
static void write_interval(std::ostream& f,I x) { f<<"{\"lo_bits\":\""<<hex64(x.lo)<<"\",\"hi_bits\":\""<<hex64(x.hi)<<"\"}"; }
static std::string escape_json(std::string s) { std::string r; for(char c:s) { if(c=='\\' || c=='\"') r+='\\'; r+=c; } return r; }
static std::string fp_control() {
#if defined(__aarch64__)
  std::uint64_t fpcr; asm volatile("mrs %0, fpcr" : "=r"(fpcr));
  if(fpcr & ((std::uint64_t(1)<<24)|(std::uint64_t(1)<<19))) throw std::runtime_error("ARM FPCR FTZ/FZ16 enabled");
  std::ostringstream s; s<<"aarch64_fpcr_0x"<<std::hex<<fpcr; return s.str();
#elif defined(__x86_64__) || defined(_M_X64)
  unsigned mxcsr; asm volatile("stmxcsr %0" : "=m"(mxcsr));
  if(mxcsr & ((1u<<15)|(1u<<6))) throw std::runtime_error("x86 MXCSR FTZ/DAZ enabled");
  std::ostringstream s; s<<"x86_mxcsr_0x"<<std::hex<<mxcsr; return s.str();
#else
#error "005b only supports an explicitly checked ARM64 or x86-64 FP control register"
#endif
}
static const char* architecture() {
#if defined(__aarch64__)
  return "aarch64";
#elif defined(__x86_64__)
  return "x86_64";
#else
  return "unknown";
#endif
}

int main(int argc,char** argv) {
 try {
  if(fegetround()!=FE_TONEAREST) throw std::runtime_error("rounding mode is not FE_TONEAREST");
  std::string control=fp_control();
  int count=argc>1?std::stoi(argv[1]):8; if(count<8||count>100) throw std::runtime_error("count must be 8..100");
  auto terms=make_terms(); if(terms.size()!=26880) throw std::runtime_error("term count");
  const int tails[8][4]={{936,673,68,917},{201,693,737,586},{901,455,74,837},{926,419,373,965},{148,785,707,948},{139,53,348,452},{399,641,776,78},{931,367,597,638}};
  std::vector<std::array<I,5>> points; points.reserve(count);
  for(auto& tail:tails) points.push_back({exact(1),input_ratio(tail[0],1009),input_ratio(tail[1],1009),input_ratio(tail[2],1009),input_ratio(tail[3],1009)});
  std::uint64_t state=0x9e3779b97f4a7c15ULL;
  while((int)points.size()<count) { std::array<I,5>x={exact(1),exact(0),exact(0),exact(0),exact(0)}; for(int j=1;j<5;++j) { state=state*6364136223846793005ULL+1442695040888963407ULL; int n=int((state>>11)%1000003ULL)+1; x[j]=input_ratio(n,1000004); } points.push_back(x); }
  const double times[5]={1.0/128,1.0/32,1.0/8,1.0/2,1}; std::vector<std::array<I,5>> values(points.size());
  auto start=std::chrono::steady_clock::now(); for(size_t i=0;i<points.size();++i) for(int t=0;t<5;++t) values[i][t]=heat(points[i],times[t],terms);
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::ofstream f("005b_interval_results.json");
  f<<"{\n  \"endpoint_encoding\":\"IEEE754-binary64-bit-pattern\",\n"
   <<"  \"assumptions\":\"IEC60559 binary64 scalar operations; FE_TONEAREST; FTZ/DAZ disabled; no fast-math or FMA contraction; every primitive expanded by nextafter\",\n"
   <<"  \"source_sha256\":\""<<CUBE5_SOURCE_SHA256<<"\",\n  \"build_command\":\""<<escape_json(CUBE5_BUILD_COMMAND)<<"\",\n"
   <<"  \"compiler\":\""<<escape_json(__VERSION__)<<"\",\n  \"architecture\":\""<<architecture()<<"\",\n  \"fp_control\":\""<<control<<"\",\n"
   <<"  \"terms\":"<<terms.size()<<",\"points\":"<<points.size()<<",\"seconds\":"<<std::setprecision(17)<<seconds<<",\n  \"times_bits\":[";
  for(int t=0;t<5;++t) f<<(t?",":"")<<"\""<<hex64(times[t])<<"\""; f<<"],\n  \"point_results\":[\n";
  for(size_t i=0;i<points.size();++i) {
    f<<"    {\"index\":"<<i<<",\"locked_002G\":"<<(i<8?"true":"false")<<",\"input_intervals\":[";
    for(int j=0;j<5;++j) { if(j) f<<","; write_interval(f,points[i][j]); }
    f<<"],\"values\":["; for(int t=0;t<5;++t) { if(t) f<<","; write_interval(f,values[i][t]); } f<<"]}"<<(i+1==points.size()?"\n":" ,\n");
  }
  f<<"  ]\n}\n";
  std::cout<<"terms="<<terms.size()<<" points="<<points.size()<<" seconds="<<seconds<<" fp="<<control<<"\n";
 } catch(const std::exception& e) { std::cerr<<"FAIL "<<e.what()<<"\n"; return 1; }
}
