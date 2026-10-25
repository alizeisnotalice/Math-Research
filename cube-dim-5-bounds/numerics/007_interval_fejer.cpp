// Reads only Arb-generated binary64 enclosures; no libm trig is used here.
#define CUBE5_BUILD_COMMAND "embedded-by-007"
#define CUBE5_SOURCE_SHA256 "see-005b-locked-source"
#define main embedded_005b_main
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#include "005b_interval_heat.cpp"
#pragma clang diagnostic pop
#undef main

#ifndef CUBE7_SOURCE_SHA256
#error "Pass CUBE7_SOURCE_SHA256"
#endif
#ifndef CUBE7_INPUT_SHA256
#error "Pass CUBE7_INPUT_SHA256"
#endif
#ifndef CUBE7_BUILD_COMMAND
#error "Pass CUBE7_BUILD_COMMAND"
#endif

struct Node { int segment, j; I time, weight; };
static double from_bits(const std::string& s) {
  std::uint64_t u=std::stoull(s,nullptr,0); return std::bit_cast<double>(u);
}
static std::vector<Node> load_nodes() {
  std::ifstream f("007_fejer_nodes.bits"); if(!f) throw std::runtime_error("cannot open 007_fejer_nodes.bits");
  std::vector<Node> nodes; int segment,j; std::string a,b,c,d;
  while(f>>segment>>j>>a>>b>>c>>d) {
    I time={from_bits(a),from_bits(b)}, weight={from_bits(c),from_bits(d)};
    if(!(std::isfinite(time.lo)&&std::isfinite(time.hi)&&std::isfinite(weight.lo)&&std::isfinite(weight.hi)&&time.lo<=time.hi&&weight.lo>0&&weight.lo<=weight.hi)) throw std::runtime_error("bad node interval");
    nodes.push_back({segment,j,time,weight});
  }
  if(nodes.size()!=175) throw std::runtime_error("expected 175 nodes");
  return nodes;
}
static I heat_interval(const std::array<I,5>& x,I time,const std::vector<Term>& terms) {
  I sum=exact(0);
  for(const auto& term:terms) {
    C product={exact(double(term.det)),exact(0)};
    for(int j=0;j<5;++j) {
      I dot=exact(0); for(int k=0;k<5;++k) dot=add(dot,mul(exact(double(term.ray[j][k])),x[k]));
      C denominator={mul(time,exact(double(term.rate[j]))),neg(dot)};
      product=cmul(product,cinv(denominator));
    }
    sum=add(sum,product.re);
  }
  I pi={down(0x1.921fb54442d18p+1),up(0x1.921fb54442d18p+1)};
  I norm=exact(1),two_pi=mul(exact(2),pi); for(int k=0;k<5;++k) norm=mul(norm,two_pi);
  return divide(sum,norm);
}
int main() {
 try {
  if(fegetround()!=FE_TONEAREST) throw std::runtime_error("rounding mode is not FE_TONEAREST");
  std::string control=fp_control(); auto nodes=load_nodes(); auto terms=make_terms();
  const int tail[4]={936,673,68,917};
  std::array<I,5> x={exact(1),input_ratio(tail[0],1009),input_ratio(tail[1],1009),input_ratio(tail[2],1009),input_ratio(tail[3],1009)};
  auto started=std::chrono::steady_clock::now(); I sum=exact(0);
  for(const auto& node:nodes) sum=add(sum,mul(node.weight,heat_interval(x,node.time,terms)));
  double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  std::ofstream f("007_heat_quadrature.json");
  f<<"{\n  \"endpoint_encoding\":\"IEEE754-binary64-bit-pattern\",\n"
   <<"  \"source_sha256\":\""<<CUBE7_SOURCE_SHA256<<"\",\n  \"input_nodes_sha256\":\""<<CUBE7_INPUT_SHA256<<"\",\n"
   <<"  \"build_command\":\""<<escape_json(CUBE7_BUILD_COMMAND)<<"\",\n"
   <<"  \"compiler\":\""<<escape_json(__VERSION__)<<"\",\n  \"architecture\":\""<<architecture()<<"\",\n  \"fp_control\":\""<<control<<"\",\n"
   <<"  \"nodes\":175,\"terms\":"<<terms.size()<<",\"seconds\":"<<std::setprecision(17)<<seconds<<",\n  \"quadrature\":";
  write_interval(f,sum); f<<",\n  \"scope\":\"Interval heat sum using Arb-enclosed Fejer input bit patterns; analytic interpolation and tail errors are added only in 007_final_verify.py.\"\n}\n";
  std::cout<<"nodes=175 seconds="<<seconds<<" fp="<<control<<"\n";
 } catch(const std::exception& e) { std::cerr<<"FAIL "<<e.what()<<"\n"; return 1; }
}
