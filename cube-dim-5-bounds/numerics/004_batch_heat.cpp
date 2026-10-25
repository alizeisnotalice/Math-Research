#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

struct Term { std::array<std::array<int,5>,5> ray; std::array<double,5> rate; double det; };
struct Eval { double value; double absolute_sum; };

static const int base[7][5][5] = {
 {{1,1,0,0,0},{1,0,0,0,0},{4,1,1,1,1},{3,1,1,1,0},{2,1,1,0,0}},
 {{1,1,0,0,0},{4,1,1,1,1},{2,1,1,1,1},{3,1,1,1,0},{2,1,1,0,0}},
 {{1,1,1,1,0},{1,1,0,0,0},{2,1,1,1,1},{3,1,1,1,0},{2,1,1,0,0}},
 {{1,1,1,1,0},{2,2,2,1,1},{1,1,0,0,0},{2,1,1,1,1},{2,1,1,0,0}},
 {{1,1,1,1,0},{2,2,2,1,1},{1,1,0,0,0},{3,3,2,2,2},{2,1,1,1,1}},
 {{1,1,1,1,1},{1,1,1,1,0},{2,2,2,1,1},{3,3,2,2,2},{2,1,1,1,1}},
 {{1,1,1,1,0},{1,1,1,0,0},{2,2,2,1,1},{1,1,0,0,0},{2,1,1,0,0}}
};
static const int dets[7] = {1,2,2,2,2,1,1};
static const int rates[7][5] = {
 {16,16,64,48,32},{16,64,36,48,32},{24,16,36,48,32},{24,48,16,36,32},
 {24,48,16,72,36},{30,24,48,72,36},{24,24,48,16,32}
};

std::vector<Term> make_terms() {
  std::vector<Term> out; out.reserve(26880);
  std::array<int,5> p={0,1,2,3,4};
  do {
    for(int mask=0; mask<32; ++mask) {
      for(int c=0;c<7;++c) {
        Term term; term.det=dets[c];
        for(int j=0;j<5;++j) { term.rate[j]=rates[c][j]; for(int row=0;row<5;++row) term.ray[j][row]=((mask>>row)&1 ? -1:1)*base[c][j][p[row]]; }
        out.push_back(term);
      }
    }
  } while(std::next_permutation(p.begin(),p.end()));
  return out;
}

Eval heat(const std::array<double,5>& x, double time, const std::vector<Term>& terms) {
  std::complex<double> sum(0,0), compensation(0,0); double absolute_sum=0;
  for(const auto& term: terms) {
    std::complex<double> product(term.det,0);
    for(int j=0;j<5;++j) { double dot=0; for(int k=0;k<5;++k) dot += term.ray[j][k]*x[k]; product /= std::complex<double>(time*term.rate[j],-dot); }
    absolute_sum += std::abs(product);
    auto y=product-compensation; auto next=sum+y; compensation=(next-sum)-y; sum=next;
  }
  constexpr double pi=3.141592653589793238462643383279502884;
  const double norm=std::pow(2*pi,5);
  return {sum.real()/norm, absolute_sum/norm};
}

int main(int argc,char** argv) {
  const int count=argc>1 ? std::stoi(argv[1]) : 1000;
  const auto terms=make_terms();
  if(terms.size()!=26880) return 2;
  std::vector<std::array<double,5>> points={
    {1,936.0/1009,673.0/1009,68.0/1009,917.0/1009}, {1,201.0/1009,693.0/1009,737.0/1009,586.0/1009},
    {1,901.0/1009,455.0/1009,74.0/1009,837.0/1009}, {1,926.0/1009,419.0/1009,373.0/1009,965.0/1009},
    {1,148.0/1009,785.0/1009,707.0/1009,948.0/1009}, {1,139.0/1009,53.0/1009,348.0/1009,452.0/1009},
    {1,399.0/1009,641.0/1009,776.0/1009,78.0/1009}, {1,931.0/1009,367.0/1009,597.0/1009,638.0/1009}
  };
  uint64_t state=0x9e3779b97f4a7c15ULL;
  while((int)points.size()<count) { std::array<double,5> x={1,0,0,0,0}; for(int j=1;j<5;++j) { state=state*6364136223846793005ULL+1442695040888963407ULL; x[j]=((state>>11)&((1ULL<<53)-1))/double(1ULL<<53); } points.push_back(x); }
  const std::array<double,5> times={1.0/128,1.0/32,1.0/8,1.0/2,1};
  std::vector<std::vector<Eval>> values(times.size(),std::vector<Eval>(points.size()));
  const auto start=std::chrono::steady_clock::now();
  for(size_t t=0;t<times.size();++t) for(size_t i=0;i<points.size();++i) values[t][i]=heat(points[i],times[t],terms);
  const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  std::ofstream f("004_cpp_results.json"); f<<std::setprecision(17);
  f<<"{\n  \"terms\": "<<terms.size()<<",\n  \"points\": "<<points.size()<<",\n  \"seconds\": "<<seconds<<",\n  \"times\": [";
  for(size_t t=0;t<times.size();++t) f<<(t?", ":"")<<times[t]; f<<"],\n  \"locked_values\": [\n";
  for(int i=0;i<8;++i) { f<<"    {\"point\": ["; for(int k=0;k<5;++k) f<<(k?", ":"")<<points[i][k]; f<<"], \"values\": ["; for(size_t t=0;t<times.size();++t) f<<(t?", ":"")<<"{\"p\": "<<values[t][i].value<<", \"abs_sum\": "<<values[t][i].absolute_sum<<"}"; f<<"]}"<<(i==7?"\n":" ,\n"); }
  f<<"  ]\n}\n";
  std::cout<<"terms="<<terms.size()<<" points="<<points.size()<<" evaluations="<<points.size()*times.size()<<" seconds="<<seconds<<"\n";
}
