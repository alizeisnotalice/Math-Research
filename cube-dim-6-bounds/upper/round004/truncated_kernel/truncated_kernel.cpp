#include <algorithm>
#include <array>
#include <chrono>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

// All node ordering/equality is exact rational arithmetic.  On this host
// long double is binary64 (verified separately), so values are exploratory.
struct Ray { long long q; std::array<long long, 6> r; };
struct Cone { long long det; std::array<Ray, 6> rays; };
struct Node { long long n, d; };
bool eq(const Node& a, const Node& b) { return (__int128)a.n * b.d == (__int128)b.n * a.d; }
bool lt(const Node& a, const Node& b) { return (__int128)a.n * b.d < (__int128)b.n * a.d; }
long double value(const Node& a) { return (long double)a.n / a.d; }
struct Kahan { long double s = 0, c = 0; void add(long double x) { long double y=x-c, t=s+y; c=(t-s)-y; s=t; } };

// F^(k)(x)/k! for F(x)=atan(x/e).  k>=1 follows by differentiating the two
// logarithms defining atan; it is used only for an exactly equal-node block.
long double jet(long double x, long double e, int k) {
  if (k == 0) return std::atan(x / e);
  using C = std::complex<long double>;
  const C ii(0, 1), a(e, x), b(e, -x);
  C numerator = std::pow(ii, k) / std::pow(a, k) - std::pow(-ii, k) / std::pow(b, k);
  C result = ((k & 1) ? 1.L : -1.L) * numerator / ((long double)(2 * k) * ii);
  return result.real();
}

long double hermite_dd(std::array<Node, 6> nodes, long double e) {
  std::sort(nodes.begin(), nodes.end(), lt);
  long double table[6][6]{};
  for (int i=0; i<6; ++i) table[i][0] = jet(value(nodes[i]), e, 0);
  for (int order=1; order<6; ++order) for (int i=order; i<6; ++i) {
    if (eq(nodes[i], nodes[i-order])) table[i][order] = jet(value(nodes[i]), e, order);
    else table[i][order] = (table[i][order-1]-table[i-1][order-1]) / (value(nodes[i])-value(nodes[i-order]));
  }
  return table[5][5];
}

int main(int argc, char** argv) {
  const std::string fan_path = argc > 1 ? argv[1] : "fan_integer_input.txt";
  const std::string point_path = argc > 2 ? argv[2] : "points.txt";
  const std::string eps_path = argc > 3 ? argv[3] : "epsilons.txt";
  std::ifstream fan(fan_path), points(point_path), eps_file(eps_path);
  int cones_n; fan >> cones_n; std::vector<Cone> cones(cones_n);
  for (auto& c: cones) { fan >> c.det; for (auto& r: c.rays) { fan >> r.q; for (auto& z:r.r) fan >> z; } }
  long long den; int points_n; points >> den >> points_n; std::vector<std::array<long long,6>> xs(points_n);
  for (auto& x:xs) for (auto& z:x) points >> z;
  int eps_n; eps_file >> eps_n; std::vector<std::pair<long long,long long>> eps(eps_n);
  for (auto& x:eps) eps_file >> x.first >> x.second;
  std::cout << "{\n  \"status\":\"exploratory_double_H_epsilon\",\n  \"LDBL_MANT_DIG\":" << LDBL_MANT_DIG << ",\n  \"sizeof_long_double\":" << sizeof(long double) << ",\n  \"evaluations\":[\n";
  bool first = true;
  for (int p=0; p<points_n; ++p) for (auto [epn, epd]:eps) {
    auto started = std::chrono::steady_clock::now(); Kahan sum, abs_sum; std::array<int,6> perm={0,1,2,3,4,5}; long long terms=0, equal_events=0;
    do { for (int mask=0; mask<64; ++mask) {
      std::array<long long,6> tx{}; for (int i=0;i<6;++i) tx[i]=((mask>>i)&1 ? -1 : 1)*xs[p][perm[i]];
      for (const auto& c:cones) {
        std::array<Node,6> b{}; long double alpha=c.det;
        for (int j=0;j<6;++j) { long long n=0; for(int i=0;i<6;++i) n+=tx[i]*c.rays[j].r[i]; b[j]={n, den*c.rays[j].q}; alpha/=c.rays[j].q; }
        auto sorted=b; std::sort(sorted.begin(), sorted.end(), lt); for(int j=1;j<6;++j) if(eq(sorted[j-1],sorted[j])) ++equal_events;
        long double term=alpha*hermite_dd(b, (long double)epn/epd); sum.add(term); abs_sum.add(std::fabs(term)); ++terms;
      }
    }} while(std::next_permutation(perm.begin(),perm.end()));
    long double h=sum.s/std::pow(2.L*std::numbers::pi_v<long double>,6);
    long double seconds=std::chrono::duration<long double>(std::chrono::steady_clock::now()-started).count();
    if (!first) std::cout << ",\n"; first=false;
    std::cout << "    {\"point_index\":" << p << ",\"numerators\":["; for(int i=0;i<6;++i){if(i)std::cout<<",";std::cout<<xs[p][i];}
    std::cout << "],\"denominator\":" << den << ",\"epsilon\":\""<<epn<<"/"<<epd<<"\",\"raw_terms\":"<<terms<<",\"equal_adjacent_events\":"<<equal_events<<",\"signed_sum_binary64\":\""<<std::setprecision(21)<<sum.s<<"\",\"absolute_sum_binary64\":\""<<abs_sum.s<<"\",\"cancellation_ratio\":\""<<(abs_sum.s/std::fabs(sum.s))<<"\",\"H_epsilon_binary64\":\""<<h<<"\",\"elapsed_seconds\":\""<<seconds<<"\"}";
  }
  std::cout << "\n  ]\n}\n";
}
