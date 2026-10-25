// Exact orbit-weighted corner pilot. Existing graph cache is never expanded.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "round012_graph_moments_gmp.hpp"
using I=round012::Integer;
using U=round012::Code;
using Poly=std::unordered_map<U,I>;
using Term=std::pair<U,I>;
constexpr const char* BACKEND="7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8";
std::unique_ptr<round012::Engine> engine;
std::unordered_set<U> available;
int edge[5][5];
std::vector<Term> f;
std::string source_sha;

uint64_t checksum(const std::string&s){uint64_t h=14695981039346656037ULL;for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;}
std::string hex(uint64_t x){std::ostringstream s;s<<std::hex<<std::setw(16)<<std::setfill('0')<<x;return s.str();}
unsigned degree(U code){return round012::degree(round012::decode(code));}
unsigned parity(U code){return round012::parity(round012::decode(code));}
U add_code(U a,U b){auto x=round012::decode(a),y=round012::decode(b);for(unsigned e=0;e<10;e++){x[e]+=y[e];if(x[e]>31)throw std::runtime_error("graph edge exponent exceeds 31");}if(round012::degree(x)>26)throw std::runtime_error("graph degree exceeds 26");return round012::encode(x);}
I B(U code){if(parity(code))return 0;U key=engine->canonical(code);if(!available.count(key))throw std::runtime_error("missing cached graph "+std::to_string(key)+" degree "+std::to_string(degree(key)));I value=engine->B(key);assert(engine->stats().computed==0);return value;}
void prune(Poly&p){for(auto it=p.begin();it!=p.end();)if(it->second==0)it=p.erase(it);else++it;}
int sign(const std::vector<int>&p){int n=0;for(size_t i=0;i<p.size();i++)for(size_t j=i+1;j<p.size();j++)n+=p[i]>p[j];return n%2?-1:1;}
Poly minor(const std::vector<int>&rows,const std::vector<int>&cols){Poly out;std::vector<int> p(rows.size());for(size_t i=0;i<p.size();i++)p[i]=i;do{round012::Counts k{};for(size_t i=0;i<p.size();i++)if(rows[i]!=cols[p[i]])++k[edge[rows[i]][cols[p[i]]]];out[round012::encode(k)]+=sign(p);}while(std::next_permutation(p.begin(),p.end()));prune(out);return out;}
Poly multiply_raw(const Poly&a,const std::vector<Term>&b){Poly out;for(const auto&x:a)for(const auto&y:b)out[add_code(x.first,y.first)]+=x.second*y.second;prune(out);return out;}
Poly multiply_orbit(const Poly&a){Poly out;out.reserve(a.size()*2+100);for(const auto&x:a)for(const auto&y:f)out[engine->canonical(add_code(x.first,y.first))]+=x.second*y.second;prune(out);return out;}
void save_power(const std::string&path,const Poly&p,int n){
 std::vector<U> keys;keys.reserve(p.size());for(const auto&t:p)keys.push_back(t.first);std::sort(keys.begin(),keys.end());
 std::string tmp=path+".partial";std::ofstream o(tmp);o<<"astra010-f-power-v1 "<<source_sha<<" "<<n<<"\n";
 for(U key:keys){std::string line=std::to_string(key)+"\t"+p.at(key).get_str();o<<line<<"\t"<<hex(checksum(line))<<"\n";}o.close();if(std::rename(tmp.c_str(),path.c_str()))throw std::runtime_error("checkpoint rename failed");
}
std::pair<Poly,int> load_power(const std::string&path){
 std::ifstream in(path);if(!in)return {Poly{{0,I(1)}},0};std::string line;std::getline(in,line);std::istringstream head(line);std::string protocol,sha;int n=-1;head>>protocol>>sha>>n;if(protocol!="astra010-f-power-v1"||sha!=source_sha||n<0||n>11)throw std::runtime_error("checkpoint fingerprint mismatch");
 Poly p;while(std::getline(in,line)){auto a=line.find('\t'),b=line.find('\t',a+1);if(a==std::string::npos||b==std::string::npos||line.substr(b+1)!=hex(checksum(line.substr(0,b))))throw std::runtime_error("power checkpoint checksum failure");U code=std::stoull(line.substr(0,a));I coefficient(line.substr(a+1,b-a-1));if(degree(code)>unsigned(2*n)||engine->canonical(code)!=code||coefficient==0||p.count(code))throw std::runtime_error("invalid power checkpoint row");p[code]=coefficient;}return {std::move(p),n};
}
mpq_class expectation_product(const Poly&orbit,const Poly&factor,int power){
 std::array<std::vector<Term>,32> buckets;for(const auto&t:factor)buckets[parity(t.first)].push_back(t);
 std::array<I,27> sums{};size_t terms=0;
 for(const auto&a:orbit)for(const auto&b:buckets[parity(a.first)]){U code=add_code(a.first,b.first);sums[degree(code)]+=a.second*b.second*B(code);++terms;}
 mpq_class result=0;for(unsigned n=0;n<=26;n++){mpq_class term(sums[n],engine->denominator(n));term.canonicalize();result+=term;}
 I scale;mpz_ui_pow_ui(scale.get_mpz_t(),880,power);result/=scale;result.canonicalize();std::cout<<" expectation factor terms "<<factor.size()<<", evaluated products "<<terms<<std::endl;return result;
}
int main(int argc,char**argv){try{
 if(argc!=5){std::cerr<<"usage: corner result.json cache.tsv power_checkpoint.tsv source_sha256\n";return 2;}source_sha=argv[4];if(source_sha.size()!=64)throw std::runtime_error("source SHA256 required");auto began=std::chrono::steady_clock::now();
 int e=0;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)edge[i][j]=edge[j][i]=e++;
 engine=std::make_unique<round012::Engine>(argv[2],BACKEND,0.0,1,1000000);
 available.insert(0);std::ifstream cache(argv[2]);std::string line;std::getline(cache,line);while(std::getline(cache,line))available.insert(std::stoull(line.substr(0,line.find('\t'))));assert(available.size()==engine->cached_moments());
 f.emplace_back(0,I(880));for(unsigned i=0;i<10;i++){f.emplace_back(U(1)<<(5*i),I(88));f.emplace_back(U(2)<<(5*i),I(27));for(unsigned j=i+1;j<10;j++){bool adjacent=round012::eu[i]==round012::eu[j]||round012::eu[i]==round012::ev[j]||round012::ev[i]==round012::eu[j]||round012::ev[i]==round012::ev[j];f.emplace_back((U(1)<<(5*i))+(U(1)<<(5*j)),I(adjacent?-22:22));}}
 assert(f.size()==66);Poly raw{{0,I(1)}},compressed{{0,I(1)}};
 for(int n=1;n<=3;n++){raw=multiply_raw(raw,f);compressed=multiply_orbit(compressed);Poly expected;for(const auto&t:raw)expected[engine->canonical(t.first)]+=t.second;prune(expected);assert(expected==compressed);std::cout<<"independent raw/orbit crosscheck power "<<n<<", raw "<<raw.size()<<", orbit "<<compressed.size()<<std::endl;}
 auto saved=load_power(argv[3]);Poly power=std::move(saved.first);int n=saved.second;std::string seventh=std::string(argv[3])+".f7";
 std::cout<<"resume power "<<n<<", orbit terms "<<power.size()<<", available graph moments "<<available.size()<<std::endl;
 for(;n<11;){power=multiply_orbit(power);++n;save_power(argv[3],power,n);if(n==7)save_power(seventh,power,n);double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();std::cout<<"saved power "<<n<<", orbit terms "<<power.size()<<", seconds "<<seconds<<std::endl;}
 auto stored7=load_power(seventh);if(stored7.second!=7)throw std::runtime_error("completed f7 checkpoint missing");
 Poly determinant=minor({0,1,2,3,4},{0,1,2,3,4}),squared,cofactor;
 for(const auto&a:determinant)for(const auto&b:determinant)squared[add_code(a.first,b.first)]+=a.second*b.second;prune(squared);
 for(int i=0;i<5;i++)for(int j=i;j<5;j++){std::vector<int> rows,cols;for(int k=0;k<5;k++){if(k!=j)rows.push_back(k);if(k!=i)cols.push_back(k);}int multiplier=((i+j)%2?-1:1)*(i==j?1:2);for(const auto&t:minor(rows,cols))cofactor[t.first]+=multiplier*t.second;}prune(cofactor);
 mpq_class Edf7=expectation_product(stored7.first,determinant,7);mpq_class Ed2f7=expectation_product(stored7.first,squared,7);mpq_class Ef7=expectation_product(stored7.first,Poly{{0,I(1)}},7);mpq_class Ebf11=expectation_product(power,cofactor,11);
 // First moment of f is also available by the closed form 1+27/616.
 Poly f_orbit=multiply_orbit(Poly{{0,I(1)}});mpq_class Ef=expectation_product(f_orbit,Poly{{0,I(1)}},1);assert(Ef==mpq_class(643,616));assert(Edf7>0&&Ed2f7>0&&Ef7>0&&Ebf11>0&&Edf7>=Ed2f7);assert(engine->stats().computed==0);
 double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();std::ofstream out(argv[1]);out<<"{\"status\":\"one_second_order_shape_necessary_corner_only\",\"source_sha256\":\""<<source_sha<<"\",\"backend_sha256\":\""<<BACKEND<<"\",\"E_d_f7\":\""<<Edf7.get_str()<<"\",\"E_d2_f7\":\""<<Ed2f7.get_str()<<"\",\"E_f7\":\""<<Ef7.get_str()<<"\",\"E_b_f11\":\""<<Ebf11.get_str()<<"\",\"E_f\":\""<<Ef.get_str()<<"\",\"f7_orbit_terms\":"<<stored7.first.size()<<",\"f11_orbit_terms\":"<<power.size()<<",\"raw_orbit_crosscheck_through_power\":3,\"existing_cache_count\":"<<available.size()<<",\"new_moments_computed\":"<<engine->stats().computed<<",\"seconds\":"<<seconds<<"}\n";std::cout<<"wrote "<<argv[1]<<", seconds "<<seconds<<std::endl;
 }catch(const std::exception&e){std::cerr<<"checkpoint/error: "<<e.what()<<std::endl;return 3;}}
