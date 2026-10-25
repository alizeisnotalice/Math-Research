#ifndef ASTRA015_GRAPH_MOMENTS_BETA_HALF_HPP
#define ASTRA015_GRAPH_MOMENTS_BETA_HALF_HPP
// Independent exact Eulerian K5 moments. B_k=D_|k| E_R r^k, beta=1/2.
// All B, D, C and binomial products are GMP integers. uint64_t is only a
// checked 50-bit graph key, never a moment accumulator.
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cerrno>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

namespace astra015 {
using Integer = mpz_class;
using Code = std::uint64_t;
using Counts = std::array<unsigned, 10>;
constexpr unsigned SHIFT=5, MASK=31, MAX_DEGREE=31;
inline constexpr std::array<unsigned,10> eu{{0,0,0,0,1,1,1,2,2,3}};
inline constexpr std::array<unsigned,10> ev{{1,2,3,4,2,3,4,3,4,4}};

class CheckpointStop : public std::runtime_error {
public: explicit CheckpointStop(const std::string &s):std::runtime_error(s){}
};
inline unsigned degree(const Counts &k){unsigned n=0;for(auto j:k)n+=j;return n;}
inline unsigned parity(const Counts &k){unsigned p=0;for(unsigned e=0;e<10;++e){
 if(k[e]&1U)p^=(1U<<eu[e])|(1U<<ev[e]);}return p;}
inline Code encode(const Counts &k){Code q=0;for(unsigned e=0;e<10;++e){
 if(k[e]>MASK)throw std::invalid_argument("edge exponent exceeds five-bit key");
 q|=Code(k[e])<<(SHIFT*e);}
 if(degree(k)>MAX_DEGREE)throw std::invalid_argument("total degree exceeds 31");return q;}
inline Counts decode(Code q){if(q>>(SHIFT*10))throw std::invalid_argument("invalid graph key high bits");
 Counts k{};for(unsigned e=0;e<10;++e)k[e]=(q>>(SHIFT*e))&MASK;
 if(degree(k)>MAX_DEGREE)throw std::invalid_argument("total degree exceeds 31");return k;}

struct Stats {
 std::size_t loaded=0,computed=0,canonical_cache_peak=0,subsets_visited=0;
 std::size_t journal_trailing_bytes_removed=0,max_B_bits=0,max_product_bits=0;
};

class Engine {
 std::array<Integer,MAX_DEGREE+1> denominator_,cumulant_;
 std::array<std::array<Integer,MAX_DEGREE+1>,MAX_DEGREE+1> choose_;
 std::vector<std::array<unsigned,10>> permutations_;
 std::unordered_map<Code,Code> canonical_cache_;
 std::unordered_map<Code,Integer> moments_;
 std::size_t canonical_limit_,max_new_,ticks_=0;
 Stats stats_;
 int journal_fd_=-1;
 std::chrono::steady_clock::time_point began_;
 double seconds_;
 std::string journal_path_;

 static std::uint64_t checksum(const std::string &s){std::uint64_t h=14695981039346656037ULL;
  for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;}
 static std::string hex(std::uint64_t x){std::ostringstream s;s<<std::hex<<std::setw(16)<<std::setfill('0')<<x;return s.str();}
 void append(const std::string &s){const char *p=s.data();std::size_t left=s.size();while(left){
  ssize_t n=::write(journal_fd_,p,left);if(n<0&&errno==EINTR)continue;
  if(n<=0)throw std::runtime_error("checkpoint write failed");p+=n;left-=std::size_t(n);}}
 void budget(){if((++ticks_&255U)!=0)return;
  if(seconds_>0&&elapsed()>=seconds_)throw CheckpointStop("time checkpoint reached");}
 void open_journal(const std::string &source_id){
  if(journal_path_.empty())return;
  if(source_id.size()!=64||source_id.find_first_not_of("0123456789abcdef")!=std::string::npos)
   throw std::invalid_argument("persistent cache requires a SHA256 source fingerprint");
  journal_fd_=::open(journal_path_.c_str(),O_RDWR|O_CREAT,0600);
  if(journal_fd_<0)throw std::runtime_error("cannot open checkpoint");
  if(::flock(journal_fd_,LOCK_EX|LOCK_NB)!=0)throw std::runtime_error("checkpoint already locked by another writer");
  std::ifstream in(journal_path_,std::ios::binary);std::string raw((std::istreambuf_iterator<char>(in)),{});
  const std::string header="astra015-K5-beta1/2-shift5-v1 "+source_id+"\n";
  if(raw.empty()){append(header);sync();return;}
  if(raw.rfind(header,0)!=0)throw std::runtime_error("checkpoint protocol/source fingerprint mismatch");
  auto end=raw.find_last_of('\n');
  if(end+1<raw.size()){stats_.journal_trailing_bytes_removed=raw.size()-end-1;
   if(::ftruncate(journal_fd_,off_t(end+1))!=0)throw std::runtime_error("checkpoint tail truncation failed");
   raw.resize(end+1);}
  std::istringstream rows(raw.substr(header.size()));std::string line;
  while(std::getline(rows,line)){
   auto a=line.find('\t'),b=line.find('\t',a==std::string::npos?0:a+1);
   if(a==std::string::npos||b==std::string::npos||line.substr(b+1)!=hex(checksum(line.substr(0,b))))
    throw std::runtime_error("checkpoint checksum/format mismatch");
   Code code=std::stoull(line.substr(0,a));Integer value(line.substr(a+1,b-a-1));auto k=decode(code);
   if(parity(k)||canonical(code)!=code||value<0||value>denominator_[degree(k)])
    throw std::runtime_error("checkpoint graph/range mismatch");
   auto found=moments_.find(code);if(found!=moments_.end()&&found->second!=value)
    throw std::runtime_error("checkpoint conflicting duplicate");
   if(found==moments_.end()){moments_.emplace(code,value);++stats_.loaded;}
  }
  if(::lseek(journal_fd_,0,SEEK_END)<0)throw std::runtime_error("checkpoint seek failed");
 }
 void save(Code code,const Integer &value){if(journal_fd_<0)return;
  std::string payload=std::to_string(code)+"\t"+value.get_str();append(payload+"\t"+hex(checksum(payload))+"\n");
  if((stats_.computed&63U)==0)sync();}
 void subsets(unsigned e,const Counts &k,Code sub,unsigned size,unsigned odd,
              const Integer &binomial,Code whole,Integer &total){
  budget();if(e==10){++stats_.subsets_visited;
   if(size&&odd==0){Integer term=cumulant_[size]*B(whole-sub)*binomial;
    stats_.max_product_bits=std::max(stats_.max_product_bits,mpz_sizeinbase(term.get_mpz_t(),2));total+=term;}
   return;}
  for(unsigned j=0;j<=k[e];++j){Integer next=binomial*choose_[k[e]][j];
   subsets(e+1,k,sub+(Code(j)<<(SHIFT*e)),size+j,
           odd^((j&1U)?((1U<<eu[e])|(1U<<ev[e])):0),next,whole,total);}
 }
public:
 explicit Engine(const std::string &journal="",const std::string &source_id="",
                 double checkpoint_seconds=120,std::size_t max_new_moments=0,
                 std::size_t canonical_cache_limit=500000):canonical_limit_(canonical_cache_limit),
  max_new_(max_new_moments),began_(std::chrono::steady_clock::now()),seconds_(checkpoint_seconds),journal_path_(journal){
  denominator_[0]=1;cumulant_[0]=0;
  for(unsigned n=1;n<=MAX_DEGREE;++n){denominator_[n]=denominator_[n-1]*(2*n-1);
   cumulant_[n]=n==1?Integer(1):Integer(cumulant_[n-1]*(2*n-3));}
  for(unsigned n=0;n<=MAX_DEGREE;++n){choose_[n][0]=choose_[n][n]=1;
   for(unsigned j=1;j<n;++j)choose_[n][j]=choose_[n-1][j-1]+choose_[n-1][j];}
  unsigned edge[5][5]{};for(unsigned e=0;e<10;++e)edge[eu[e]][ev[e]]=edge[ev[e]][eu[e]]=e;
  std::array<unsigned,5> p{{0,1,2,3,4}};do{std::array<unsigned,10> map{};
   for(unsigned e=0;e<10;++e)map[e]=edge[p[eu[e]]][p[ev[e]]];permutations_.push_back(map);
  }while(std::next_permutation(p.begin(),p.end()));
  moments_.emplace(0,Integer(1));canonical_cache_[0]=0;
  try{open_journal(source_id);}catch(...){if(journal_fd_>=0)::close(journal_fd_);journal_fd_=-1;throw;}
 }
 Engine(const Engine&)=delete;Engine&operator=(const Engine&)=delete;
 ~Engine(){if(journal_fd_>=0){::fsync(journal_fd_);::flock(journal_fd_,LOCK_UN);::close(journal_fd_);}}
 void sync(){if(journal_fd_>=0&&::fsync(journal_fd_)!=0)throw std::runtime_error("checkpoint fsync failed");}
 double elapsed()const{return std::chrono::duration<double>(std::chrono::steady_clock::now()-began_).count();}
 const Stats &stats()const{return stats_;}
 std::size_t cached_moments()const{return moments_.size();}
 const Integer &denominator(unsigned n)const{if(n>MAX_DEGREE)throw std::invalid_argument("denominator degree>31");return denominator_[n];}
 Code canonical(Code code){auto found=canonical_cache_.find(code);if(found!=canonical_cache_.end())return found->second;
  auto k=decode(code);Code best=code;for(const auto &p:permutations_){Code q=0;
   for(unsigned e=0;e<10;++e)q|=Code(k[e])<<(SHIFT*p[e]);best=std::min(best,q);}
  if(canonical_limit_&&canonical_cache_.size()>=canonical_limit_){canonical_cache_.clear();canonical_cache_[0]=0;}
  canonical_cache_[code]=best;stats_.canonical_cache_peak=std::max(stats_.canonical_cache_peak,canonical_cache_.size());return best;}
 Integer B(Code code){auto original=decode(code);if(parity(original))return Integer(0);code=canonical(code);
  auto found=moments_.find(code);if(found!=moments_.end())return found->second;
  if(max_new_&&stats_.computed>=max_new_)throw CheckpointStop("completed-submoment checkpoint reached");
  budget();auto k=decode(code);unsigned n=degree(k);Integer total=0;
  subsets(0,k,0,0,0,Integer(1),code,total);
  if(total<0||total>denominator_[n])throw std::runtime_error("exact Eulerian moment outside [0,D_n]");
  moments_.emplace(code,total);++stats_.computed;
  stats_.max_B_bits=std::max(stats_.max_B_bits,mpz_sizeinbase(total.get_mpz_t(),2));save(code,total);
  if(max_new_&&stats_.computed>=max_new_)throw CheckpointStop("completed-submoment checkpoint reached");
  return total;}
 Integer B(const Counts &k){return B(encode(k));}
 mpq_class expectation(const Counts &k){mpq_class q(B(k),denominator(degree(k)));q.canonicalize();return q;}
};
} // namespace astra015
#endif
