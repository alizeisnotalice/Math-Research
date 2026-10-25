#include "round012_graph_moments_gmp.hpp"
#include <cstdio>
#include <iostream>
#include <utility>

int main(int argc,char **argv){
 try{
  if(argc!=7){std::cerr<<"usage: runner REQUESTS JOURNAL SOURCE_SHA256 SECONDS MAX_NEW OUTPUT_JSON\n";return 2;}
  std::ifstream input(argv[1]);if(!input)throw std::runtime_error("cannot open requests");
  std::vector<std::pair<std::string,round012::Counts>> requests;std::string id;
  while(input>>id){if(id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)
    throw std::invalid_argument("request id must be alphanumeric, underscore or hyphen");
   round012::Counts k{};for(auto &j:k)if(!(input>>j))throw std::runtime_error("malformed request");
   if(round012::parity(k))throw std::invalid_argument("pilot expects Eulerian requests");
   round012::encode(k);requests.push_back({id,k});}
  round012::Engine engine(argv[2],argv[3],std::stod(argv[4]),std::stoull(argv[5]));
  std::vector<std::string> results;std::string status="completed",reason;
  try{for(const auto &request:requests){const auto &k=request.second;auto B=engine.B(k);auto D=engine.denominator(round012::degree(k));
    mpq_class moment(B,D);moment.canonicalize();std::ostringstream row;
    row<<"{\"id\":\""<<request.first<<"\",\"degree\":"<<round012::degree(k)<<",\"counts\":[";
    for(unsigned j=0;j<10;++j){if(j)row<<",";row<<k[j];}
    row<<"],\"canonical_code\":\""<<engine.canonical(round012::encode(k))<<"\",\"B\":\""<<B
       <<"\",\"D\":\""<<D<<"\",\"moment\":\""<<moment.get_str()<<"\",\"B_bits\":"<<mpz_sizeinbase(B.get_mpz_t(),2)
       <<",\"D_bits\":"<<mpz_sizeinbase(D.get_mpz_t(),2)<<"}";results.push_back(row.str());}
  }catch(const round012::CheckpointStop &stop){status="checkpoint";reason=stop.what();}
  engine.sync();auto s=engine.stats();std::string destination=argv[6],temporary=destination+".tmp";
  std::ofstream out(temporary);if(!out)throw std::runtime_error("cannot write result");
  out<<"{\"status\":\""<<status<<"\",\"checkpoint_reason\":\""<<reason<<"\",\"GMP_version\":\""<<gmp_version
     <<"\",\"elapsed_seconds\":"<<engine.elapsed()<<",\"completed_requests\":"<<results.size()<<",\"results\":[";
  for(std::size_t i=0;i<results.size();++i){if(i)out<<",";out<<results[i];}
  out<<"],\"cache\":{\"loaded\":"<<s.loaded<<",\"computed\":"<<s.computed
     <<",\"cached_moments\":"<<engine.cached_moments()<<",\"canonical_cache_peak\":"<<s.canonical_cache_peak
     <<",\"subsets_visited\":"<<s.subsets_visited<<",\"trailing_bytes_removed\":"<<s.journal_trailing_bytes_removed
     <<",\"max_B_bits\":"<<s.max_B_bits<<",\"max_product_bits\":"<<s.max_product_bits<<"}}\n";
  out.close();if(std::rename(temporary.c_str(),destination.c_str())!=0)throw std::runtime_error("atomic result rename failed");
  std::cout<<status<<": "<<results.size()<<" requests, "<<s.computed<<" new canonical submoments, "<<engine.elapsed()<<" seconds\n";
  return status=="completed"?0:75;
 }catch(const std::exception &e){std::cerr<<"round012 error: "<<e.what()<<"\n";return 1;}
}
