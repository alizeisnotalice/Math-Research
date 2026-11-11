// Exhaust all subsets of private cells. Input lists rectangle masks with integer masses.
#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
struct Box { uint32_t cells; unsigned mass; };
int main(int argc,char**argv){
 if(argc!=3)return 2;
 const auto begin=std::stoull(argv[1]),end=std::stoull(argv[2]);
 unsigned cells,count;if(!(std::cin>>cells>>count)||cells>32||cells==0||end>(uint64_t{1}<<cells)||begin>=end)return 2;
 std::vector<Box> boxes(count);for(auto&box:boxes)if(!(std::cin>>box.cells>>box.mass))return 2;
 std::sort(boxes.begin(),boxes.end(),[](Box a,Box b){return a.mass>b.mass;});
 unsigned best=0;uint32_t arg=0;uint64_t ties=0,checksum=0;
 for(uint64_t mask=begin;mask<end;++mask){
   uint32_t remaining=static_cast<uint32_t>(mask);unsigned score=0;
   for(const auto box:boxes){
     const auto hit=box.cells&mask;
     if((hit&remaining)!=0 && std::has_single_bit(hit)){
       score+=box.mass;remaining^=static_cast<uint32_t>(hit);if(!remaining)break;
     }
   }
   if(remaining)return 3;
   checksum+=score;
   if(score>best){best=score;arg=static_cast<uint32_t>(mask);ties=1;}else if(score==best)++ties;
 }
 std::cout<<begin<<' '<<end<<' '<<best<<' '<<arg<<' '<<ties<<' '<<checksum<<'\n';
}
