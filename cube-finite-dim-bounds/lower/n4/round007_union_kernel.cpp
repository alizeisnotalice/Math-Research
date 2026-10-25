// Numerical complete4D box union, compressed endpoint event sweeps.
// Not an exact certificate: endpoints and measures are IEEE doubles.
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>
struct Node {int lo,hi,child;bool operator<(const Node&o) const{return std::tie(lo,hi,child)<std::tie(o.lo,o.hi,o.child);}};
struct Hash {size_t operator()(const std::vector<int>&v)const{size_t h=1469598103934665603ull;for(int x:v){h^=size_t(x+1);h*=1099511628211ull;}return h;}};
struct Event {int pos,id,delta;bool operator<(const Event&o)const{return std::tie(pos,id,delta)<std::tie(o.pos,o.id,o.delta);}};
class Union4 {
public:
 std::array<std::vector<double>,4> coords;
 std::array<std::vector<Node>,4> nodes;
 std::array<std::unordered_map<std::vector<int>,double,Hash>,4> cache;
 long long states=0;
 double measure(int dim,const std::vector<int>&ids){
  if(ids.empty())return 0;
  if(dim==3){int lo=nodes[3][ids[0]].lo,hi=nodes[3][ids[0]].hi;double value=0;
   for(size_t k=1;k<ids.size();++k){const auto&b=nodes[3][ids[k]];if(b.lo>hi){value+=coords[3][hi]-coords[3][lo];lo=b.lo;hi=b.hi;}else hi=std::max(hi,b.hi);}return value+coords[3][hi]-coords[3][lo];}
  if(dim>0){auto it=cache[dim].find(ids);if(it!=cache[dim].end())return it->second;}
  ++states;std::vector<Event> events;events.reserve(2*ids.size());
  for(int id:ids){const auto&b=nodes[dim][id];events.push_back({b.lo,b.child,1});events.push_back({b.hi,b.child,-1});}
  std::sort(events.begin(),events.end());std::vector<int> counts(nodes[dim+1].size(),0);std::set<int> active;double value=0,compensation=0;size_t k=0;
  while(k<events.size()){
   int pos=events[k].pos;
   while(k<events.size()&&events[k].pos==pos){auto e=events[k++];int before=counts[e.id];counts[e.id]+=e.delta;if(before==0&&counts[e.id]>0)active.insert(e.id);else if(before>0&&counts[e.id]==0)active.erase(e.id);}
   if(k==events.size()||active.empty())continue;
   int next=events[k].pos;if(next<=pos)continue;
   std::vector<int> child(active.begin(),active.end());double term=(coords[dim][next]-coords[dim][pos])*measure(dim+1,child);
   double corrected=term-compensation,total=value+corrected;compensation=(total-value)-corrected;value=total;
  }
  if(dim>0)cache[dim].emplace(ids,value);return value;
 }
 double run(const double*box,int n){
  for(int d=0;d<4;++d){for(int i=0;i<n;++i){coords[d].push_back(box[i*8+2*d]);coords[d].push_back(box[i*8+2*d+1]);}std::sort(coords[d].begin(),coords[d].end());coords[d].erase(std::unique(coords[d].begin(),coords[d].end()),coords[d].end());}
  std::vector<std::array<int,8>> endpoints(n);for(int i=0;i<n;++i)for(int d=0;d<4;++d)for(int side=0;side<2;++side)endpoints[i][2*d+side]=int(std::lower_bound(coords[d].begin(),coords[d].end(),box[i*8+2*d+side])-coords[d].begin());
  std::vector<int> child(n,-1);
  for(int d=3;d>=0;--d){std::map<Node,int> unique;for(int i=0;i<n;++i)unique.emplace(Node{endpoints[i][2*d],endpoints[i][2*d+1],child[i]},0);int number=0;for(auto&entry:unique){entry.second=number++;nodes[d].push_back(entry.first);}for(int i=0;i<n;++i)child[i]=unique.at(Node{endpoints[i][2*d],endpoints[i][2*d+1],child[i]});}
  std::vector<int> root;for(int i=0;i<int(nodes[0].size());++i)root.push_back(i);return measure(0,root);
 }
};
extern "C" double union4_volume(const double*box,int n,long long*states){if(n==0){*states=0;return 0;}Union4 u;double value=u.run(box,n);*states=u.states;return value;}
