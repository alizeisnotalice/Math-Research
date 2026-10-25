#define CUBE5_BUILD_COMMAND "embedded-by-009"
#define CUBE5_SOURCE_SHA256 "see-005b-locked-source"
#define main embedded_005b_main
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#include "005b_interval_heat.cpp"
#pragma clang diagnostic pop
#undef main
#ifndef CUBE9_SOURCE_SHA256
#error "Pass CUBE9_SOURCE_SHA256"
#endif
#ifndef CUBE9_TIME_SHA256
#error "Pass CUBE9_TIME_SHA256"
#endif
#ifndef CUBE9_SPACE_SHA256
#error "Pass CUBE9_SPACE_SHA256"
#endif
#ifndef CUBE9_BUILD_COMMAND
#error "Pass CUBE9_BUILD_COMMAND"
#endif
struct Node{I t,w;}; struct Box{int rid,cid;std::array<I,5>x;};
static double val(const std::string&s){return std::bit_cast<double>(std::uint64_t(std::stoull(s,nullptr,0)));}
static std::vector<Node> nodes(){std::ifstream f("009_time_nodes.bits");if(!f)throw std::runtime_error("time input");int p,j;std::string a,b,c,d;std::vector<Node>r;while(f>>p>>j>>a>>b>>c>>d){I t={val(a),val(b)},w={val(c),val(d)};if(!(t.lo<=t.hi&&w.lo>0&&w.lo<=w.hi))throw std::runtime_error("bad node");r.push_back({t,w});}if(r.size()!=225)throw std::runtime_error("need225");return r;}
static std::vector<Box> boxes(int wanted){std::ifstream f("009_space_boxes.bits");if(!f)throw std::runtime_error("space input");std::vector<Box>r;int rid,cid;while(f>>rid>>cid){std::array<I,5>x;for(int k=0;k<5;++k){std::string a,b,c,d;f>>a>>b>>c>>d;x[k]={val(a),val(d)};if(!(std::isfinite(x[k].lo)&&std::isfinite(x[k].hi)&&x[k].lo<=x[k].hi))throw std::runtime_error("bad box");}if(rid==wanted)r.push_back({rid,cid,x});}if(r.size()!=16)throw std::runtime_error("need16boxes");return r;}
static I heat(const std::array<I,5>&x,I t,const std::vector<Term>&terms){I s=exact(0);for(const auto&z:terms){C q={exact(double(z.det)),exact(0)};for(int j=0;j<5;++j){I dot=exact(0);for(int k=0;k<5;++k)dot=add(dot,mul(exact(double(z.ray[j][k])),x[k]));q=cmul(q,cinv({mul(t,exact(double(z.rate[j]))),neg(dot)}));}s=add(s,q.re);}I pi={down(0x1.921fb54442d18p+1),up(0x1.921fb54442d18p+1)},n=exact(1),tp=mul(exact(2),pi);for(int k=0;k<5;++k)n=mul(n,tp);return divide(s,n);}
int main(int argc,char**argv){try{if(argc!=2)throw std::runtime_error("usage region");int rid=std::stoi(argv[1]);if(rid<0||rid>3)throw std::runtime_error("region");if(fegetround()!=FE_TONEAREST)throw std::runtime_error("rounding");auto ns=nodes();auto bs=boxes(rid);auto ts=make_terms();std::string fp=fp_control();std::array<I,16>total;auto st=std::chrono::steady_clock::now();for(int i=0;i<16;++i){total[i]=exact(0);for(const auto&n:ns)total[i]=add(total[i],mul(n.w,heat(bs[i].x,n.t,ts)));}double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-st).count();std::ofstream f("009_region_"+std::to_string(rid)+"_heat.json");if(!f)throw std::runtime_error("openout");f<<"{\n\"endpoint_encoding\":\"IEEE754-binary64-bit-pattern\",\n\"region\":"<<rid<<",\"source_sha256\":\""<<CUBE9_SOURCE_SHA256<<"\",\n\"time_nodes_sha256\":\""<<CUBE9_TIME_SHA256<<"\",\n\"space_boxes_sha256\":\""<<CUBE9_SPACE_SHA256<<"\",\n\"build_command\":\""<<escape_json(CUBE9_BUILD_COMMAND)<<"\",\n\"fp_control\":\""<<fp<<"\",\n\"cells\":[";for(int i=0;i<16;++i){f<<(i?",":"")<<"{\"cell\":"<<bs[i].cid<<",\"input_box\":[";for(int k=0;k<5;++k){if(k)f<<",";write_interval(f,bs[i].x[k]);}f<<"],\"quadrature\":";write_interval(f,total[i]);f<<"}";}f<<"],\n\"nodes_per_cell\":225,\"terms\":"<<ts.size()<<",\"seconds\":"<<std::setprecision(17)<<sec<<"\n}\n";f.flush();if(!f)throw std::runtime_error("writeout");std::cout<<"region="<<rid<<" cells=16 seconds="<<sec<<" fp="<<fp<<"\n";}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}}
