from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,resource,os,sys,collections
sys.set_int_max_str_digits(0); st=time.monotonic();cpu=time.process_time();p=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round052/compute');names=['r052_integer.py','results.json','registration.json','equivalence.json','verify_top64.py','verify_top64.json','R052_REPORT.md'];raw={x:(p/x).read_bytes() for x in names};r=json.loads(raw['results.json']);reg=json.loads(raw['registration.json']);eq=json.loads(raw['equivalence.json']);P=F(3579,4096);eps=F(1,1024);den=2**20;T=[1,2,4,8,16,32,64,128];sha=lambda b:hashlib.sha256(b).hexdigest()
assert len(r['top64'])==64 and len(r['workers'])==15 and r['chunks_completed']==1800
assert reg['code_sha256']==sha(raw['r052_integer.py']) and reg['equivalence_sha256']==sha(raw['equivalence.json']) and r['registration_sha256']==sha(raw['registration.json'])
assert P*(1-F(1,512))*(1-eps)*(1-F(3,1024))>F(1,2)
assert (72*F(1,256)**2)<F(1,2) and (1-F(3,4))*P/(2*(1+eps))-2*eps>0
valid=sum(min(4096,(den-1)//(4*A))*len(T) for A in range(32,257));total=225*4096*len(T);assert r['totals']=={'tuples':total,'valid':valid,'d_guard':total-valid,'sigma_guard':0,'positive_deficit_guard':0}
assert sum(w['chunks'] for w in r['workers'].values())==1800 and sum(w['tuples'] for w in r['workers'].values())==total
rows=[];occupancy=collections.Counter();rankkeys=[]
for item in r['top64']:
 a,rho,s,d=map(F,[item['a'],item['rho'],item['s'],item['d']]);A=a*8;k=d*den;t=rho*16;assert all(x.denominator==1 for x in (A,k,t));A,k,t=map(int,(A,k,t))
 assert 32<=A<=256 and 1<=k<=4096 and t in T and a>3 and s==2*a+rho and 0<d<=F(1,16) and d<F(1,32)/a
 u=P/(2*(1+eps));v=1/(P*(1-d/2)*(1-eps));mass=F(63,64)*d/(2*(1+3*d/4));M=max(F(4),2*v/(1-3*d/4));br=(1-3/a)*u-2*eps;sig=s*d*d;assert M==4 and 0<sig<F(1,2) and br>0
 lim=mass*br/M;delta=(1-eps)*lim;pen=sig/(1-sig);value=delta-pen
 for key,x in {'u_lo':u,'u_hi':v,'mass':mass,'M':M,'bracket':br,'delta_lim':lim,'delta':delta,'sigma':sig,'penalty':pen,'c_raw':value}.items():assert F(item[key])==x
 rankkeys.append((-value,A,t,k));occupancy[(A,(k-1)//512)]+=1;rows.append({'a':str(a),'rho':str(rho),'s':str(s),'d':str(d),'c_raw':str(value)})
assert rankkeys==sorted(rankkeys) and len(set(rankkeys))==64 and max(occupancy.values())==21
assert rows[0]['a']=='73/8' and rows[0]['s']=='293/16' and rows[0]['d']=='1023/1048576' and F(rows[0]['c_raw'])>F(1,2**16)
# Independent boundary witnesses and exact integer reconstruction: no producer module imported.
tests=[]
for A,k,t in [(33,1,128),(255,1027,32),(73,1023,1),(64,4095,4),(256,1023,128),(256,1024,1)]:
 good=4*A*k<den
 if not good:tests.append({'A':A,'k':k,'t':t,'guard_reject':True});continue
 a=F(A,8);d=F(k,den);s=2*a+F(t,16);mass=F(63,64)*d/(2*(1+3*d/4));br=(1-3/a)*P/(2*(1+eps))-2*eps;c=(1-eps)*mass*br/4-s*d*d/(1-s*d*d)
 fn=512*3579*(A-24)-8200*A;pn=(4*A+t)*k*k;pd=16*den**2-pn;N=64449*k*fn;D=131072*(4*den+3*k)*4198400*A
 assert c==F(N*pd-pn*D,D*pd);tests.append({'A':A,'k':k,'t':t,'guard_reject':False,'c':str(c)})
out={'status':'PASS retained-point and counting audit; global completeness conditional on producer chunk execution','input_sha256':{k:sha(v) for k,v in raw.items()},'totals':r['totals'],'valid_formula':'8 sum_{A=32}^{256} min(4096, floor((2^20-1)/(4A)))','chunks':r['chunks_completed'],'tuples_per_chunk':512*8,'max_top64_chunk_occupancy':max(occupancy.values()),'top64':rows,'independent_boundary_checks':tests,'worker_cpu_seconds':sum(w['cpu_seconds'] for w in r['workers'].values()),'reported_wall_seconds':r['wall_seconds'],'monitor_bytes':(p/'monitor.jsonl').stat().st_size,'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'script_sha256':sha(Path(__file__).read_bytes())}
Path(__file__).with_name('audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('top64','input_sha256')},indent=2))
