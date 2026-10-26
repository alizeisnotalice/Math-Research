#!/usr/bin/env python3
"""Finite four-class pilot. Read-only040 arithmetic, one jet/center, exact clips."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod, comb
import argparse, importlib.util, gzip, hashlib, json, os, sys, time, resource
import flint
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'): sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,o):
 t=p.with_suffix('.tmp'); t.write_text(json.dumps(o,sort_keys=True,separators=(',',':'))+'\n'); os.replace(t,p)
def rq(x):
 v=F(str(x)); return fmpq(v.numerator,v.denominator)
def ep(x,upper=False): return str((x.upper() if upper else x.lower()).fmpq())
def enclosure(x): return {'ball':str(x),'lower':ep(x),'upper':ep(x,True)}
def display(s): return float(F(s))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True);ap.add_argument('--resume',action='store_true');args=ap.parse_args()
 began=time.monotonic();regpath=HERE/'registration.json';reg=json.loads(regpath.read_text());ctx.prec=reg['precision_bits']; assert ctx.prec==384
 for s in reg['sources']: assert sha(s['source'])==s['sha256']
 source=Path(reg['sources'][1]['source']); spec=importlib.util.spec_from_file_location('locked040',source);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 table=json.loads(gzip.decompress(Path(reg['sources'][0]['source']).read_bytes()));vertices=table['vertices'];simplices=table['simplices'];assert len(vertices)==280 and len(simplices)==3504
 work=args.work_dir.resolve();assert not any((p/'.git').exists() for p in (work,*work.parents));work.mkdir(parents=True,exist_ok=True)
 fingerprint=hashlib.sha256(json.dumps({'code':sha(__file__),'reg':sha(regpath),'sources':reg['sources']},sort_keys=True).encode()).hexdigest()
 if (work/'run_receipt.json').exists() and not args.resume: raise ValueError('Existing production requires --resume; do not rerun centers')
 receipt={'pid':os.getpid(),'began_unix':time.time(),'fingerprint':fingerprint,'source_sha256':sha(__file__),'registration_sha256':sha(regpath),'python':sys.executable,'flint':flint.__version__,'centers_completed':[]}
 atomic(work/'run_receipt.json',receipt);atomic(HERE/'run_receipt.json',receipt)
 alphas=sorted((x for x in product(range(4),repeat=4) if sum(x)<=3),key=lambda x:(sum(x),x));index={x:i for i,x in enumerate(alphas)};parents=[[(j,index[tuple(v-int(j==t) for t,v in enumerate(alpha))]) for j in range(4) if alpha[j]] for alpha in alphas]
 final=[]
 for row in reg['anchors']:
  name=row['name'];out=HERE/(name+'_results.json');cache=work/(name+'_checkpoint.json')
  if out.exists():
   if not args.resume: raise ValueError('Existing center result')
   r=json.loads(out.read_text());assert r['fingerprint']==fingerprint;final.append(r);continue
  cinfo=row['candidates_in_order'][row['selected_candidate']];center=list(map(F,cinfo['point']));anchor=list(map(F,row['nominal_anchor']));D=cinfo['D'];nodes=[sum(int(x*D)*z for x,z in zip(center,v)) for v in vertices];assert len(set(nodes))==280
  a=rq(144*D*F(reg['eps']));b=[rq(ep(arb(a*a+k*k).sqrt())) for k in nodes];assert min(b)>=a
  scans=[];qvs=[]
  for exponent in reg['halfwidth_exponents']:
   h=F(1,2**exponent);lo=[max(F(0),p-h) for p in anchor[1:]];hi=[min(F(1),p+h) for p in anchor[1:]];ul=[x-c for x,c in zip(lo,center[1:])];uh=[x-c for x,c in zip(hi,center[1:])];assert all(l<=0<=u for l,u in zip(ul,uh))
   radii=[max(abs(l),abs(u)) for l,u in zip(ul,uh)];numer=[]
   for v in vertices:
    low=D*sum(z*(l if z>=0 else u) for z,l,u in zip(v[1:],ul,uh)); high=D*sum(z*(u if z>=0 else l) for z,l,u in zip(v[1:],ul,uh));numer.append(max(abs(low),abs(high)))
   ratios2=[n*n/F(str(a*a+k*k)) for n,k in zip(numer,nodes)];vi=max(range(280),key=lambda i:ratios2[i]);passed=ratios2[vi]<1
   s={'exponent':exponent,'h':str(h),'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),'u_lower':list(map(str,ul)),'u_upper':list(map(str,uh)),'radii':list(map(str,radii)),'volume':str(prod(u-l for l,u in zip(lo,hi))),'guard_passed':passed,'maximum_exact_q_squared':str(ratios2[vi]),'worst_vertex_index':vi,'worst_vertex':vertices[vi],'guard_failure_is_exact_not_rounding':not passed}
   if passed:
    qv=[rq(n)/bb for n,bb in zip(numer,b)];assert max(qv)<1;s['max_q_rational_upper']=str(max(qv));s['max_q_display']=str(arb(max(qv)));qvs.append(qv)
   else:qvs.append(None)
   scans.append(s)
  atomic(HERE/(name+'_guards.json'),{'fingerprint':fingerprint,'pid':os.getpid(),'scans':scans})
  stencil=[{} for _ in alphas];next_cone=0;segments=[];center_began=time.monotonic()
  if cache.exists():
   assert args.resume;saved=json.loads(cache.read_text());assert saved['fingerprint']==fingerprint;next_cone=saved['next_simplex'];segments=saved['segments'];stencil=[{int(k):rq(v) for k,v in r.items()} for r in saved['stencil']]
  def save(done):
   atomic(cache,{'fingerprint':fingerprint,'center':name,'next_simplex':done,'stencil':[{str(k):str(v) for k,v in r.items()} for r in stencil],'segments':segments+[{'pid':os.getpid(),'seconds':time.monotonic()-center_began,'end_simplex':done}]})
  for ci in range(next_cone,3504):
   cone=simplices[ci];det=int(cone['det_numerator']);jet,keys=old.simplex_jet(cone['vertices'],vertices,nodes,D,a,alphas,parents)
   for ai,r in enumerate(jet):
    for sy,value in enumerate(r):
     if not value:continue
     key=0 if sy==0 else keys[sy-1];v=stencil[ai].get(key,fmpq(0))-det*value
     if v:stencil[ai][key]=v
     else:stencil[ai].pop(key,None)
   done=ci+1
   if done%128==0 or done==3504:save(done);print(name,done,'elapsed',round(time.monotonic()-began,3),flush=True)
   if time.monotonic()-began>=120 and done<3504:save(done);raise SystemExit(75)
  logs={k:arb(k).log() for r in stencil for k in r if k};factor=arb(D)**4/(144*(2*arb.pi())**5);coefs=[factor*sum((arb(v)*(logs[k] if k else 1) for k,v in r.items()),arb(0)) for r in stencil];assert coefs[0]>0
  # One all-time absolute integral per cone, shared by every registered h.
  time_cache={};rems=[arb(0) for _ in scans]
  for cone in simplices:
   ids=cone['vertices'];det=int(cone['det_numerator']);key=tuple(sorted(abs(nodes[i]) for i in ids))
   if key not in time_cache:time_cache[key]=old.time_absolute_upper(a,[b[i] for i in ids])
   base=det*time_cache[key]
   for j,qv in enumerate(qvs):
    if qv is None:continue
    q=max(qv[i] for i in ids);tail=(1-q)**(-5)-sum((fmpq(comb(n+4,4))*q**n for n in range(4)),fmpq(0));assert tail>=0;rems[j]+=arb(base*tail)
  for s,rem in zip(scans,rems):
   if not s['guard_passed']:continue
   rem=factor*arb(rq(ep(rem,True)));radii=list(map(rq,s['radii']));ul=list(map(rq,s['u_lower']));uh=list(map(rq,s['u_upper']));variation=sum((abs(c)*arb(prod(r**n for r,n in zip(radii,alpha))) for c,alpha in zip(coefs[1:],alphas[1:])),arb(0))
   rawlo=coefs[0]-variation-rem;rawhi=coefs[0]+variation+rem;L=max(fmpq(0),rq(ep(rawlo)));U=rq(ep(rawhi,True));assert U>L
   moments=[prod((u**(n+1)-l**(n+1))/((n+1)*(u-l)) for l,u,n in zip(ul,uh,alpha)) for alpha in alphas];mean=sum((c*arb(m) for c,m in zip(coefs,moments)),arb(0));ml=max(L,rq(ep(mean-rem)));mu=min(U,rq(ep(mean+rem,True)));assert L<=ml<=mu<=U
   q=arb(5)/4;slope=(arb(U)**q-arb(L)**q)/arb(U-L);chord=arb(rq(s['volume']))*(arb(L)**q+slope*arb(mu-L));constant=arb(rq(s['volume']))*arb(U)**q
   s.update({'H_range_lower':str(L),'H_range_upper':str(U),'raw_range_lower':enclosure(rawlo),'variation':enclosure(variation),'all_time_resolvent_remainder':enclosure(rem),'mean_Taylor':enclosure(mean),'mean_H_lower':str(ml),'mean_H_upper':str(mu),'exact_mean_moments':list(map(str,moments)),'chord_integral_upper':ep(chord,True),'chord_integral_expression':str(chord),'constant_range_integral_upper':ep(constant,True),'positive_lower_certificate':L>0,'remainder_le_one_percent_center':rq(ep(rem,True))*100<=rq(ep(coefs[0])),'relative_remainder_upper':ep(rem/arb(rq(ep(coefs[0]))),True)})
  r={'fingerprint':fingerprint,'center_name':name,'selected_center':cinfo,'nominal_anchor':row['nominal_anchor'],'eps':reg['eps'],'full3504_jet_calls':1,'all_time_interval':'[1/2048,infinity)','coefficients':[{'alpha':list(alpha),'enclosure':enclosure(c)} for alpha,c in zip(alphas,coefs)],'center_H':enclosure(coefs[0]),'scans':scans,'max_guard_h':next((s['h'] for s in scans if s['guard_passed']),None),'max_positive_lower_h':next((s['h'] for s in scans if s.get('positive_lower_certificate')),None),'max_one_percent_h':next((s['h'] for s in scans if s.get('remainder_le_one_percent_center')),None),'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-center_began,'total_seconds':time.monotonic()-began,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},'checkpoint':{'path':str(cache),'sha256':sha(cache),'bytes':cache.stat().st_size}}
  atomic(out,r);final.append(r);receipt['centers_completed'].append(name);atomic(HERE/'run_receipt.json',receipt);print('COMPLETE',name,r['max_guard_h'],r['max_positive_lower_h'],r['max_one_percent_h'],str(coefs[0]),flush=True)
  if time.monotonic()-began>=120 and len(final)<4:raise SystemExit(75)
 eps=arb(rq(reg['eps']));E=arb(5)/6*(arb(14455)/(31104*arb.pi()**5))**(arb(1)/5)*arb.pi()**(-arb(4)/5)*eps**(arb(6)/5);goal=arb(rq('1/25000'));gate=(goal**(arb(4)/5)-E)**(arb(5)/4);tube=arb(rq('1/1000000'));gate_tube=((goal-tube)**(arb(4)/5)-E)**(arb(5)/4);cost=arb(5)**(arb(5)/4)*goal/arb(rq('1/10000'))**(arb(5)/4)
 budget={'eps':reg['eps'],'E_eps_Lq_norm':enclosure(E),'I_target':'1/25000','g_conditional_target':'1/10000','UH_gate_zero_tube':enclosure(gate),'UH_gate_with_1e6_tube':enclosure(gate_tube),'I_if_UH_equals_target_zero_tube':enclosure((goal**(arb(4)/5)+E)**(arb(5)/4)),'conditional_cost':enclosure(cost),'exact_cost_less_30_fourth_power_check':3125*F(1,25000)**4<30**4*F(1,10000)**5,'none_of_UH_I_g_certified_globally':True}
 atomic(HERE/'budget.json',budget);receipt.update({'status':'complete_four_centers','elapsed_seconds':time.monotonic()-began,'new_global_bound':False});atomic(HERE/'run_receipt.json',receipt);atomic(work/'run_receipt.json',receipt)
 print(json.dumps({'budget':budget,'runtime':receipt},indent=2))
if __name__=='__main__':main()
