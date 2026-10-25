#!/usr/bin/env python3
"""One registered rational radius pair per pending planar case; integer cells only."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import hashlib,json,os,resource,sys,threading,time,datetime
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
D=100000
INPUTS=['lower_round018/K8_n2_results.json','lower_round017/append_last_two/checkpoint_10of10.json','lower_round016/delivery007/case_evidence.json','lower_round016/delivery007/profile_setup_v4.json','lower_round016/delivery007/radius_certificate.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def save(p,d):p.write_text(json.dumps(d,separators=(',',':'))+'\n')
def root(n):
 den=2**80;lo=0;hi=8*den
 while lo+1<hi:
  m=(lo+hi)//2
  if m**5<=n*den**5:lo=m
  else:hi=m
 assert lo**5<=n*den**5<hi**5
 return Q(lo,den),Q(hi,den)
def runtime(start):
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def prepare():
 pending=load(PROJECT/INPUTS[0])['remaining_unresolved_cases'];assert len(pending)==16
 old=load(PROJECT/INPUTS[2]);profile=load(PROJECT/INPUTS[3]);assert Q(profile['critical_points'][0]['interval'][0])==4
 ul,uh=root(2);plans=[]
 for N in sorted({K*n for K,n in pending}):
  nl,nh=root(N);plo,phi=nl/(1+uh),nh/(1+ul)
  low=(plo*D).numerator//(plo*D).denominator;high=-((-phi*D).numerator//(-phi*D).denominator)
  if Q(low,D)==plo:low-=1
  if Q(high,D)==phi:high+=1
  assert Q(low,D)<plo<=phi<Q(high,D) and Q(high-low,D)<=Q(1,10000)
  assert Q(10*low,D).denominator!=1 and Q(10*high,D).denominator!=1
  assert nh<4 and 0<low<high<2*D
  plans.append({'N':N,'cases':[[K,n] for K,n in pending if K*n==N],'rho_grid_numerators':[low,high],'rho_grid_denominator':D,'rho_pair':[str(Q(low,D)),str(Q(high,D))],'strict_algebraic_rho_enclosure':[str(plo),str(phi)],'side_enclosure':[str(nl),str(nh)],'same_algebraic_radius_old_cases':[{'K':r['K'],'n':r['n'],'status':r['status'],'old_maximum':r.get('maximum_outer_overlap_if_available'),'original_result_sha256':r['original_result_sha256']} for r in old['cases'] if r['K']*r['n']==N]})
 assert len(plans)==14 and len({t for p in plans for t in p['rho_pair']})==28
 reg={'stage':'preregistered_16_pending_cases_14_unique_radius_pairs','registered_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pending_cases':pending,'plans':plans,'angle':[['4/5','-3/5'],['3/5','4/5']],'lattice':'s*R*Z^2, s=(1+2^(1/5))/2','common_radius':'(K*n)^(1/5)/2','normalized_radius':'(K*n)^(1/5)/(1+2^(1/5))','limits':{'seconds_per_arrangement':30,'candidate_cells_and_sweep_events_per_arrangement':10000,'maximum_arrangements':28,'stop_after_registered_pairs':True},'resources':{'nice_minimum':10,'numeric_threads':1,'single_process':True},'inputs':[{'path':p,'sha256':sha(PROJECT/p)} for p in INPUTS],'source_sha256':sha(Path(__file__)),'old_profile_scope':'previous completed radius profile has side>=4; every new side<4; old exact-radius evidence is cited, never rerun; new rational endpoints are distinct and shared by N','retained_old_case_only':'N384 had a completed old n12 capture calculation; new rational radii differ and diagnose the requested n8 threshold','stop_conditions':'save after each endpoint; on30 seconds or10000 nodes save resume checkpoint and exit; no additional pair, union, angle or a1 optimization','scope':'only rational planar min/max and >=n area bounds for the16 remaining finite cases; partial area brackets are not solved cases'}
 save(HERE/'registration.json',reg);print(json.dumps({'plans':[{k:p[k] for k in ('N','cases','rho_pair')} for p in plans],'maximum_arrangements':28,'registration_sha256':sha(HERE/'registration.json')},indent=2))
def make_atoms():
 reps=[(U,V,(4*U+3*V)//25,(-3*U+4*V)//25) for U,V in product(range(25),repeat=2) if (4*U+3*V)%25==0 and (-3*U+4*V)%25==0]
 assert len(reps)==25;atoms=[]
 for U,V,m,n in reps:
  for h,k in product((-1,0,1),repeat=2):
   mm,nn=m+4*h+3*k,n-3*h+4*k;x,y=(U+25*h)*(D//5),(V+25*k)*(D//5)
   assert Q(x,D)==Q(4*mm-3*nn,5) and Q(y,D)==Q(3*mm+4*nn,5)
   atoms.append({'id':len(atoms),'m':mm,'n':nn,'x_grid':x,'y_grid':y})
 assert len(atoms)==len({(p['m'],p['n']) for p in atoms})==225
 return atoms
class Limit(Exception):pass
def arrangement(N,mode,rnum,atoms):
 start=time.monotonic();nodes=0;dom=5*D
 def check():
  nonlocal nodes
  nodes+=1
  if nodes>10000 or time.monotonic()-start>30:raise Limit('30sec_or10000nodes')
 bs=[]
 for a in atoms:
  b=max(0,a['x_grid']-rnum),min(dom,a['x_grid']+rnum),max(0,a['y_grid']-rnum),min(dom,a['y_grid']+rnum)
  if b[0]<b[1] and b[2]<b[3]:bs.append((a['id'],b))
 xs=sorted({0,dom}|{x for i,b in bs for x in b[:2]});ys=sorted({0,dom}|{y for i,b in bs for y in b[2:]})
 xm=[sum(1<<i for i,b in bs if b[0]<=a and z<=b[1]) for a,z in zip(xs,xs[1:])]
 ym=[sum(1<<i for i,b in bs if b[2]<=a and z<=b[3]) for a,z in zip(ys,ys[1:])]
 hist={};trace=hashlib.sha256();minw=maxw=None;minc=225;maxc=0
 for xi,(xl,xh) in enumerate(zip(xs,xs[1:])):
  for yi,(yl,yh) in enumerate(zip(ys,ys[1:])):
   check();mask=xm[xi]&ym[yi];count=mask.bit_count();hist[count]=hist.get(count,0)+(xh-xl)*(yh-yl)
   trace.update(f'{xl},{xh},{yl},{yh}:{mask:x}\n'.encode())
   if count<minc:minc=count;minw=[xl,xh,yl,yh,mask]
   if count>maxc:maxc=count;maxw=[xl,xh,yl,yh,mask]
 # Independent exact y-event sweep checks the entire area histogram.
 hist2={}
 for xl,xh in zip(xs,xs[1:]):
  events={0:0,dom:0}
  for i,b in bs:
   if b[0]<=xl and xh<=b[1]:events[b[2]]=events.get(b[2],0)+1;events[b[3]]=events.get(b[3],0)-1
  count=0;last=0
  for y,delta in sorted(events.items()):
   check()
   if y>last:hist2[count]=hist2.get(count,0)+(xh-xl)*(y-last)
   count+=delta;last=y
  assert count==0
 assert hist==hist2 and sum(hist.values())==dom**2
 assert Q(10*rnum,D).denominator!=1
 def witness(w):
  xl,xh,yl,yh,mask=w
  return {'cell_grid':[xl,xh,yl,yh],'covering_atom_ids':[i for i in range(225) if mask>>i&1]}
 return {'status':'PASS_exact_rational_planar_arrangement','N':N,'mode':mode,'rho':str(Q(rnum,D)),'grid_denominator':D,'period_grid':[0,dom,0,dom],'atom_count':225,'active_windows':len(bs),'segments':[len(xs)-1,len(ys)-1],'elementary_cells':(len(xs)-1)*(len(ys)-1),'minimum_closed_period_capture':minc,'maximum_closed_period_capture':maxc,'closed_maximum_guard':'2rho not in(1/5)Z, so no simultaneous genuine start/end per coordinate; each closed event stratum is dominated by one neighboring open quadrant. The periodic seam has no genuine event. Hence open-cell and closed-period maxima agree.','count_area_grid_units':{str(k):v for k,v in sorted(hist.items())},'normalized_area_by_count':{str(k):str(Q(v,D*D)) for k,v in sorted(hist.items())},'min_witness':witness(minw),'max_witness':witness(maxw),'all_cell_mask_trace_sha256':trace.hexdigest(),'independent_y_event_histogram_agrees':True,'nodes':nodes,'runtime':runtime(start),'scope':'exact planar arrangement only at this one rational radius; boundaries included in extrema, area boundaries have measurezero'}
def run():
 start=time.monotonic();reg=load(HERE/'registration.json');assert sha(Path(__file__))==reg['source_sha256']
 for row in reg['inputs']:assert sha(PROJECT/row['path'])==row['sha256']
 print(json.dumps({'stage':'before_arrangements','runtime':runtime(start),'registered_cases':16,'unique_radius_pairs':14}),flush=True)
 if os.environ.get('LOWER019_HANDOFF')=='1':
  deadline=time.monotonic()+20
  while not (HERE/'start.flag').exists() and time.monotonic()<deadline:time.sleep(.05)
 atoms=make_atoms();ap=HERE/'atoms.json'
 if ap.exists():assert load(ap)==atoms
 else:save(ap,atoms)
 completed=[]
 for plan in reg['plans']:
  N=plan['N']
  for mode,rnum in zip(('inner','outer'),plan['rho_grid_numerators']):
   q=HERE/f'N{N}_{mode}.json'
   if q.exists():
    r=load(q);assert r['rho']==str(Q(rnum,D)) and r['status']=='PASS_exact_rational_planar_arrangement'
   else:
    try:r=arrangement(N,mode,rnum,atoms)
    except Limit as e:
     save(HERE/'checkpoint.json',{'status':'CHECKPOINT_LIMIT','next':[N,mode],'completed':completed,'reason':str(e),'runtime':runtime(start),'registration_sha256':sha(HERE/'registration.json')});return
    r['atom_manifest_sha256']=sha(ap);r['source_sha256']=sha(Path(__file__));r['registration_sha256']=sha(HERE/'registration.json');save(q,r)
   completed.append({'N':N,'mode':mode,'path':q.name,'sha256':sha(q)})
   save(HERE/'checkpoint.json',{'status':'saved_after_one_arrangement','completed':completed,'runtime':runtime(start),'registration_sha256':sha(HERE/'registration.json')})
   print(json.dumps({'N':N,'mode':mode,'rho':r['rho'],'min':r['minimum_closed_period_capture'],'max':r['maximum_closed_period_capture'],'cells':r['elementary_cells'],'nodes':r['nodes'],'seconds':r['runtime']['wall_seconds']}),flush=True)
 ul,uh=root(2);sl,sh=(1+ul)/2,(1+uh)/2;cases=[]
 for K,n in reg['pending_cases']:
  N=K*n;a=load(HERE/f'N{N}_inner.json');b=load(HERE/f'N{N}_outer.json')
  area=lambda row:sum((Q(v) for k,v in row['normalized_area_by_count'].items() if int(k)>=n),Q(0))
  low,high=area(a),area(b);assert 0<=low<=high<=25
  if b['maximum_closed_period_capture']<n:status='CERTIFIED_EMPTY_capture'
  elif a['minimum_closed_period_capture']>=n:status='CERTIFIED_FULL_PLANE_capture'
  else:status='PARTIAL_AREA_BRACKET_only'
  cases.append({'K':K,'n':n,'N':N,'rho_pair':[a['rho'],b['rho']],'inner_min_max':[a['minimum_closed_period_capture'],a['maximum_closed_period_capture']],'outer_min_max':[b['minimum_closed_period_capture'],b['maximum_closed_period_capture']],'normalized_capture_area_interval':[str(low),str(high)],'physical_capture_area_interval':[str(low*sl*sl),str(high*sh*sh)],'status':status,'is_case_resolved':status!='PARTIAL_AREA_BRACKET_only','radius_formula':'(K*n)^(1/5)/2 independent ofa1','evidence':[{'path':f'N{N}_{m}.json','sha256':sha(HERE/f'N{N}_{m}.json')} for m in ('inner','outer')]})
 out={'status':'COMPLETE_16_case_rational_pair_diagnostics','registered_case_count':16,'unique_mass_count':14,'completed_arrangements':len(completed),'cases':cases,'resolved_empty_cases':[[c['K'],c['n']] for c in cases if c['status']=='CERTIFIED_EMPTY_capture'],'resolved_full_cases':[[c['K'],c['n']] for c in cases if c['status']=='CERTIFIED_FULL_PLANE_capture'],'remaining_partial_cases':[[c['K'],c['n']] for c in cases if c['status']=='PARTIAL_AREA_BRACKET_only'],'s_interval':[str(sl),str(sh)],'atom_completeness':'25 distinct representatives ofRZ² modulo5Z², each with h,k in{-1,0,1}; every radius<2<5, so omitted translates cannot meet[0,5]²','distinct_mass':'a planar selectedn-tuple has distinct integer indices; product withK distinct unrotated points givesKn distinct atoms','same_common_radius':'all planar and unrotated inclusion tests use r=(K*n)^(1/5)/2; rho_inner*s<r<rho_outer*s','a1_dependence':'none in these planar captures; no new a1 event or objective was evaluated','source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),'arrangement_evidence':completed,'total_nodes':sum(load(HERE/r['path'])['nodes'] for r in completed),'runtime':runtime(start),'scope':'16 finite planar-case diagnostics only; partial brackets do not resolve exact capture cases; no five-dimensional union, a1 optimization, allscale optimum orC5 claim'}
 save(HERE/'results.json',out);save(HERE/'checkpoint.json',{'status':'COMPLETE_registered_batch_only','completed':completed,'result_sha256':sha(HERE/'results.json'),'runtime':runtime(start),'no_further_task_started':True})
 print(json.dumps({k:v for k,v in out.items() if k not in ('cases','arrangement_evidence')},indent=2),flush=True)
if __name__=='__main__':
 assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
 prepare() if len(sys.argv)>1 and sys.argv[1]=='--prepare' else run()
