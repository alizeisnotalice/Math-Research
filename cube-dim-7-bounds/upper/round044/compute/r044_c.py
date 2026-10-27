#!/usr/bin/env python3
"""R044-C: finite mixed-integer support proposals; exact primal payment checks."""
import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess,warnings
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint,OptimizeWarning
from scipy.sparse import lil_matrix
D=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
def register():
 summary=json.loads((D/'summary.json').read_bytes());rows=[];B=json.loads((D/'B_summary.json').read_bytes());bindex={r['file'].removeprefix('B_'):r for r in B['results']}
 for r in summary['results']:
  p=D/r['file'];raw=p.read_bytes();assert sha(raw)==r['sha256'];d=json.loads(raw);L=Q(d['L_initial']);G=Q(d['global_proposed_load']);gap=(L-G)/L
  rows.append({'file':p.name,'sha256':sha(raw),'relative_gap':str(gap),'B_status':bindex.get(p.name,{}).get('stop','not_in_B')})
 priority=[x for x in rows if x['B_status']=='512_limit']+ [x for x in rows if x['B_status']=='first_exact_improvement']
 assert len(priority)==15;priority.sort(key=lambda r:(r['B_status']!='512_limit',-Q(r['relative_gap']),r['file']));other=sorted([x for x in rows if x['B_status']=='not_in_B'],key=lambda r:(-Q(r['relative_gap']),r['file']))
 out={'batch':'R044-C','source_sha256':sha(Path(__file__).read_bytes()),'workers':15,'node_limit_per_case':50000,'mip_rel_gap':0,'wall_time_limit':None,'priorities':'13 B no-result then 2 B improvements, then remaining61 iff at least one first result returns within 5 s; no case is repeated','model':'min number changed q_i subject total payment, full maximal-clique loads <= midpoint(initial, known feasible global), q>=0; binary bigM uses valid bound q_i<=target since every vertex is in a clique','caveat':'MILP numerical optimality and dual bounds are diagnostic only; exact candidate verified separately','monitor':'0.5s initial, then 5s','cases':priority+other,'first_wave_count':15,'registered_unix_time':time.time()}
 (D/'C_registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'registered_C','cases_max':len(rows),'sha256':sha((D/'C_registration.json').read_bytes()),'first_wave':15}),flush=True)

def work(r):
 start=time.monotonic();cpu=time.process_time();raw=(D/r['file']).read_bytes();assert sha(raw)==r['sha256'];d=json.loads(raw);h=list(map(Q,d['h']));q=list(map(Q,d['q_initial']));g=list(map(Q,d['global_proposed_q']));C=d['maximal_cliques'];U=Q(d['U']);n=len(q);L=Q(d['L_initial']);G=Q(d['global_proposed_load']);assert sum((h[i]*q[i] for i in range(n)),Q())==U and sum((h[i]*g[i] for i in range(n)),Q())==U
 result=None;msg='no_strict_certified_gap';node_count=None;num_bound=None;num_gap=None;support_float=None;solver_status=None
 if G<L:
  target=(L+G)/2; q0=np.array([float(v/L) for v in q]);upper=float(target/L);hmax=max(h);hn=np.array([float(v/hmax) for v in h]);pay=float(U/(hmax*L))
  # A tiny stricter target only assists rationalizing a numerical solution; acceptance checks the exact requested target.
  cap=float((target-(L-G)/1000000)/L);rows=len(C)+2*n+1;A=lil_matrix((rows,2*n));lb=np.full(rows,-np.inf);ub=np.full(rows,np.inf)
  for j,cl in enumerate(C):A[j,cl]=1;ub[j]=cap
  for i in range(n):
   M=max(q0[i],max(0.,upper-q0[i]))+1e-14
   A[len(C)+2*i,i]=1;A[len(C)+2*i,n+i]=-M;ub[len(C)+2*i]=q0[i]
   A[len(C)+2*i+1,i]=-1;A[len(C)+2*i+1,n+i]=-M;ub[len(C)+2*i+1]=-q0[i]
  A[-1,:n]=hn;lb[-1]=pay;ub[-1]=pay
  with warnings.catch_warnings():
   warnings.simplefilter('ignore');sol=milp(np.r_[np.zeros(n),np.ones(n)],integrality=np.r_[np.zeros(n),np.ones(n)],bounds=Bounds(np.zeros(2*n),np.r_[np.full(n,upper),np.ones(n)]),constraints=LinearConstraint(A.tocsr(),lb,ub),options={'node_limit':50000,'mip_rel_gap':0,'threads':1})
  solver_status=int(sol.status);msg=str(sol.message);node_count=getattr(sol,'mip_node_count',None);num_bound=getattr(sol,'mip_dual_bound',None);num_gap=getattr(sol,'mip_gap',None)
  if sol.x is not None:
   support=[i for i in range(n) if sol.x[n+i]>.5];support_float=len(support)
   if support:
    anchor=max(support,key=lambda i:h[i]);v=q.copy()
    for i in support:
     if i!=anchor:v[i]=Q(max(0,round(float(sol.x[i])*10**12)),10**12)*L
    remainder=U-sum((h[i]*v[i] for i in range(n) if i!=anchor),Q());v[anchor]=remainder/h[anchor]
    if v[anchor]>=0 and all(x>=0 for x in v):
     actual=max(sum((v[i] for i in cl),Q()) for cl in C)
     if actual<=target:
      assert sum((h[i]*v[i] for i in range(n)),Q())==U;changed=[i for i in range(n) if v[i]!=q[i]]
      result={'q_feasible':list(map(str,v)),'exact_changed_indices':changed,'exact_support_count':len(changed),'exact_load':str(actual),'exact_target':str(target),'exact_payment_conserved':True,'exact_load_below_target':True}
 out={'source':r['file'],'source_sha256':r['sha256'],'status':msg,'solver_status':solver_status,'node_count_numerical':int(node_count) if node_count is not None else None,'numerical_dual_bound':float(num_bound) if num_bound is not None else None,'numerical_mip_gap':float(num_gap) if num_gap is not None else None,'numerical_support':support_float,'exact_candidate':result,'scope':'MILP optimality is numerical only; exact candidate verifies upper bound on support, never a proved minimum','pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dest=D/('C_'+r['file']);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'case_C','source':r['file'],'status':solver_status,'candidate':bool(result),'support':len(result['exact_changed_indices']) if result else None,'pid':out['pid'],'wall':out['wall_seconds']}),flush=True);return {'file':dest.name,'sha256':sha(dest.read_bytes()),'solver_status':solver_status,'candidate':bool(result),'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes']}
def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (D/'C_monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    txt=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in txt.splitlines():
     x=line.strip().split(None,5)
     if len(x)==6 and int(x[1])==parent and 'spawn_main' in x[5]:ws.append({'pid':int(x[0]),'cpu_percent':float(x[2]),'rss_kib':int(x[3]),'cpu_time':x[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(w['cpu_percent'] for w in ws),'total_rss_kib':sum(w['rss_kib'] for w in ws),'workers':ws}
   except Exception as e:row={'wall_seconds':time.monotonic()-start,'monitor_error':repr(e)}
   f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((D/'C_registration.json').read_bytes());assert sha(Path(__file__).read_bytes())==reg['source_sha256'] and all(sha((D/x['file']).read_bytes())==x['sha256'] for x in reg['cases']);start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start_C','pid':os.getpid(),'workers':15,'registered_cases':76,'registration_sha256':sha((D/'C_registration.json').read_bytes())}),flush=True)
 out=[];auto=False
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:
   futures={ex.submit(work,r) for r in reg['cases'][:15]};first=next(concurrent.futures.as_completed(futures));futures.remove(first);first_result=first.result();out.append(first_result);auto=first_result['wall_seconds']<=5
   if auto:futures.update(ex.submit(work,r) for r in reg['cases'][15:])
   for f in concurrent.futures.as_completed(futures):out.append(f.result())
 finally:stop.set();th.join(timeout=10)
 summary={'batch':'R044-C','registration_sha256':sha((D/'C_registration.json').read_bytes()),'additional_61_submitted':auto,'wall_seconds':time.monotonic()-start,'results':out};(D/'C_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'event':'complete_C','cases':len(out),'additional_61_submitted':auto,'wall_seconds':summary['wall_seconds']}),flush=True)
