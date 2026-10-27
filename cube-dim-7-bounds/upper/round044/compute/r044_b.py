#!/usr/bin/env python3
"""R044-B: registered new triples, floating lex proposals with exact load/payment checks."""
import os,sys,time,json,hashlib,resource,concurrent.futures,itertools,threading,subprocess,warnings
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
import numpy as np
from scipy.optimize import linprog,minimize,OptimizeWarning
D=Path(__file__).resolve().parent;BASE=D.parent.parent/'round038/compute';sys.path.insert(0,str(BASE))
from joint_core import load
sha=lambda b:hashlib.sha256(b).hexdigest()

def register():
 import gzip
 compact=json.loads(gzip.decompress((BASE/'phasea_compact.json.gz').read_bytes()));stops={c['raw_filename']:c['stop'] for c in compact['cases']};rows=[]
 for p in D.glob('case_phasea_*.json'):
  raw=p.read_bytes();d=json.loads(raw);L=Q(d['L_initial']);G=Q(d['global_proposed_load'])
  if stops[d['source_raw']]=='stop_pair_coordinate' and G<L:
   retry=D/('retry_'+p.name);rr=json.loads(retry.read_bytes());assert rr['source_sha256']==sha(raw)
   rows.append({'file':p.name,'sha256':sha(raw),'retry_file':retry.name,'retry_sha256':sha(retry.read_bytes()),'relative_certified_feasible_gap':str((L-G)/L),'source_stop':'R038 nonedge coordinate stop; not R039 all-pair stationary'})
 rows.sort(key=lambda z:(-Q(z['relative_certified_feasible_gap']),z['file']));rows=rows[:15]
 out={'batch':'R044-B','selection':'15 largest relative exact gap to explicit feasible global candidate among R038 nonedge coordinate stops','source_sha256':sha(Path(__file__).read_bytes()),'phasea_bundle_sha256':sha((BASE/'phasea_compact.json.gz').read_bytes()),'workers_max':15,'max_new_triples_per_case':512,'old_32_triples':'excluded using fixed retry ledger','stop':'first exact strict load reduction, or 512 new triples exhausted','method':'8 clique-class float LP load minimization then float sum-q and quadratic energy tie proposals; all accepted payments and load comparisons exact; no proof of subproblem lex optimality','monitor':'0.5 second startup then 5 second snapshots','cases':rows,'registered_unix_time':time.time()};(D/'B_registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'registered_B','cases':len(rows),'sha256':sha((D/'B_registration.json').read_bytes()),'source_sha256':out['source_sha256']}),flush=True)

def worker(r):
 start=time.monotonic();cpu=time.process_time();raw=(D/r['file']).read_bytes();assert sha(raw)==r['sha256'];d=json.loads(raw);retryraw=(D/r['retry_file']).read_bytes();assert sha(retryraw)==r['retry_sha256'];retry=json.loads(retryraw);old={tuple(x['triple']) for x in retry['triple_logs']};h=list(map(Q,d['h']));q=list(map(Q,d['q_initial']));globalq=list(map(Q,d['global_proposed_q']));n=len(q);U=Q(d['U']);C=[sum((1<<i for i in c),0) for c in d['maximal_cliques']];L=Q(d['L_initial']);assert load(q,C)==L and load(globalq,C)<L and sum((h[i]*q[i] for i in range(n)),Q())==U
 # Rank coordinate influence by exact absolute global-candidate payment displacement.
 order=sorted(range(n),key=lambda i:(-abs(h[i]*(globalq[i]-q[i])),i));tri_iter=(tuple(sorted(t)) for t in itertools.combinations(order,3));logs=[];answer=None
 for tri in tri_iter:
  if tri in old:continue
  if len(logs)>=512:break
  scale=max(h[i] for i in tri);hf=np.array([float(h[i]/scale) for i in tri]);pay=sum((h[i]*q[i] for i in tri),Q());b=float(pay/(scale*L));ix=set(tri);groups={}
  for mask in C:
   tag=tuple((mask>>i)&1 for i in tri);outside=sum((q[i] for i in range(n) if mask>>i&1 and i not in ix),Q())
   if tag not in groups or outside>groups[tag]:groups[tag]=outside
  A=np.array([list(tag)+[-1.] for tag in groups]);rhs=np.array([-float(v/L) for v in groups.values()]);eq=np.array([list(hf)+[0.]])
  with warnings.catch_warnings():
   warnings.simplefilter('ignore',OptimizeWarning);s1=linprog([0,0,0,1],A_ub=A,b_ub=rhs,A_eq=eq,b_eq=[b],bounds=[(0,None)]*4,method='highs',options={'threads':1})
  entry={'triple':list(tri),'class_count':len(groups),'load_solver_status':int(s1.status)}
  if not s1.success:entry['status']='load_lp_failed';entry['message']=str(s1.message);logs.append(entry);continue
  # The requested target is lexicographic; floating tolerances only generate proposals.
  if s1.fun>=1-1e-12:entry['status']='no_float_strict_load_proposal';logs.append(entry);continue
  cap=min(1-1e-12,float(s1.fun)+1e-11);A2=np.vstack([A,[0,0,0,1]]);rhs2=np.r_[rhs,cap]
  with warnings.catch_warnings():
   warnings.simplefilter('ignore',OptimizeWarning);s2=linprog([1,1,1,0],A_ub=A2,b_ub=rhs2,A_eq=eq,b_eq=[b],bounds=[(0,None)]*4,method='highs',options={'threads':1})
  vec=s2.x if s2.success else s1.x;entry['sum_solver_status']=int(s2.status)
  if s2.success:
   sumcap=sum(vec[:3])+1e-11
   cons=[{'type':'eq','fun':lambda v:eq[0]@v-b},{'type':'ineq','fun':lambda v:rhs2-A2@v},{'type':'ineq','fun':lambda v:sumcap-sum(v[:3])}]
   s3=minimize(lambda v:float(np.dot(hf,v[:3]**2)),vec,method='SLSQP',bounds=[(0,None)]*4,constraints=cons,options={'maxiter':100,'ftol':1e-12});entry['energy_solver_success']=bool(s3.success)
   if s3.success:vec=s3.x
  anchor=max(range(3),key=lambda k:h[tri[k]]);vals=[None]*3
  for k in range(3):
   if k!=anchor:vals[k]=Q(max(0,round(float(vec[k])*10**12)),10**12)*L
  vals[anchor]=(pay-sum((h[tri[k]]*vals[k] for k in range(3) if k!=anchor),Q()))/h[tri[anchor]]
  if vals[anchor]<0:entry['status']='rationalization_infeasible';logs.append(entry);continue
  new=q.copy()
  for k,i in enumerate(tri):new[i]=vals[k]
  assert all(v>=0 for v in new) and sum((h[i]*new[i] for i in range(n)),Q())==U;Ln=load(new,C);entry['status']='exact_strict_improvement' if Ln<L else 'exact_no_load_improvement';entry['L_new']=str(Ln);logs.append(entry)
  if Ln<L:answer={'triple':list(tri),'q_new':list(map(str,new)),'L_new':str(Ln),'exact_payment_conserved':True,'strict_load_improvement':True};break
 out={'source_file':r['file'],'source_sha256':r['sha256'],'source_retry_sha256':r['retry_sha256'],'scope':'R038 nonedge stationary, not R039 all-pair stationarity; floating lex proposals, exact accepted load/payment','attempted':len(logs),'old_triples_excluded':len(old),'logs':logs,'result':answer,'stop':'first_exact_improvement' if answer else ('512_limit' if len(logs)==512 else 'candidate_domain_exhausted'),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dest=D/('B_'+r['file']);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'case_B','name':r['file'],'pid':out['pid'],'attempted':out['attempted'],'stop':out['stop'],'wall':out['wall_seconds']}),flush=True);return {'file':dest.name,'sha256':sha(dest.read_bytes()),'attempted':out['attempted'],'stop':out['stop'],'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes']}

def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (D/'B_monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    text=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in text.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws),'workers':ws}
   except Exception as e:row={'wall_seconds':time.monotonic()-start,'monitor_error':repr(e)}
   f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((D/'B_registration.json').read_bytes());assert sha(Path(__file__).read_bytes())==reg['source_sha256'] and all(sha((D/r['file']).read_bytes())==r['sha256'] and sha((D/r['retry_file']).read_bytes())==r['retry_sha256'] for r in reg['cases']);start=time.monotonic();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();print(json.dumps({'event':'start_B','pid':os.getpid(),'workers':min(15,len(reg['cases'])),'cases':len(reg['cases']),'registration_sha256':sha((D/'B_registration.json').read_bytes())}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:out=list(ex.map(worker,reg['cases']))
 finally:stop.set();th.join(timeout=10)
 summary={'batch':'R044-B','registration_sha256':sha((D/'B_registration.json').read_bytes()),'wall_seconds':time.monotonic()-start,'results':out};(D/'B_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'event':'complete_B','cases':len(out),'wall_seconds':summary['wall_seconds']}),flush=True)
