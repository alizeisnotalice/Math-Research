#!/usr/bin/env python3
"""Retry 76x32 failed three-coordinate proposals; no full LP is rerun."""
import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess,warnings
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
import numpy as np
from scipy.optimize import linprog,OptimizeWarning
D=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
def register():
 summary=json.loads((D/'summary.json').read_bytes());rows=summary['results'];assert len(rows)==76
 data={'batch':'R044 three-coordinate retry','reason':'initial LP initialized HiGHS default global scheduler; subsequent triple option threads=1 returned status 4, not mathematical infeasibility','no_full_LP_rerun':True,'source_sha256':sha(Path(__file__).read_bytes()),'inputs':[{'file':r['file'],'sha256':r['sha256']} for r in rows],'workers':15,'triples_per_case':'exactly the pre-registered list from original case; max32','registered_unix_time':time.time()};(D/'triple_retry_registration.json').write_text(json.dumps(data,indent=2)+'\n');print('retry registered SHA',sha((D/'triple_retry_registration.json').read_bytes()),flush=True)

def work(row):
 start=time.monotonic();cpu=time.process_time();raw=(D/row['file']).read_bytes();assert sha(raw)==row['sha256'];d=json.loads(raw);h=list(map(Q,d['h']));q=list(map(Q,d['q_initial']));L=Q(d['L_initial']);C=[sum((1<<i for i in c),0) for c in d['maximal_cliques']];U=Q(d['U']);logs=[];bestL=L;bestq=q;besttri=None
 assert all(v['status']=='lp_failed' for v in d['triple_proposals'])
 for item in d['triple_proposals']:
  tri=item['triple'];scale=max(h[i] for i in tri);pay=sum((h[i]*q[i] for i in tri),Q());hf=[float(h[i]/scale) for i in tri];b=float(pay/(scale*L));S=set(tri);A=[];rhs=[]
  for C0 in d['maximal_cliques']:
   A.append([float(i in C0) for i in tri]+[-1.0]);rhs.append(-float(sum((q[i] for i in C0 if i not in S),Q())/L))
  with warnings.catch_warnings():
   warnings.simplefilter('ignore',OptimizeWarning);sol=linprog([0,0,0,1],A_ub=np.array(A),b_ub=rhs,A_eq=[hf+[0]],b_eq=[b],bounds=[(0,None)]*4,method='highs',options={'threads':1})
  row0={'triple':tri,'solver_status':int(sol.status),'solver_message':str(sol.message)}
  if sol.success:
   anchor=max(range(3),key=lambda k:h[tri[k]]);v=[None]*3
   for k in range(3):
    if k!=anchor:v[k]=Q(max(0,round(float(sol.x[k])*10**12)),10**12)*L
   residual=pay-sum((h[tri[k]]*v[k] for k in range(3) if k!=anchor),Q());v[anchor]=residual/h[tri[anchor]]
   if v[anchor]>=0:
    qq=q.copy()
    for k,i in enumerate(tri):qq[i]=v[k]
    assert sum((h[i]*qq[i] for i in range(len(h))),Q())==U
    Ln=max(sum((qq[i] for i in c),Q()) for c in d['maximal_cliques']);row0.update({'status':'exact_feasible','load':str(Ln)})
    if Ln<bestL:bestL,bestq,besttri=Ln,qq,tri
   else:row0['status']='rationalization_infeasible'
  else:row0['status']='lp_failed'
  logs.append(row0)
 out={'source_file':row['file'],'source_sha256':row['sha256'],'no_full_LP_rerun':True,'triple_logs':logs,'best_triple':besttri,'best_q':list(map(str,bestq)),'L_initial':str(L),'L_best':str(bestL),'strict_improvement':bestL<L,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'floating proposals; exact feasibility/load only; not an optimality or operator certificate'}
 dest=D/('retry_'+row['file']);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'retry','source':row['file'],'pid':out['pid'],'valid':sum(x['status']=='exact_feasible' for x in logs),'improved':out['strict_improvement'],'wall':out['wall_seconds']}),flush=True);return {'file':dest.name,'sha256':sha(dest.read_bytes()),'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes'],'strict_improvement':out['strict_improvement'],'valid':sum(x['status']=='exact_feasible' for x in logs)}

def monitor(stop,parent):
 start=time.monotonic()
 with (D/'retry_resource_monitor.jsonl').open('w') as f:
  while not stop.is_set():
   try:
    text=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);workers=[]
    for line in text.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:workers.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(workers),'total_cpu_percent':sum(w['cpu_percent'] for w in workers),'total_rss_kib':sum(w['rss_kib'] for w in workers),'workers':workers}
   except Exception as e:row={'wall_seconds':time.monotonic()-start,'monitor_error':repr(e)}
   f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();stop.wait(5)
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((D/'triple_retry_registration.json').read_bytes());assert reg['source_sha256']==sha(Path(__file__).read_bytes());assert all(sha((D/x['file']).read_bytes())==x['sha256'] for x in reg['inputs']);start=time.monotonic();stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);thread.start();print(json.dumps({'event':'start_retry','pid':os.getpid(),'workers':15,'cases':76,'registration_sha256':sha((D/'triple_retry_registration.json').read_bytes())}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:result=list(ex.map(work,reg['inputs']))
 finally:stop.set();thread.join(timeout=10)
 out={'batch':'R044 triple retry','registration_sha256':sha((D/'triple_retry_registration.json').read_bytes()),'wall_seconds':time.monotonic()-start,'results':result};(D/'triple_retry_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete_retry','cases':len(result),'wall_seconds':out['wall_seconds']}),flush=True)
