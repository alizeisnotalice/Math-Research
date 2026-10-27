#!/usr/bin/env python3
"""R044: 76 cached final states, full mass LP and bounded three-coordinate probes."""
import os,sys,time,json,hashlib,resource,concurrent.futures,gzip,threading,subprocess,itertools,warnings,math
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
import numpy as np
from scipy.optimize import linprog,OptimizeWarning
D=Path(__file__).resolve().parent;BASE=D.parent.parent/'round038/compute';sys.path.insert(0,str(BASE))
from joint_core import rational_bounds,maskids,load
from replay_compact import step
sha=lambda b:hashlib.sha256(b).hexdigest()
canonical=lambda x:json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()
EXCLUDE={'phasea_C_2_05_full_disjoint.json','phasea_C_1_06_full_disjoint.json','phasea_A_3_15_partial.json'}

def register():
 packed=(BASE/'phasea_compact.json.gz').read_bytes();man=json.loads((BASE/'compact_manifest.json').read_bytes());assert sha(packed)==man['compressed_sha256'];raw=gzip.decompress(packed);assert sha(raw)==man['canonical_content_sha256'];cases=json.loads(raw)['cases'];chosen=[x for x in cases if x['raw_filename'] not in EXCLUDE];assert len(chosen)==76
 out={'batch':'R044','objective':'for each cached final payment perform one full mass LP proposal and rational primal/dual verification, plus at most 32 distinct genuine three-coordinate LP proposals and exact payment/load checks; no operator-bound claim','case_count':76,'workers':15,'blas_threads':1,'full_mass_LP_pairs':1,'max_triple_proposals':32,'stop':'finite list; no grid/iteration or idle load','excluded_three':sorted(EXCLUDE),'bundle_sha256':sha(packed),'manifest_sha256':sha((BASE/'compact_manifest.json').read_bytes()),'replay_source_sha256':sha((BASE/'replay_compact.py').read_bytes()),'core_source_sha256':sha((BASE/'joint_core.py').read_bytes()),'code_sha256':sha(Path(__file__).read_bytes()),'cases':[{'raw_filename':x['raw_filename'],'raw_sha256':x['raw_sha256'],'q_after_canonical_sha256':x['q_after_canonical_sha256']} for x in chosen],'registered_unix_time':time.time()}
 (D/'registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'registration_sha256':sha((D/'registration.json').read_bytes()),'code_sha256':out['code_sha256'],'cases':76,'workers':15},indent=2),flush=True)

def reconstruct(c):
 h=list(map(Q,c['witness_budgets']));q=list(map(Q,c['q_before']));C=[sum((1<<i for i in v),0) for v in c['maximal_cliques']];U=Q(c['U'])
 for i,j in c['accepted_pairs']:
  x,y,_,_=step(q,h,C,i,j);q[i]=x;q[j]=y
 assert sha(canonical(list(map(str,q))))==c['q_after_canonical_sha256'] and sum((h[i]*q[i] for i in range(len(h))),Q())==U and load(q,C)==Q(c['L_after'])
 return h,q,C,U

def proposals(q,h,C,L):
 active=sorted((i for i in range(len(q)) if q[i]>0),key=lambda i:(-h[i]*q[i],i))
 recipients=sorted(range(len(q)),key=lambda i:(-h[i],i));seen=set();out=[]
 for i in active[:8]:
  for j in list(dict.fromkeys(active[:8]+recipients[:16])):
   for k in recipients[:16]:
    tri=tuple(sorted((i,j,k)))
    if len(set(tri))!=3 or tri in seen:continue
    seen.add(tri);out.append(tri)
    if len(out)==32:return out
 return out

def triple(q,h,C,L,tri):
 # True 3-coordinate restricted minimax LP; all floats are proposals only.
 scale=max(h[i] for i in tri);hf=np.array([float(h[i]/scale) for i in tri]);qscale=L;pay=sum((h[i]*q[i] for i in tri),Q());b=float(pay/(scale*qscale));ix=set(tri)
 A=[];rhs=[]
 for mask in C:
  A.append([1.0 if mask>>i&1 else 0.0 for i in tri]+[-1.0]);other=sum((q[i] for i in maskids(mask) if i not in ix),Q());rhs.append(-float(other/qscale))
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',OptimizeWarning);sol=linprog([0,0,0,1],A_ub=np.array(A),b_ub=np.array(rhs),A_eq=[list(hf)+[0]],b_eq=[b],bounds=[(0,None)]*4,method='highs',options={'threads':1})
 if not sol.success:return None,'lp_failed'
 # Scale exact float proposals onto a rational grid; solve payment exactly for largest h.
 anchor=max(range(3),key=lambda k:h[tri[k]]);v=[None]*3
 for k in range(3):
  if k!=anchor:v[k]=Q(max(0,round(float(sol.x[k])*10**12)),10**12)*qscale
 residual=pay-sum((h[tri[k]]*v[k] for k in range(3) if k!=anchor),Q());v[anchor]=residual/h[tri[anchor]]
 if v[anchor]<0:return None,'rationalization_infeasible'
 qq=q.copy()
 for k,i in enumerate(tri):qq[i]=v[k]
 assert sum((h[i]*qq[i] for i in tri),Q())==pay
 Lnew=load(qq,C);return (Lnew,qq),'valid'

def work(c):
 start=time.monotonic();cpu=time.process_time();h,q,C,U=reconstruct(c);n=len(q);L=load(q,C)
 Mlo,Mhi,primal,dual=rational_bounds(h,C);assert Mlo>0 and all(sum((primal[j] for j,m in enumerate(C) if m>>i&1),Q())>=h[i] for i in range(n)) and all(sum((dual[i] for i in maskids(m)),Q())<=1 for m in C)
 globalq=[U*y/Mlo for y in dual];assert sum((h[i]*globalq[i] for i in range(n)),Q())==U;globalL=load(globalq,C);assert globalL<=U/Mlo
 triples=proposals(q,h,C,L);bestL=L;bestq=q;besttri=None;logs=[]
 for tri in triples:
  value,status=triple(q,h,C,L,tri);entry={'triple':list(tri),'status':status}
  if value is not None:
   Ln,qn=value;entry['load']=str(Ln)
   if Ln<bestL:bestL,bestq,besttri=Ln,qn,tri
  logs.append(entry)
 out={'source_raw':c['raw_filename'],'source_raw_sha256':c['raw_sha256'],'source_final_q_sha256':c['q_after_canonical_sha256'],'geometry':c['geometry'],'n':c['n'],'a':c['a'],'h':list(map(str,h)),'U':str(U),'maximal_cliques':c['maximal_cliques'],'q_initial':list(map(str,q)),'L_initial':str(L),'mass_primal':[str(x) for x in primal],'mass_dual':[str(x) for x in dual],'M_lower':str(Mlo),'M_upper':str(Mhi),'mass_primal_dual_feasible':True,'global_proposed_q':list(map(str,globalq)),'global_proposed_load':str(globalL),'global_improves':globalL<L,'triple_proposals':logs,'triple_count':len(logs),'best_triple':list(besttri) if besttri else None,'best_triple_q':list(map(str,bestq)),'best_triple_load':str(bestL),'triple_improves':bestL<L,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'finite payment candidate, not C_n bound; three-coordinate LP uses floating proposal with exact postcheck'}
 dest=D/('case_'+c['raw_filename']);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'case','name':c['raw_filename'],'pid':os.getpid(),'triples':len(logs),'global_improves':out['global_improves'],'triple_improves':out['triple_improves'],'wall':out['wall_seconds']}),flush=True);return {'file':dest.name,'sha256':sha(dest.read_bytes()),'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes'],'global_improves':out['global_improves'],'triple_improves':out['triple_improves']}

def monitor(stop,parent):
 path=D/'resource_monitor.jsonl';start=time.monotonic()
 with path.open('w') as f:
  while not stop.is_set():
   try:
    text=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);rows=[]
    for line in text.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:rows.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    record={'wall_seconds':time.monotonic()-start,'worker_count':len(rows),'total_cpu_percent':sum(x['cpu_percent'] for x in rows),'total_rss_kib':sum(x['rss_kib'] for x in rows),'workers':rows}
   except Exception as e:record={'wall_seconds':time.monotonic()-start,'monitor_error':repr(e)}
   f.write(json.dumps(record,separators=(',',':'))+'\n');f.flush();stop.wait(5)
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((D/'registration.json').read_bytes());packed=(BASE/'phasea_compact.json.gz').read_bytes();assert sha(Path(__file__).read_bytes())==reg['code_sha256'] and sha(packed)==reg['bundle_sha256'] and sha((BASE/'joint_core.py').read_bytes())==reg['core_source_sha256'] and sha((BASE/'replay_compact.py').read_bytes())==reg['replay_source_sha256'];cases=[x for x in json.loads(gzip.decompress(packed))['cases'] if x['raw_filename'] not in EXCLUDE];assert len(cases)==76
 start=time.monotonic();stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);thread.start();print(json.dumps({'event':'start','parent_pid':os.getpid(),'workers':15,'cases':76,'monitor':str(D/'resource_monitor.jsonl'),'registration_sha256':sha((D/'registration.json').read_bytes())}),flush=True);results=[]
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:
   futures=[ex.submit(work,c) for c in cases]
   for fut in concurrent.futures.as_completed(futures):results.append(fut.result())
 finally:stop.set();thread.join(timeout=10)
 summary={'batch':'R044','registration_sha256':sha((D/'registration.json').read_bytes()),'wall_seconds':time.monotonic()-start,'results':results};(D/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'event':'complete','cases':len(results),'wall_seconds':summary['wall_seconds']}),flush=True)
