#!/usr/bin/env python3
"""R038-A: balance 79 locked historical failure patterns, not new geometry search."""
import os,sys,time,json,hashlib,resource,concurrent.futures,importlib.util,glob
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
from joint_core import maximal_cliques,load
from joint_geometry import geometry_fingerprint
from pair_balance import balance
P=Path(__file__).resolve();D=P.parent;OLD=D.parent.parent/'round037/compute';GEOM=OLD/'joint_geometry.py'
spec=importlib.util.spec_from_file_location('r037_geometry_frozen',GEOM);legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inputs():
 rows=[]
 for f in sorted(OLD.glob('r037_[ABC]_[0-4].jsonl')):
  for line in f.read_text().splitlines():
   r=json.loads(line)
   for v in r.get('patterns',[]):
    if v['status']=='H_FAIL_CERTIFIED':rows.append({'template':r['template'],'stratum':r['stratum'],'caseindex':r['caseindex'],'pattern':v['pattern'],'source_jsonl':f.name,'q_before':v['q_before'],'U':r['union_volume'],'M_lower_history':r['M_lower'],'M_upper_history':r['M_upper'],'L0_history':v['L0']})
 assert len(rows)==79 and len({(r['template'],r['stratum'],r['caseindex'],r['pattern']) for r in rows})==79
 return rows

def register():
 rows=inputs();oldfiles=sorted(OLD.glob('r037_[ABC]_[0-4].jsonl'));paths=[P,D/'pair_balance.py',D/'joint_core.py',D/'joint_geometry.py',GEOM]
 reg={'batch':'R038 Phase A','authorization':'79 existing R037 H_FAIL pattern inputs; no new geometry search or baseline recomputation','workers':15,'tasks':79,'accepted_step_limit':128,'pair_policy':'full nonedge lex(L,deltaE,i,j); plateau minimum energy; skip zero-total pairs','graph_cap_nodes':2**20,'graph_cap_cliques':65536,'source_hashes':{str(p):sha(p) for p in paths},'input_jsonl_hashes':{p.name:sha(p) for p in oldfiles},'identities':[{k:r[k] for k in ('template','stratum','caseindex','pattern','source_jsonl')} for r in rows],'history_bounds':'informational only, not rigorous H conclusions without new separately verified LP certificate','source_reconstruction':'R037 exact historical generator only; no atomization, selection, flow, or baseline LP'}
 (D/'phasea_registration.json').write_text(json.dumps(reg,indent=2)+'\n');print(json.dumps({'registered':79,'source_hashes':reg['source_hashes'],'registration_sha256':sha(D/'phasea_registration.json')},indent=2))

def task(r):
 start=time.monotonic();cpu=time.process_time();template=r['template'];si=r['stratum'];j=r['caseindex'];key=f"{template}_{si}_{j:02d}_{r['pattern']}";n,a,N,_=legacy.STRATA[si];cubes,retry=legacy.generate(template,si,j);assert cubes is not None and len(cubes)==N
 adj=legacy.witness_graph(cubes,a);cliques,nodes=maximal_cliques(adj);h=legacy.volumes(n,a,cubes);q=list(map(Q,r['q_before']));U=Q(r['U']);assert len(q)==N and sum((h[i]*q[i] for i in range(N)),Q())==U
 oldL=load(q,cliques);assert oldL==Q(r['L0_history'])
 qnew,steps,stop,evals=balance(q,h,adj,cliques);newL=load(qnew,cliques);assert newL<=oldL and sum((h[i]*qnew[i] for i in range(N)),Q())==U
 out={'batch':'R038 Phase A','historical_identity':{k:r[k] for k in ('template','stratum','caseindex','pattern','source_jsonl')},'n':n,'a':a,'N':N,'geometry_fingerprint':geometry_fingerprint(cubes,n),'geometry': [{'center':list(map(str,z)),'radius':str(v)} for z,v in cubes],'generation_retry':retry,'witness_budgets':[str(v) for v in h],'maximal_cliques':[[i for i in range(N) if C>>i&1] for C in cliques],'graph_nodes':nodes,'U':str(U),'q_before':[str(v) for v in q],'q_after':[str(v) for v in qnew],'L_before':str(oldL),'L_after':str(newL),'steps':steps,'stop':stop,'pair_evaluations':evals,'history_M_lower_unverified':r['M_lower_history'],'history_M_upper_unverified':r['M_upper_history'],'status':'completed; historical LP interval requires separate independent certificate','pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dest=D/('phasea_'+key+'.json');dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'pattern','id':key,'stop':stop,'steps':len(steps),'pair_evaluations':evals,'pid':os.getpid(),'wall':out['wall_seconds']}),flush=True)
 return {'id':key,'output':dest.name,'sha256':sha(dest),'L_before':str(oldL),'L_after':str(newL),'steps':len(steps),'stop':stop,'pid':os.getpid(),'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes']}

if __name__=='__main__':
 if sys.argv[1:] == ['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((D/'phasea_registration.json').read_text());assert all(sha(Path(p))==digest for p,digest in reg['source_hashes'].items());assert all(sha(OLD/p)==digest for p,digest in reg['input_jsonl_hashes'].items());rows=inputs();start=time.monotonic();print(json.dumps({'event':'start','parent_pid':os.getpid(),'workers':15,'patterns':len(rows),'registration_sha256':sha(D/'phasea_registration.json')}),flush=True)
 results=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=15) as pool:
  futures=[pool.submit(task,r) for r in rows]
  for f in concurrent.futures.as_completed(futures):results.append(f.result())
 out={'batch':'R038 Phase A','registration_sha256':sha(D/'phasea_registration.json'),'wall_seconds':time.monotonic()-start,'results':results};(D/'phasea_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','cases':len(results),'wall':out['wall_seconds']}),flush=True)
