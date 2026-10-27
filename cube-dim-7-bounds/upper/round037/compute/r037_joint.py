#!/usr/bin/env python3
"""R037 registered finite geometric/payment experiment; exact outcomes only."""
import os,sys,time,json,hashlib,resource,concurrent.futures,traceback
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
from joint_geometry import TEMPLATES,STRATA,seed,generate,witness_graph,volumes
from joint_atoms import atoms,initial_flow,GeometryCensored
from joint_core import maximal_cliques,rational_bounds,load,Censored
from joint_compress import compress
PATH=Path(__file__).resolve();OUT=PATH.parent;K=(Q(1,2),Q(1),Q(3,2),Q(2));HASHES={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [PATH,OUT/'joint_geometry.py',OUT/'joint_atoms.py',OUT/'joint_core.py',OUT/'joint_compress.py']}

def orders(cubes,template,si,j):
 n=len(cubes)
 for v in range(8):
  def key(i):return (-cubes[i][1],seed(template,si,j,'tie-'+str(v)+'-'+str(i)))
  yield sorted(range(n),key=key)

def examine(pattern,q,h,U,adj,cliques,Mlo,Mhi):
 L0=load(q,cliques);lo=U/Mhi;hi=U/Mlo
 if L0-hi<Q(1,10000)*L0:return {'pattern':pattern,'status':'gap_below_threshold_or_unresolved','L0':str(L0),'opt_lower':str(lo),'opt_upper':str(hi)}
 qp,moves,stop=compress(q,h,adj,cliques);Lalg=load(qp,cliques)
 if Lalg<=(L0+lo)/2:status='H_pass_certified'
 elif Lalg>(L0+hi)/2:status='H_FAIL_CERTIFIED'
 else:status='H_interval_unresolved'
 return {'pattern':pattern,'status':status,'L0':str(L0),'Lalg':str(Lalg),'opt_lower':str(lo),'opt_upper':str(hi),'q_before':[str(x) for x in q] if status=='H_FAIL_CERTIFIED' else None,'q_after':[str(x) for x in qp] if status=='H_FAIL_CERTIFIED' else None,'moves':moves,'stop':stop}

def task(template,si):
 n,a,N,levels=STRATA[si];output=OUT/('r037_'+template+'_'+str(si)+'.jsonl');start=time.monotonic();cpu=time.process_time();stats={}
 with output.open('w') as f:
  for j in range(16):
   begin=time.monotonic();record={'template':template,'stratum':si,'caseindex':j,'seed_hex':format(seed(template,si,j,'geometry'),'064x'),'n':n,'a':a,'N':N,'levels':levels,'pid':os.getpid()}
   try:
    cubes,retries=generate(template,si,j);record['generation_retry']=retries
    if cubes is None:record['status']='geometry_unavailable'
    else:
     assert len(cubes)==N
     adj=witness_graph(cubes,a)
     if any(C.bit_count()==N for C in [sum((1<<i for i in range(N)),0)] if all((adj[i]|(1<<i))==C for i in range(N))):record['status']='common_star_skip'
     else:
      cliques,nodes=maximal_cliques(adj);record.update({'clique_count':len(cliques),'graph_nodes':nodes});h=volumes(n,a,cubes)
      if template=='C':
       U=sum(((2*r)**n for _,r in cubes),Q());q=[Q(1,a**n)]*N;Mlo,Mhi=rational_bounds(h,cliques);patterns=[examine('full_disjoint',q,h,U,adj,cliques,Mlo,Mhi)]
      else:
       cells,potential=atoms(cubes,n,a);record.update({'elementary_cells':potential,'atom_classes':len(cells)});U=sum(cells.values(),Q());best=None
       for tie,order in enumerate(orders(cubes,template,si,j)):
        selected,ps,q,tot,_=initial_flow(cells,cubes,n,a,order);assert tot==U and selected==(1<<N)-1
        L=load(q,cliques)
        if best is None or L>best[0]:best=(L,tie,order,q)
       Mlo,Mhi=rational_bounds(h,cliques);record['best_tie_order']=best[1];patterns=[examine('full',best[3],h,U,adj,cliques,Mlo,Mhi)]
       tau=K[j%4]/a**n;selected,ps,q,total,atoms_used=initial_flow(cells,cubes,n,a,best[2],tau);assert total==U
       record['partial_selected']=selected.bit_count();record['partial_rejected_atom_classes']=atoms_used;patterns.append(examine('partial',q,h,U,adj,cliques,Mlo,Mhi))
      record.update({'status':'completed','M_lower':str(Mlo),'M_upper':str(Mhi),'union_volume':str(U),'patterns':patterns})
   except (Censored,GeometryCensored) as ex:record['status']=str(ex)
   except Exception as ex:record.update({'status':'error','error':repr(ex),'traceback':traceback.format_exc()})
   record['case_wall_seconds']=time.monotonic()-begin;stats[record['status']]=stats.get(record['status'],0)+1
   f.write(json.dumps(record,separators=(',',':'))+'\n');f.flush()
   print(json.dumps({'event':'case','template':template,'stratum':si,'case':j,'status':record['status'],'pid':os.getpid(),'wall':record['case_wall_seconds']},separators=(',',':')),flush=True)
 return {'template':template,'stratum':si,'pid':os.getpid(),'stats':stats,'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'output':str(output)}
if __name__=='__main__':
 started=time.monotonic();tasks=[(a,i) for a in TEMPLATES for i in range(5)];print(json.dumps({'event':'start','batch':'R037','pid':os.getpid(),'workers':15,'tasks':15,'cases':240,'hashes':HASHES}),flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=15) as p:results=list(p.map(task,*zip(*tasks)))
 out={'batch':'R037','geometry_cases':240,'max_patterns':400,'workers':15,'hashes':HASHES,'wall_seconds':time.monotonic()-started,'results':results};(OUT/'r037_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','wall':out['wall_seconds'],'stats':[r['stats'] for r in results]}),flush=True)
