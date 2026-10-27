#!/usr/bin/env python3
"""Reproduce one selected failure and export an independently checkable LP dual."""
import os,sys,json,time,hashlib,resource,math
from fractions import Fraction as Q
from pathlib import Path
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from joint_geometry import generate,STRATA,witness_graph,volumes
from joint_atoms import atoms,initial_flow
from joint_core import maximal_cliques,maskids,load
from joint_compress import compress,indices
from r037_joint import orders,K
p=Path(__file__).resolve();start=time.monotonic();cpu=time.process_time();template='A';si=0;j=6;n,a,N,_=STRATA[si]
old=next(json.loads(line) for line in (p.parent/'r037_A_0.jsonl').read_text().splitlines() if json.loads(line)['caseindex']==j);pattern=next(x for x in old['patterns'] if x['pattern']=='partial');assert pattern['status']=='H_FAIL_CERTIFIED'
cubes,retry=generate(template,si,j);adj=witness_graph(cubes,a);cs,nodes=maximal_cliques(adj);h=volumes(n,a,cubes);cell,potential=atoms(cubes,n,a);U=sum(cell.values(),Q())
best=None
for tie,order in enumerate(orders(cubes,template,si,j)):
 _,_,q,vol,_=initial_flow(cell,cubes,n,a,order);val=load(q,cs)
 if best is None or val>best[0]:best=(val,tie,order)
assert best[1]==old['best_tie_order'];tau=K[j%4]/a**n;sel,ps,q,total,natoms=initial_flow(cell,cubes,n,a,best[2],tau);assert total==U and [str(v) for v in q]==pattern['q_before'];L0=load(q,cs)
qp,moves,stop=compress(q,h,adj,cs);Lalg=load(qp,cs);assert [str(v) for v in qp]==pattern['q_after'] and moves==pattern['moves']
rows=[];cols=[]
for col,C in enumerate(cs):
 for i in maskids(C):rows.append(col);cols.append(i)
mat=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(cs),N));hf=np.array([float(v) for v in h]);hf/=max(hf)
sol=linprog(-hf,A_ub=mat,b_ub=np.ones(len(cs)),bounds=(0,None),method='highs');assert sol.success
den=10**12;y=[Q(max(0,math.floor(float(z)*den)),den) for z in sol.x];worst=max(sum((y[i] for i in maskids(C)),Q()) for C in cs)
if worst>1:y=[z/worst for z in y]
assert all(sum((y[i] for i in maskids(C)),Q())<=1 for C in cs)
Mlo=sum((h[i]*y[i] for i in range(N)),Q());assert Mlo>0;optup=U/Mlo
assert L0-optup>=L0/10000 and 2*Lalg>L0+optup
out={'status':'exact_H_failure_witness_pending_independent_audit','template':template,'stratum':si,'caseindex':j,'pattern':'partial','n':n,'a':a,'N':N,'generation_retry':retry,'tie_order_index':best[1],'tie_order':best[2],'tau':str(tau),'selected_indices':list(maskids(sel)),'union_volume':str(U),'h':[str(x) for x in h],'q':[str(x) for x in q],'maximal_cliques':[list(maskids(C)) for C in cs],'clique_count':len(cs),'bronkerbosch_nodes':nodes,'M_dual_feasible_vector':[str(v) for v in y],'M_lower_exact':str(Mlo),'Lopt_upper_exact':str(optup),'L0_exact':str(L0),'Lalg_exact':str(Lalg),'H_fail_margin_exact':str(2*Lalg-L0-optup),'H_gap_margin_exact':str(L0-optup-L0/10000),'accepted_moves':moves,'stopping_reason':stop,'elementary_cells':potential,'atom_classes':len(cell),'rejected_atom_classes':natoms,'source_hash':hashlib.sha256(p.read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
p.with_name('r037_failure_witness.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('status','template','stratum','caseindex','selected_indices','clique_count','stopping_reason','pid','wall_seconds')},indent=2))
