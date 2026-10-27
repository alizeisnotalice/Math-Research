"""Export actual selected radii, exact increments and rejected-cell flows for witness A/0/6."""
import os,sys,json,time,hashlib,resource
from pathlib import Path
from fractions import Fraction as Q
from joint_geometry import generate,STRATA,volumes
from joint_atoms import atoms,initial_flow
from r037_joint import K
sys.dont_write_bytecode=True
p=Path(__file__).resolve();d=json.loads(p.with_name('r037_failure_witness.json').read_text());start=time.monotonic();cpu=time.process_time();n,a,N,_=STRATA[0];cubes,_=generate('A',0,6);h=volumes(n,a,cubes);cells,potential=atoms(cubes,n,a);tau=K[6%4]/a**n
selected,ps,q,U,na=initial_flow(cells,cubes,n,a,d['tie_order'],tau);assert q==list(map(Q,d['q'])) and str(U)==d['union_volume']
reject={}
for (qm,em),vol in cells.items():
 if qm&selected:continue
 sig=em&selected;reject[sig]=reject.get(sig,Q())+vol
flow=[]
for sig,vol in sorted(reject.items()):
 denom=sum((v for i,v in ps.items() if sig>>i&1),Q());assert denom>tau
 flow.append({'eligible_indices':[i for i in ps if sig>>i&1],'volume':str(vol),'paying_p_sum':str(denom),'strictly_above_tau':True,'z':{str(i):str(vol*v/denom) for i,v in ps.items() if sig>>i&1}})
assert sum((h[i]*q[i] for i in range(N)),Q())==U
out={'case':'A/0/6/partial','cube_centers_and_radii':[{'center':list(map(str,z)),'radius':str(r)} for z,r in cubes],'h':[str(x) for x in h],'selected_indices':sorted(ps),'selected_p':{str(i):str(v) for i,v in ps.items()},'selected_increment_v':{str(i):str(v*h[i]) for i,v in ps.items()},'rejected_cell_flows':flow,'union_volume':str(U),'all_witnesses':'f=(11/10)*1_{union of all V_i}; h_i=|V_i|, Qi subset Vi','strict_payment_identity':True,'script_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'witness_sha256':hashlib.sha256(p.with_name('r037_failure_witness.json').read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
p.with_name('r037_failure_geometry.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS selected',sorted(ps),'rejected classes',len(flow),'sum payments exact')
