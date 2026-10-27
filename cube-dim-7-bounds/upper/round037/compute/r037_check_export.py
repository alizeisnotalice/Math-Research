"""Check rational frozen witness ledger without scipy or search."""
import sys,json,time,hashlib,os,resource
from fractions import Fraction as Q
from pathlib import Path
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
p=Path(__file__).resolve();w=json.loads(p.with_name('r037_failure_witness.json').read_text());g=json.loads(p.with_name('r037_failure_geometry.json').read_text());start=time.monotonic();cpu=time.process_time()
h=list(map(Q,w['h']));q=list(map(Q,w['q']));y=list(map(Q,w['M_dual_feasible_vector']));M=[tuple(c) for c in w['maximal_cliques']];U=Q(w['union_volume']);tau=Q(w['tau']);assert all(v>=0 for v in q+y)
assert all(sum((y[i] for i in C),Q())<=1 for C in M);Mlo=sum((h[i]*y[i] for i in range(len(h))),Q());assert Mlo==Q(w['M_lower_exact'])
assert sum((h[i]*q[i] for i in range(len(h))),Q())==U;L0=max(sum((q[i] for i in C),Q()) for C in M);assert L0==Q(w['L0_exact'])==Q(w['Lalg_exact']);assert L0-U/Mlo>=L0/10000 and 2*L0>L0+U/Mlo
pvals={int(i):Q(v) for i,v in g['selected_p'].items()};inc={int(i):Q(v) for i,v in g['selected_increment_v'].items()};assert all(pvals[i]*h[i]==inc[i] for i in pvals)
charges={i:inc[i] for i in pvals}
for a in g['rejected_cell_flows']:
 inds=a['eligible_indices'];v=Q(a['volume']);P=sum((pvals[i] for i in inds),Q());assert P>tau and P==Q(a['paying_p_sum']);assert sum((Q(x) for x in a['z'].values()),Q())==v
 for i in inds:
  z=Q(a['z'][str(i)]);assert z==v*pvals[i]/P;charges[i]+=z
assert all(charges.get(i,Q())==h[i]*q[i] for i in range(len(h)))
assert all(h[i]==max(h) for i in range(len(h)) if q[i]>0)
out={'status':'PASS','scope':'exact exported payments, dual feasibility, gap/H arithmetic and zero-saving sufficient structure; independent geometric atomization and clique completeness still require audit','witness_sha256':hashlib.sha256(p.with_name('r037_failure_witness.json').read_bytes()).hexdigest(),'geometry_sha256':hashlib.sha256(p.with_name('r037_failure_geometry.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
p.with_name('r037_check_export.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
