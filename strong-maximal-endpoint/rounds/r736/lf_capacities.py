#!/usr/bin/env python3
"""Certified inner solve at every independently verified LF vertex."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[name]='1'
import json,time,sys,hashlib,importlib.util
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
HERE=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'r732'
spec=importlib.util.spec_from_file_location('capacity_inner',BASE/'capacity.py');cap=importlib.util.module_from_spec(spec);spec.loader.exec_module(cap)
def run(task):
 base,vertex,pitext=task;pi=list(map(F,pitext));mass=list(map(F,base['atom_masses']));vol=list(map(F,base['row_volumes']));masks=base['atom_masks'];inc=[[int(s>>i&1) for s in masks] for i in range(len(vol))]
 assert cap.hall_check(pi,mass,inc,masks)
 assert cap.testing_price(pi,mass,inc,vol,F(base['target']))==F(base['exact_NP_price'])
 oldw,olds=cap.lf_weights,cap.primal_source
 # These dependencies are process-local; the shared exact primal is replayed independently.
 cap.lf_weights=lambda *args:(pi,'verified LF vertex')
 cap.primal_source=lambda *args:(list(map(F,base['source_f'])),list(map(F,base['source_ell'])))
 try:
  rec=dict(boxes=base['boxes'],p=base['target'],source_envelope=base['exact_NP_price'],dimension=base['dimension'])
  out=cap.run((1000*base['case']+vertex,rec))
 finally:cap.lf_weights,cap.primal_source=oldw,olds
 out.update(source_case=base['case'],vertex=vertex);return out
if __name__=='__main__':
 start=time.monotonic();raw=(BASE/'capacity.json').read_bytes();bases=json.loads(raw)['results'];vr=(HERE/'lf_vertices.json').read_bytes();enum=json.loads(vr);verified=json.loads((HERE/'verify_face_wolfram.json').read_text())
 assert hashlib.sha256(raw).hexdigest()==enum['source_sha256']==verified['source_sha256']
 assert all(x['listed_vertices_exact_LF'] and x['all_LF_pi_in_listed_convex_hull'] for x in verified['results'])
 tasks=[(bases[x['case']],j,pi) for x in enum['results'] for j,pi in enumerate(x['unique_pi'])]
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,tasks))
 summaries=[]
 for b in bases:
  rr=[x for x in rows if x['source_case']==b['case']];lo=max(F(x['capacity_lower']) for x in rr);hi=max(F(x['capacity_upper']) for x in rr)
  assert hi>=F(b['capacity_lower'])
  summaries.append(dict(case=b['case'],target=b['target'],LF_vertex_count=len(rr),LF_capacity_lower=cap.decimal_bound(lo),LF_capacity_upper=cap.decimal_bound(hi,True),LF_capacity_interval_width_upper=cap.decimal_bound(hi-lo,True),improvement_over_original_fixed_pi_lower=cap.decimal_bound(lo-F(b['capacity_upper']))))
 out=dict(status='complete',scope='Complete LF-face finite catalogue capacities; not generalCn asymptotics',source_sha256=hashlib.sha256(raw).hexdigest(),vertices_sha256=hashlib.sha256(vr).hexdigest(),workers=12,log_series_terms=48,wall_seconds=time.monotonic()-start,vertex_certificates=rows,LF_capacity_summaries=summaries)
 (HERE/'lf_capacities.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(dict(status='complete',vertex_cases=len(rows),whole_faces=len(summaries),max_outer_width=max(float(x['LF_capacity_interval_width_upper']) for x in summaries),strict_outer_improvements=[(x['case'],x['improvement_over_original_fixed_pi_lower']) for x in summaries if F(x['improvement_over_original_fixed_pi_lower'])>0],wall_seconds=out['wall_seconds'])))
