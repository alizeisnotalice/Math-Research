import json,math
from pathlib import Path
from fractions import Fraction as Q
import run,halo13
P=Path(__file__).resolve().parent
report={'files':0,'candidates':0,'finite':True,'exact_halo_agrees':True,'score_agrees':True,'rejected':[],'max_float_fraction_diff':0.0}
for f in sorted((P/'checkpoints').glob('*.json')):
 r=json.loads(f.read_text()); report['files']+=1
 for z in r['results']:
  report['candidates']+=1
  if 'fraction_score' not in z or not all(math.isfinite(z[k]) for k in ('fraction_score','float_score')):
   report['finite']=False;report['rejected'].append([f.name,z['segments'],'nonfinite or invalid geometry']);continue
  E=[tuple(map(Q,x)) for x in z['intervals']]; H=run.exact_halo(E,1-Q(1,2**r['k'])); H2=halo13.halo(E,1-Q(1,2**r['k']),True)
  if H!=H2 or H!=[tuple(map(Q,x)) for x in z['halo']]:report['exact_halo_agrees']=False;report['rejected'].append([f.name,z['segments'],'halo disagreement'])
  value,_=run.score(E,2**(-r['k']),r['tau'],True)
  if abs(value-z['fraction_score'])>1e-13:report['score_agrees']=False;report['rejected'].append([f.name,z['segments'],'score disagreement'])
  report['max_float_fraction_diff']=max(report['max_float_fraction_diff'],abs(z['float_score']-value))
(P/'audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

# Certificate acceptance requires the complete archived batch, not a vacuous pass.
inputs=[json.loads(f.read_text()) for f in sorted((P/'checkpoints').glob('*.json'))]
expected={(tau,k,rep) for tau in (.25,.5,.75) for k in (8,12,16) for rep in range(4)}
assert len(inputs)==36 and {(r['tau'],r['k'],r['rep']) for r in inputs}==expected
assert all(r['status']=='completed' for r in inputs)
assert report['files']==36 and report['candidates']==72
assert report['finite'] and report['exact_halo_agrees'] and report['score_agrees'] and not report['rejected']
for r in inputs:
 for z in r['results']:
  E=[tuple(map(Q,pair)) for pair in z['intervals']]
  assert len(E)==z['segments'] and all(a<b for a,b in E)
  assert all(E[i][1]<E[i+1][0] for i in range(len(E)-1))
print('ACCEPTED: complete batch and admissible disjoint interval inputs')
