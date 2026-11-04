import json, math
from fractions import Fraction as Q
from pathlib import Path
from run import exact_halo, difference
P=Path(__file__).resolve().parent
files=sorted((P/'checkpoints').glob('*.json')); records=[json.loads(x.read_text()) for x in files]; winners=[]; statements=[]
for tau in [.25,.5,.75]:
 for k in [8,12,16]:
  rs=[x for x in records if x['tau']==tau and x['k']==k]
  vals=[(z['fraction_score'],r,z) for r in rs for z in r['results'] if 'fraction_score' in z and math.isfinite(z['fraction_score']) and math.isfinite(z['float_score'])]
  if not vals:continue
  v,r,z=max(vals,key=lambda x:x[0]); base=r['baseline']; chosen={'tau':tau,'k':k,'tasks_finished':len(rs),'baseline':base,'raw_best':v,'envelope':max(v,base),'delta':v-base,'seed':r['seed'],'segments':z['segments'],'intervals':z['intervals'],'halo':z['halo']}; winners.append(chosen)
  E=[tuple(map(Q,x)) for x in z['intervals']];H=exact_halo(E,1-Q(1,2**k));assert H==[tuple(map(Q,x)) for x in z['halo']]
  D=difference(H,E); q=Q(str(tau));t=Q(1,2**k)
  def m(S):return '+'.join(f'(Sign[{b}] Abs[{b}]^({q})-Sign[{a}] Abs[{a}]^({q}))' for a,b in S) or '0'
  A='(('+m(D)+')/('+m(E)+')/('+str(t)+')^('+str(q)+'))'; B='(1/((1+'+str(t)+')^('+str(q)+')-('+str(t)+')^('+str(q)+')))'
  statements.append(f'Print["tau={tau},k={k}, numeric delta=", N[{A}-{B},50], ", exact greater=", FullSimplify[{A}>{B}]];')
(P/'summary.json').write_text(json.dumps({'completed':len(records),'evaluations':sum(z['evaluations'] for r in records for z in r['results']),'winners':winners},indent=2))
(P/'exact_comparisons.wls').write_text('\n'.join(statements)+'\n')
print('completed:',len(records),'evaluations:',sum(z['evaluations'] for r in records for z in r['results']))
for x in winners:print(x['tau'],x['k'],x['raw_best'],x['baseline'],x['delta'],x['seed'],x['segments'])
