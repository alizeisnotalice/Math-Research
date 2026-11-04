import json
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
records=json.loads((P/'candidates.json').read_text())
lines=[]
for tau in [.25,.5,.75]:
 r=next(x for x in records if x['tau']==tau and x['k']==4 and x['segments']==3)
 E=[tuple(map(Q,x)) for x in r['intervals']]; H=[tuple(map(Q,x)) for x in r['halo']]; q=Q(str(tau)); t=Q(1,16)
 def mass(S):return '+'.join(f'(Sign[{b}] Abs[{b}]^({q})-Sign[{a}] Abs[{a}]^({q}))' for a,b in S)
 candidate='((('+mass(H)+')/('+mass(E)+'))-1)/('+str(t)+')^('+str(q)+')'
 base='1/((1+'+str(t)+')^('+str(q)+')-('+str(t)+')^('+str(q)+'))'
 lines.append('Print["tau='+str(tau)+' ", FullSimplify[('+candidate+') > ('+base+')]];')
(P/'algebraic_check.wls').write_text('\n'.join(lines)+'\n')
# Known explicit symmetric single-hole construction is always feasible.
envelope=[]
for tau in [.25,.5,.75]:
 for k in [4,8,12,16]:
  rr=[x for x in records if x['tau']==tau and x['k']==k]
  b=rr[0]['baseline_score']; raw=max(x['rational_score'] for x in rr)
  envelope.append({'tau':tau,'k':k,'raw_DE_3':rr[0]['rational_score'],'raw_DE_4':rr[1]['rational_score'],'single_hole_baseline':b,'feasible_envelope':max(b,raw),'attained_by':'raw_DE_3' if rr[0]['rational_score']>max(b,rr[1]['rational_score']) else ('raw_DE_4' if rr[1]['rational_score']>b else 'single_hole_baseline')})
(P/'envelope.json').write_text(json.dumps(envelope,indent=2))
