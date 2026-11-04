import json, math
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent

def verify(E,alpha):
 ends=sum(([a,b] for a,b in E),[]); ell=sum(b-a for a,b in E); lo=ends[0]-ell/(2*alpha); hi=ends[-1]+ell/(2*alpha)
 pieces=list(E)
 for e in ends:
  breaks=sorted({lo,hi,*((e+f)/2 for f in ends if lo<(e+f)/2<hi)})
  def G(x):
   l=min(e,2*x-e); r=max(e,2*x-e)
   return sum(max(Q(0),min(r,b)-max(l,a)) for a,b in E)-2*alpha*abs(x-e)
  for a,b in zip(breaks,breaks[1:]):
   ga,gb=G(a),G(b)
   if ga>0 and gb>0: pieces.append((a,b))
   elif ga>0: pieces.append((a,a+(b-a)*ga/(ga-gb)))
   elif gb>0: pieces.append((a+(b-a)*(-ga)/(gb-ga),b))
 out=[]
 for a,b in sorted(pieces):
  if b<=a:continue
  if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
  else: out.append((a,b))
 return out

records=json.loads((P/'candidates.json').read_text()); lines=[]
for r in records:
 E=[tuple(map(Q,x)) for x in r['intervals']]; H=verify(E,1-Q(1,2**r['k'])); expected=[tuple(map(Q,x)) for x in r['halo']]
 assert H==expected, (r['seed'],H,expected)
for tau in [.25,.5,.75]:
 for k in [4,16]:
  r=next(x for x in records if x['tau']==tau and x['k']==k and x['segments']==3)
  E=[tuple(map(Q,x)) for x in r['intervals']]; H=verify(E,1-Q(1,2**k)); q=Q(str(tau)); t=Q(1,2**k)
  def wlmass(S):return '+'.join(f'(Sign[{b}] Abs[{b}]^({q})-Sign[{a}] Abs[{a}]^({q}))' for a,b in S)
  lines.append('N[((('+wlmass(H)+')/('+wlmass(E)+'))-1)/('+str(t)+')^('+str(q)+'),60]')
(P/'wolfram_check.wls').write_text('Print /@ {'+','.join(lines)+'};\n')
print('Independent exact verifier passed',len(records),'candidates')
# Fixed geometry: hold outer endpoints fixed and scale every hole endpoint about 0 by t/t0.
rows=[]
for tau in [.25,.5,.75]:
 for n in [3,4]:
  r=next(x for x in records if x['tau']==tau and x['k']==16 and x['segments']==n)
  E=[tuple(map(Q,x)) for x in r['intervals']]; ends=sum(([a,b] for a,b in E),[])
  for k in [16,20,24,28]:
   t=Q(1,2**k); scale=t*2**16; ee=[ends[0]]+[v*scale for v in ends[1:-1]]+[ends[-1]]; ee=list(zip(ee[::2],ee[1::2])); h=verify(ee,1-t)
   def w(S):return sum(math.copysign(abs(float(b))**tau,float(b))-math.copysign(abs(float(a))**tau,float(a)) for a,b in S)
   rows.append({'tau':tau,'segments':n,'k':k,'score':(w(h)/w(ee)-1)/float(t)**tau})
(P/'fixed_shapes.json').write_text(json.dumps(rows,indent=2))
