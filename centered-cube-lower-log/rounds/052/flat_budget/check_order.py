"""Exact finite-atom checks of the common-order weak functional lemma."""
from fractions import Fraction as F
from pathlib import Path
import json,random
rng=random.Random(202610295202)
def W(values,masses):
 return max((v*sum(p for x,p in zip(values,masses) if x>=v) for v in set(values)),default=F(0))
rows=[]
for _ in range(64):
 n=rng.randrange(2,9);r=rng.randrange(1,6)
 masses=[F(rng.randrange(1,10)) for i in range(n)];z=sum(masses)+rng.randrange(0,5);masses=[p/z for p in masses]
 weights=[F(rng.randrange(1,10)) for i in range(r)];z=sum(weights);weights=[p/z for p in weights]
 phi=[sorted(F(rng.randrange(0,20),rng.randrange(1,6)) for i in range(n)) for a in range(r)]
 mixed=[sum(w*row[i] for w,row in zip(weights,phi)) for i in range(n)]
 lhs=W(mixed,masses);rhs=sum(w*W(row,masses) for w,row in zip(weights,phi));assert lhs<=rhs
 rows.append({'states':n,'prior_atoms':r,'lhs':str(lhs),'rhs':str(rhs)})
m=[F(1,2)]*2;f=[F(2),F(1)];g=f[::-1]
ce={'f':list(map(str,f)),'g':list(map(str,g)),'W_f':str(W(f,m)),'W_g':str(W(g,m)),'W_average':str(W([(x+y)/2 for x,y in zip(f,g)],m))}
assert ce['W_average']=='3/2'
out={'common_order_exact_cases':64,'checks':rows,'noncommon_order_counterexample':ce,'scope':'Finite-atom exact rational checks; common-order theorem is proved separately. Counterexample does not assert a true cube source realizes these two profiles.'}
(Path(__file__).resolve().parent/'ORDER_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'cases':64,'counterexample':ce},indent=2))
