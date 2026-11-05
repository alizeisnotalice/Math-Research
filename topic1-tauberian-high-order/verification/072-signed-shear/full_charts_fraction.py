from fractions import Fraction as F
from pathlib import Path
import json
Z=(0,0,0,0)
def C(a):return {Z:F(a)} if a else {}
def V(i):
 p=[0]*4;p[i]=1;return {tuple(p):F(1)}
def add(*ps):
 r={}
 for p in ps:
  for k,v in p.items():r[k]=r.get(k,F(0))+v
 return {k:v for k,v in r.items() if v}
def sc(p,a):return {k:v*a for k,v in p.items() if v*a}
def mul(p,q):
 r={}
 for i,x in p.items():
  for j,y in q.items():
   k=tuple(a+b for a,b in zip(i,j));r[k]=r.get(k,F(0))+x*y
 return {k:v for k,v in r.items() if v}
def power(p,n):
 o=C(1)
 for _ in range(n):o=mul(o,p)
 return o
def integ(p,i,L,H):
 out={}
 for exp,c in p.items():
  n=exp[i]+1;rem=list(exp);rem[i]=0
  x=mul({tuple(rem):c/F(n)},add(power(H,n),sc(power(L,n),-1)))
  out=add(out,x)
 return out
def calc(chart,term,e,f,w):
 s,y,z,x=(V(i) for i in range(4));E=C(e)
 if chart==1:d1,d2=s,x
 elif chart==3:d1,d2=add(s,y),x
 else:d1,d2=add(s,y),s
 a=add(E,sc(d1,-1));b=add(E,sc(d2,-1))
 p=(b,a,sc(mul(a,b),1/w))[term-1]
 if chart in (1,3):
  p=integ(p,3,s,E);p=integ(p,0,C(0),sc(add(E,sc(y,-1)),e/f));p=integ(p,1,C(0),E)
 else:
  bound=sc(mul(add(E,sc(y,-1)),add(E,sc(z,-1))),1/f)
  p=integ(p,0,C(0),bound);p=integ(p,2,C(0),E);p=integ(p,1,C(0),E)
 assert all(k==Z for k in p)
 return p.get(Z,F(0))
base=Path(__file__).resolve().parent
ref=json.loads((base/'full_charts_samples.json').read_text())
params=[(F(1,16),F(1),F(1)),(F(1,32),F(3,4),F(2))]
rows=[]
for x in ref:
 e,f,w=params[x['sample']-1];actual=calc(x['chart'],x['term'],e,f,w);expected=F(x['value'].replace(' ',''))
 rows.append({'chart':x['chart'],'term':x['term'],'sample':x['sample'],'value':str(actual),'matches_wl':actual==expected})
print(json.dumps({'parameters':[[str(z) for z in p] for p in params],'cases':len(rows),'all_match':all(r['matches_wl'] for r in rows),'rows':rows},indent=2))
