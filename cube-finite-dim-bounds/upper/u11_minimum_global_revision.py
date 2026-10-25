"""Complete the SAME star revision: 20 true-vertex terms, 3 physical sign cases."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from itertools import combinations,product
from math import gcd,lcm
from functools import reduce
import hashlib,json,time
from flint import fmpq_mpoly_ctx,fmpq
ctx=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');xs=ctx.gens()
ROOT=Path(__file__).parent
stars=json.loads((ROOT/'u11_minimum_star_results.json').read_text())['stars']
def pol(dic):return ctx.from_dict({tuple(map(int,e.split(','))):fmpq(c) for e,c in dic.items()})
def canon(z):
 den=lcm(*(q.denominator for q in z));nn=[int(q*den) for q in z];div=reduce(gcd,(abs(x) for x in nn if x));pp=tuple(x//div for x in nn);sc=Q(div,den)
 if next(x for x in pp if x)<0:pp=tuple(-x for x in pp);sc=-sc
 return pp,sc

def term(star,ind,sg,kind,ray):
 n=pol(star['numerator_dict']);d=pol(star['denominator_dict']);const,fs=d.factor();co=-Q(1,16)/Q(str(const));ff=Counter()
 subs=[int(s)*xs[i] for i,s in zip(ind,sg)]
 for f,m in fs:
  ft=f.compose(*subs);vec=[]
  td=ft.to_dict()
  for i in range(4):e=[0]*4;e[i]=1;vec.append(Q(str(td.get(tuple(e),0))))
  pr,sc=canon(vec);ff[pr]+=m;co/=sc**m
 return {'n':n.compose(*subs),'co':co,'f':ff,'kind':kind,'ray':ray}
terms=[]
for i,j in combinations(range(4),2):
 rem=[k for k in range(4) if k not in (i,j)]
 for sg in (-1,1):
  ray=[0]*4;ray[i]=1;ray[j]=sg
  terms.append(term(stars[0],[i,j]+rem,[1,sg,1,1],'B',tuple(ray)))
for ss in product((-1,1),repeat=3):
 ray=(1,)+ss;terms.append(term(stars[1],list(range(4)),ray,'D',ray))
common=Counter()
for t in terms:common|=t['f']
degree=sum(common.values());print('GLOBAL common degree',degree,'factors',len(common),flush=True)
if degree>96:raise RuntimeError('registered degree cap exceeded')
lin={f:sum((ctx.constant(a)*x for a,x in zip(f,xs)),ctx.constant(0)) for f in common}
den=ctx.constant(1)
for f,m in common.items():den*=lin[f]**m
for t in terms:
 pp=t['n']*ctx.constant(fmpq(t['co'].numerator,t['co'].denominator))
 for f,m in common.items():pp*=lin[f]**(m-t['f'][f])
 t['fullnum']=pp
samples={'I':(10,3,2,1),'II':(4,3,2,Q(1,2)),'III':(5,3,2,1)}
rows=[];start=time.time()
for kind,point in samples.items():
 n=ctx.constant(0)
 for t in terms:
  v=sum(a*b for a,b in zip(point,t['ray']));sg=1 if v>0 else -1
  if t['kind']=='B':assert sg==1
  n+=sg*t['fullnum']
 print('case',kind,'built numerator',len(n.to_dict()),flush=True)
 cd=n.gcd(den);nr=n//cd;dr=den//cd
 assert nr*cd==n and dr*cd==den
 # Fix denominator positive at the chosen strict chamber point.
 ev=dr(*[fmpq(str(q)) for q in point])
 if ev<0:nr=-nr;dr=-dr
 print('case',kind,'reduced degrees',nr.total_degree(),dr.total_degree(),'terms',len(nr.to_dict()),len(dr.to_dict()),flush=True)
 print('D factors',dr.factor(),flush=True)
 rows.append({'case':kind,'numerator_dict':{','.join(map(str,e)):str(c) for e,c in nr.to_dict().items()},
              'denominator_dict':{','.join(map(str,e)):str(c) for e,c in dr.to_dict().items()},
              'denominator_factorization':str(dr.factor()),'numerator_degree':int(nr.total_degree()),'denominator_degree':int(dr.total_degree()),
              'numerator_terms':len(nr.to_dict()),'denominator_terms':len(dr.to_dict())})
out={'revision_count':1,'global_common_degree':degree,'global_common_factors':len(common),'physical_cases':rows,
     'not_yet_g_certificate':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(ROOT/'u11_minimum_global_results.json').write_text(json.dumps(out,indent=2)+'\n')
print('GLOBAL REVISION DONE',time.time()-start,flush=True)
