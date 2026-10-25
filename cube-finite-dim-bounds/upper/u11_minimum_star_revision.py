"""Only permitted representation revision: two exact spectral-vertex stars."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from functools import reduce
from math import gcd,lcm,prod
import hashlib,json,time
from flint import fmpq_mpoly_ctx,fmpq
import u11_minimum_dd as dd

g=dd.g;ctx=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');vs=ctx.gens()
def canonical(z):
 den=lcm(*(q.denominator for q in z));nn=[int(q*den) for q in z]
 div=reduce(gcd,(abs(x) for x in nn if x));prim=tuple(x//div for x in nn);scale=Q(div,den)
 if next(x for x in prim if x)<0:prim=tuple(-x for x in prim);scale=-scale
 return prim,scale

def star(ray):
 rid=g.RAYS.index(ray);q=dd.QBODY;terms=[];common=Counter()
 for ids,det in g.CONES:
  if rid not in ids:continue
  fs=Counter();co=Q(det)/prod(q[i] for i in ids)
  for j in ids:
   if j==rid:continue
   z=tuple(Q(ray[k])/q[rid]-Q(g.RAYS[j][k])/q[j] for k in range(4))
   fac,scale=canonical(z);fs[fac]+=1;co/=scale
  common|=fs;terms.append((co,fs))
 degree=sum(common.values());print('star',ray,'terms',len(terms),'common degree',degree,flush=True)
 if degree>32:return {'ray':ray,'status':'stopped common degree>32','degree':degree}
 lin={f:sum((ctx.constant(int(a))*v for a,v in zip(f,vs)),ctx.constant(0)) for f in common}
 den=ctx.constant(1)
 for f,k in common.items():den*=lin[f]**k
 num=ctx.constant(0)
 for co,fs in terms:
  term=ctx.constant(fmpq(co.numerator,co.denominator))
  for f,k in common.items():term*=lin[f]**(k-fs[f])
  num+=term
 print('built exact numerator',flush=True)
 cancel=num.gcd(den);nred=num//cancel;dred=den//cancel
 assert nred*cancel==num and dred*cancel==den
 nf=nred.factor();df=dred.factor()
 print('N factor',nf,'D factor',df,flush=True)
 return {'ray':ray,'status':'exact star cancellation PASS','terms':len(terms),'common_degree':degree,
         'numerator':str(nred),'denominator':str(dred),'numerator_factorization':str(nf),'denominator_factorization':str(df),
         'numerator_dict':{','.join(map(str,e)):str(c) for e,c in nred.to_dict().items()},
         'denominator_dict':{','.join(map(str,e)):str(c) for e,c in dred.to_dict().items()},
         'exact_cancellation_identity':'num*dred==nred*den','gcd_degree':int(cancel.total_degree())}
start=time.time();rows=[star((1,1,0,0)),star((1,1,1,1))]
out={'revision_count':1,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'stars':rows,'not_a_g_certificate':True}
Path(__file__).with_name('u11_minimum_star_results.json').write_text(json.dumps(out,indent=2)+'\n')
print('REVISION DONE',time.time()-start,flush=True)
