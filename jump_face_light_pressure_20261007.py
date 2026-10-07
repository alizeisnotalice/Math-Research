"""Exact posterior-algebra guard; no certification of actual FIRST geometry."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
rows=[]
for n in (4,16,64):
 d=F(1,n)
 profiles=[('uniform',[d]*n,[d]*n),('swapped',[d,1-d]+[F(0)]*(n-2),[1-d,d]+[F(0)]*(n-2)),('opposed', [F(i+1,n*(n+1)//2) for i in range(n)], [F(n-i,n*(n+1)//2) for i in range(n)])]
 for label,p,h in profiles:
  for masked in (False,True):
   def gate(i,j): return F((i+2*j)%5,4) if masked else F(1)
   diag=sum((p[i]*h[i]*gate(i,i) for i in range(n) if p[i]<=d or h[i]<=d),F(0))
   off=sum((p[i]*h[j]*gate(i,j) for i in range(n) for j in range(n) if i!=j and (p[i]<=d or h[j]<=d)),F(0))
   assert diag<=2*d
   assert sum(p)==sum(h)==1
   if label=='uniform' and not masked: assert off==1-d and off>2*d
   if label=='swapped' and not masked: assert diag==2*d*(1-d)
   rows.append(dict(n=n,profile=label,masked=masked,delta=str(d),diagonal_light=str(diag),offdiagonal_light=str(off),diagonal_bound=str(2*d)))
path=Path(__file__)
out=path.with_name(path.stem+'_results.json')
assert not out.exists(), 'Preserve prior result; do not overwrite.'
out.write_text(json.dumps(dict(status='PASS_EXACT_POSTERIOR_ALGEBRA_ONLY',rows=rows,scope='Subprobability posterior products only; gates in [0,1]. Does not certify FIRST, original kernels, or short-shell budget.',script_sha256=hashlib.sha256(path.read_bytes()).hexdigest()),indent=2)+'\n')
print(len(rows),'exact configurations passed')
