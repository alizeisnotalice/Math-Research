from pathlib import Path
from fractions import Fraction as F
from math import comb,isqrt
import json,hashlib
D=Path(__file__).parent
src=D/'actual_future_capacity_bridge_exact_guard_results_20261007.json'
data=json.loads(src.read_text()); checks=0

def check(v):
 global checks
 checks+=1
 assert v

def digest(a,b):
 chunks=[]
 for x in (a,b):
  z=abs(x).to_bytes(max(1,(abs(x).bit_length()+7)//8),'big')
  chunks.append((b'-' if x<0 else b'+')+len(z).to_bytes(8,'big')+z)
 return hashlib.sha256(b''.join(chunks)).hexdigest()

for rnd in data['rounds']:
 n=rnd['n']; ar=rnd['source_arity']; rn=isqrt(n)
 for row in rnd['models']:
  sig,g,p=map(F,(row['sigma'],row['gamma'],row['p_C']))
  k=row['certified_integer_cutoff']; U=row['u_integer_upper']; rr=row['sqrt_integer_upper']
  check((2**U)*p>=1); check(F(rr*rr)>=2*n*sig*U)
  check(k>=n*sig+rr+F(2*U,3)); check(p*F(ar-1,ar)==F(1,rn*ar))
  z=F(row['rational_tilt_z']); B=1-sig+sig*z
  pn,pd=B.numerator**n*z.denominator**k,B.denominator**n*z.numerator**k
  check(digest(pn,pd)==row['prefactor_certificate']['framed_integer_pair_sha256']); check(pn*p.denominator<=pd*p.numerator)
  # Independently sum upper tail DOWNWARD from its top coefficient.
  prob=sig*g/(1-sig+sig*g); a,b=prob.numerator,prob.denominator-prob.numerator
  den=prob.denominator**n; term=a**n; tail=gate=0
  for j in range(n,k-1,-1):
   tail+=term; gate+=(6 if j%2 else 4)*term
   if j>k:
    nt=term*j*b; dv=(n-j+1)*a
    check(nt%dv==0); term=nt//dv
  check(term==comb(n,k)*a**k*b**(n-k))
  check(digest(tail,den)==row['exact_plain_tail']['framed_integer_pair_sha256'])
  check(digest(gate,12*den)==row['exact_gate_tail']['framed_integer_pair_sha256'])
  check(0<=gate<=12*tail); check(tail*p.denominator<=den*p.numerator)
  check(rn*(ar-1)*gate<=ar*12*den)
 for row in rnd['endpoints']:
  sig,p=F(row['sigma']),F(row['p']); k=row['k']
  if sig in (0,1):
   K=n*sig; tail=int(K>=k)
  else:
   check(p==0); tail=0
  check(tail==row['tail_mass']); check(tail<=p)
check(hashlib.sha256((D/'actual_future_capacity_bridge_exact_guard_20261007.py').read_bytes()).hexdigest()==data['script_sha256'])
out={'status':'PASS_INDEPENDENT_SAVED_RECONSTRUCTION','checks':checks,'models':18,'endpoint_models':15,'method':'upper-tail integer sum descending from k=n; author code neither imported nor run','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'coefficient/capacity only, not actual FIRST or geometry'}
(D/'root_future_capacity_saved_review_20261007.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out))
