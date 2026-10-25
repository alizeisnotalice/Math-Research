"""Eight new checks by confluent residue series, no DD table/interpolation solve."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from functools import lru_cache
from math import comb,factorial
import json,time,hashlib,os
import sympy as S
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent

def B(x):x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def ends(v):return Q(str(v.lower().fmpq())),Q(str(v.upper().fmpq()))
def stored(d):
 p=Q(2)**int(d['common_binary_exponent']);return int(d['lower_mantissa'])*p,int(d['upper_mantissa'])*p

def weights(nodes):
 counts=Counter(nodes);out=[]
 for x,m in sorted(counts.items()):
  coeff=[Q(1)]+[Q(0)]*(m-1)
  for y,n in counts.items():
   if x==y:continue
   b=[Q((-1)**r*comb(n+r-1,r),(x-y)**(n+r)) for r in range(m)]
   coeff=[sum(coeff[j]*b[r-j] for j in range(r+1)) for r in range(m)]
  for k in range(m):out.append((x,k,coeff[m-1-k]))
 return out

def main():
 assert not(P/'scalar_results.json').exists()
 clock=json.loads((P/'input_manifest.json').read_text());src=Path(clock['source_directory']);tic=time.monotonic()
 ctx.prec=1024
 x,u=S.symbols('x u');forms=[S.diff(u/(u*u+x*x),x,k-1)/factorial(k) for k in range(1,8)]
 @lru_cache(None)
 def jet(v,t,k):
  if k==0:return (B(v)/B(t)).atan()
  exact=S.cancel(forms[k-1].subs({x:S.Rational(v.numerator,v.denominator),u:S.Rational(t.numerator,t.denominator)}))
  return B(Q(str(exact)))
 sets=[[Q(0)]*6,list(map(Q,[-2,-1,0,0,1,2])),list(map(Q,[-1,0,0,0,0,2])),list(map(Q,[3,-2,-3,1,-1,2])),[Q(1,3)]*6,[Q(-1,7),Q(-1,7),Q(2,9),Q(2,9),Q(2,9),Q(4,5)]]
 cases=[(0,0,0),(0,0,1),(1,2,3),(2,1,1),(3,0,0),(3,0,1),(4,0,4),(5,1,5)]
 saved=json.loads((src/'scalar_native.json').read_text())['scalar_Hessian_checks'];lookup={(r['set'],r['i'],r['j']):r['hessian'] for r in saved};rows=[]
 for case in cases:
  k,i,j=case;z=sets[k];ws=weights(z+[z[i],z[j]])
  # Residue identity annihilates polynomials degree<7, and returns1 on degree7.
  assert all(sum(w*(0 if d<n else comb(d,n)*v**(d-n)) for v,n,w in ws)==int(d==7) for d in range(8))
  ans=sum((B(w)*(jet(v,Q(1,250),n)-jet(v,Q(6,125),n)) for v,n,w in ws),arb(0))*(2 if i==j else 1)
  lo,hi=ends(ans);a,b=stored(lookup[case]);assert a<=lo<=hi<=b,case
  rows.append({'set':k,'i':i,'j':j,'residue_weights':[[str(v),n,str(w)] for v,n,w in ws],'multiplier':2 if i==j else 1,'reference_lower_exact':str(lo),'reference_upper_exact':str(hi),'native_contains_reference':True})
 out={'status':'PASS_EIGHT_NEW_SCALAR_RELATIONS','algorithm':'exact confluent pole-product series weights plus symbolic rational derivatives; no author DD recurrence or Hermite interpolation solve','precision_bits':1024,'cases':rows,'pid':os.getpid(),'runtime_seconds':time.monotonic()-tic,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'full_native_calls':0,'author_anchor_calls':0}
 (P/'scalar_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','pid','runtime_seconds']}))
if __name__=='__main__':main()
