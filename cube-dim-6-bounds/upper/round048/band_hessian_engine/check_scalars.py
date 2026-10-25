"""Sixteen scalar Hessian checks by independent Hermite interpolation."""
import sys
sys.dont_write_bytecode=True
import json,hashlib,math,time
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent
SETS=[[Q(0)]*6,list(map(Q,[-2,-1,0,0,1,2])),list(map(Q,[-1,0,0,0,0,2])),list(map(Q,[3,-2,-3,1,-1,2])),[Q(1,3)]*6,[Q(-1,7),Q(-1,7),Q(2,9),Q(2,9),Q(2,9),Q(4,5)]]
CASES=[(0,0,0),(0,0,1),(1,2,2),(1,2,3),(1,0,2),(2,1,1),(2,1,2),(2,0,5),(3,0,0),(3,0,1),(3,3,5),(4,0,0),(4,0,4),(5,0,1),(5,2,2),(5,1,5)]
def B(q):q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def bounds(d):
 e=Q(2)**int(d['common_binary_exponent']);return int(d['lower_mantissa'])*e,int(d['upper_mantissa'])*e
def ends(v):
 l,e=v.lower().man_exp();u,f=v.upper().man_exp();e,f=int(e),int(f);g=min(e,f)
 return dict(lower_mantissa=str(int(l)*2**(e-g)),upper_mantissa=str(int(u)*2**(f-g)),common_binary_exponent=str(g))
def cmul(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def cpow(x,k):
 z=(Q(1),Q(0))
 for _ in range(k):z=cmul(z,x)
 return z
def jet(x,u,k):
 if not k:return (B(x)/B(u)).atan()
 den=u*u+x*x;q=(u/den,x/den)
 return B(cmul(cpow((Q(0),Q(1)),k-1),cpow(q,k))[0]/k)
def hermite_weights(nodes):
 conditions=[(x,k) for x,m in sorted(Counter(nodes).items()) for k in range(m)]
 assert len(conditions)==8
 V=[[Q(0) if d<k else Q(math.comb(d,k))*x**(d-k) for d in range(8)] for x,k in conditions]
 # Solve V^T w=e_7 using exact Fraction elimination; no DD recurrence.
 A=[[V[j][i] for j in range(8)]+[Q(i==7)] for i in range(8)]
 for k in range(8):
  pivot=next(j for j in range(k,8) if A[j][k]);A[k],A[pivot]=A[pivot],A[k]
  q=A[k][k];A[k]=[x/q for x in A[k]]
  for j in range(8):
   if j!=k:
    q=A[j][k];A[j]=[a-q*b for a,b in zip(A[j],A[k])]
 w=[A[j][8] for j in range(8)]
 assert all(sum(w[j]*V[j][d] for j in range(8))==Q(d==7) for d in range(8))
 return conditions,w
def main():
 if (P/'scalar_validation.json').exists():raise SystemExit('Refuse scalar replay')
 t=time.monotonic();ctx.prec=1024;n=json.loads((P/'scalar_native.json').read_text());assert n['parser_selftests_passed'] and len(n['scalar_Hessian_checks'])==16
 rows=[];a,T=Q(1,250),Q(6,125)
 for (s,i,j),v in zip(CASES,n['scalar_Hessian_checks']):
  assert (v['set'],v['i'],v['j'])==(s,i,j)
  z=SETS[s];cond,w=hermite_weights(z+[z[i],z[j]])
  ref=sum((B(q)*(jet(x,a,k)-jet(x,T,k)) for (x,k),q in zip(cond,w)),arb(0))*(2 if i==j else 1)
  rd=ends(ref);lo,hi=bounds(v['hessian']);rl,rh=bounds(rd);assert lo<=rl<=rh<=hi,(s,i,j)
  rows.append(dict(set=s,i=i,j=j,conditions=[[str(x),k] for x,k in cond],polynomial_leading_weights_exact=[str(x) for x in w],diagonal_multiplier=2 if i==j else 1,reference_dyadic=rd,pass_containment=True))
 out=dict(status='PASS',scalar_Hessian_checks=16,reference_precision_bits=1024,algorithm='exact degree7 Hermite interpolation weights; exact complex rational normalized atan jets plus Arb atan values; no finite differences',rows=rows,elapsed_seconds=time.monotonic()-t,native_sha256=hashlib.sha256((P/'scalar_native.json').read_bytes()).hexdigest())
 (P/'scalar_validation.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 16 scalar second derivatives by independent Hermite polynomial method')
if __name__=='__main__':main()
