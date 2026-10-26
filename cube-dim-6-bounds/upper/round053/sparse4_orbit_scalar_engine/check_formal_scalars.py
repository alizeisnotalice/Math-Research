"""16 FORMAL checks; T and epsilon are indeterminates, no physical values."""
import json,math,time
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import sympy as S
P=Path(__file__).resolve().parent
SETS=[[Q(0)]*6,list(map(Q,[-2,-1,0,0,1,2])),list(map(Q,[-1,0,0,0,0,2])),list(map(Q,[-3,-2,-1,1,2,3])),[Q(1,3)]*6,[Q(-1,7),Q(-1,7),Q(2,9),Q(2,9),Q(2,9),Q(4,5)],list(map(Q,[0,0,0,1,1,1])),list(map(Q,[-1,-1,-1,0,0,2]))]
def leading_weights(z):
 c=[(x,k) for x,m in sorted(Counter(z).items()) for k in range(m)];V=[[Q(0) if j<k else Q(math.comb(j,k))*x**(j-k) for j in range(6)] for x,k in c]
 A=[[V[j][i] for j in range(6)]+[Q(i==5)] for i in range(6)]
 for k in range(6):
  p=next(i for i in range(k,6) if A[i][k]);A[k],A[p]=A[p],A[k];q=A[k][k];A[k]=[a/q for a in A[k]]
  for i in range(6):
   if i!=k:
    q=A[i][k];A[i]=[a-q*b for a,b in zip(A[i],A[k])]
 return c,[A[j][6] for j in range(6)]
def laurent_jet(x,k):
 v=[Q(0)]*4
 if x and k==0:v[0]=Q(1 if x>0 else -1)
 if not x and k in (1,3,5):v[(k+1)//2]=Q((-1)**((k-1)//2),k)
 return v
def dd_vec(z):
 z=sorted(z);d=[laurent_jet(x,0) for x in z]
 for k in range(1,6):d=[laurent_jet(z[j],k) if z[j]==z[j+k] else [(a-b)/(z[j+k]-z[j]) for a,b in zip(d[j+1],d[j])] for j in range(6-k)]
 return d[0]
def main():
 if(P/'formal_results.json').exists():raise SystemExit('Refuse formal replay')
 tick=time.monotonic();T,y=S.symbols('T y',positive=True);f=S.atan(y/T);jets=[S.diff(f,y,k)/S.factorial(k) for k in range(6)];rows=[]
 for case,z in enumerate(SETS):
  # H formula at symbolic common T, including genuine Hermite zero jets.
  zz=sorted(z);d=[S.atan(S.Rational(x.numerator,x.denominator)/T) for x in zz]
  for k in range(1,6):d=[jets[k].subs(y,S.Rational(zz[j].numerator,zz[j].denominator)) if zz[j]==zz[j+k] else (d[j+1]-d[j])/S.Rational((zz[j+k]-zz[j]).numerator,(zz[j+k]-zz[j]).denominator) for j in range(6-k)]
  e=[Q(1)]+[Q(0)]*6
  for x in z:
   for k in range(6,0,-1):e[k]+=x*e[k-1]
  numerator=T**6-S.Rational(e[2].numerator,e[2].denominator)*T**4+S.Rational(e[4].numerator,e[4].denominator)*T**2-S.Rational(e[6].numerator,e[6].denominator)
  product=numerator/S.prod(T*T+S.Rational(x.numerator,x.denominator)**2 for x in z)
  residual=S.factor(S.diff(d[0],T)+product);assert residual==0
  rows.append(dict(configuration=2*case+1,kind='H symbolic derivative versus independent real resolvent product',nodes=[str(x) for x in z],residual='0',T_substituted=False))
  v=dd_vec(z);conditions,w=leading_weights(z);ind=[sum(q*laurent_jet(x,k)[j] for q,(x,k) in zip(w,conditions)) for j in range(4)]
  assert v==ind==dd_vec([-x for x in z])
  rows.append(dict(configuration=2*case+2,kind='FORMAL G Laurent constant and all negative powers versus Hermite polynomial weights',nodes=[str(x) for x in z],pi_over2_coefficient=str(v[0]),epsilon_negative_coefficients={str(-k):str(v[(k+1)//2]) for k in [1,3,5]},weights_exact=[str(x) for x in w],conditions=[[str(x),k] for x,k in conditions],negative_pair_parity=True,physical_G_eligibility_claim=False))
 assert len(rows)==16
 out=dict(status='PASS_FORMAL_ONLY_G_PHYSICAL_INTERPRETATION_PENDING_B052',configurations=rows,positive_T_evaluations=0,full_kernel_calls=0,sympy=S.__version__,seconds=time.monotonic()-tick)
 (P/'formal_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],formal_checks=16,seconds=out['seconds'],physical_values=0)))
if __name__=='__main__':main()
