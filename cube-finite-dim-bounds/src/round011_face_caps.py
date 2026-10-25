"""Finite exact frame enumeration and numerical display of proved cap budgets."""
from itertools import product,combinations,permutations
from fractions import Fraction as F
from pathlib import Path
import math,json,hashlib
import sympy as s
from flint import arb,ctx
ctx.prec=192
def ball(q):
 q=F(q);return arb(q.numerator)/q.denominator
def power(x,a,b):return (x.log()*a/b).exp()

def determinant(columns):
 n=len(columns);ans=0
 for perm in permutations(range(n)):
  inv=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
  ans+=(-1)**inv*math.prod(columns[j][perm[j]] for j in range(n))
 return ans

def model(n,tau):
 body=[(1,)+v for v in product((-1,1),repeat=n-1)]
 face=[]
 for i,j in combinations(range(n),2):
  for sign in (-1,1):
   v=[0]*n;v[i]=1;v[j]=sign;face.append(tuple(v))
 vectors=body+(face if tau else []); weights=[F(1)]*len(body)+([tau]*len(face) if tau else [])
 m=len(vectors);pair_sum=F(0);hist={};basis=[]
 for i,v in enumerate(vectors):
  for j,w in enumerate(vectors):
   if i==j:continue
   best=max((abs(determinant((v,w)+rest)),rest) for rest in combinations(body,n-2))
   det,rest=best;assert det>0
   pair_sum+=weights[i]*weights[j]/det;hist[str(det)]=hist.get(str(det),0)+1
   basis.append({'i':i,'j':j,'completion':rest,'abs_det':det})
 # a=1/(2*m*sqrt(n)); all completion Cauchy rates are one.
 Cnum=4*m*m*n*pair_sum;C=float(Cnum)/math.pi**n
 cmin=float(min(weights));D=1/((4 if n==3 else 16)*math.pi**n*(n-1))
 K=4*m*m*n*n*cmin*cmin
 cap=[]
 for eps in (1e-6,1e-8):
  if n==3:
   A=D+C/(2*cmin*cmin)*math.log1p(K);B=C/(cmin*cmin);L=math.log(1/eps)
   value=math.sqrt(2)*math.pi*eps**2*(A**1.5+B**1.5*(L*L+L+.5))
  else:
   B=C*math.pi*m*n/cmin
   value=2**(1/3)*4*math.pi*(D**(4/3)*eps**3/3+B**(4/3)*3*eps**(5/3)/5)
  er=ball(F(str(eps)));cr=ball(Cnum)/arb.pi()**n;cm=ball(min(weights));dr=1/(arb((4 if n==3 else 16)*(n-1))*arb.pi()**n)
  if n==3:
   ar=dr+cr/(2*cm**2)*(1+4*m*m*n*n*cm**2).log();br=cr/cm**2;lr=(1/er).log()
   cert=arb(2).sqrt()*arb.pi()*er**2*(power(ar,3,2)+power(br,3,2)*(lr**2+lr+ball(F(1,2))))
  else:
   br=cr*arb.pi()*m*n/cm
   cert=power(arb(2),1,3)*4*arb.pi()*(power(dr,4,3)*er**3/3+ball(F(3,5))*power(br,4,3)*power(er,5,3))
  cap.append({'epsilon':eps,'epsilon_exact':str(F(str(eps))),'analytic_budget_interval':[str(cert.lower().fmpq()),str(cert.upper().fmpq())],'one_face_one_cap_upper_display':value,'scope':'floating display of analytic upper expression; no rounded strict certificate'})
 return {'n':n,'tau':str(tau),'directions':vectors,'m':m,'ordered_pairs':len(basis),'det_histogram':hist,'pair_sum_exact':str(pair_sum),'C_numerator_over_pi_n':str(Cnum),'cmin_exact':str(min(weights)),'tail_D':'1/('+str((4 if n==3 else 16)*(n-1))+'*pi^'+str(n)+')','cap_displays':cap,'basis_certificate':basis}

def main():
 r=s.symbols('r',positive=True)
 primitive=r*r/2*(s.log(1/r)**2+s.log(1/r)+s.Rational(1,2))
 assert s.simplify(s.diff(primitive,r)-r*s.log(1/r)**2)==0
 assert s.integrate(r**s.Rational(2,3),r)==s.Rational(3,5)*r**s.Rational(5,3)
 cases=[model(3,F(t)) for t in ['1/10','1/3','1','3']]+[model(4,F(0))]
 out={'status':'primary analytic lemma with exact finite determinant checks; no new bound; independent audit pending','cases':cases,'shared_libraries':{'SymPy':s.__version__,'Fraction':'stdlib'},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 target=Path(__file__).resolve().parents[1]/'results/round011_face_caps.json';target.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps([{k:v for k,v in case.items() if k not in ['basis_certificate','directions']} for case in cases],indent=2))
if __name__=='__main__':main()
