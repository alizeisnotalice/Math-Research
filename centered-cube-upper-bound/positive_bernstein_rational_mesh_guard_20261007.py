from pathlib import Path
from fractions import Fraction as F
from math import comb,isqrt
from bisect import bisect_left
import hashlib,json
P=Path(__file__).resolve().parent
REG=P/'positive_bernstein_rational_mesh_registration_20261007.json'
reg=json.loads(REG.read_text()); checks=0

def check(ok,tag):
 global checks
 checks+=1
 if not ok: raise AssertionError(tag)

def mesh(M):
 q=[F(j*j,j*j+(M-j)**2) for j in range(M+1)]
 s=[F(j*(M-j),j*j+(M-j)**2) for j in range(M+1)]
 radii=[]
 for j in range(M):
  ratio=(q[j+1]-q[j])/(s[j]+s[j+1])
  check(ratio==F(M,(M-j)*(M-j-1)+j*(j+1)),'exact neighbor rational identity')
  radii.append(ratio*ratio)
  # Exact equal-distance point, including the degenerate edge limit.
  p=(q[j]*s[j+1]+q[j+1]*s[j])/(s[j]+s[j+1])
  check(q[j]<=p<=q[j+1],'covering crossover inside interval')
  if s[j]: check((p-q[j])**2/s[j]**2==radii[-1],'left radius')
  if s[j+1]: check((q[j+1]-p)**2/s[j+1]**2==radii[-1],'right radius')
 R=F(4,M*M) if M%2==0 else F(4*M*M,(M*M-1)**2)
 check(max(radii)==R,'closed-form exact global radius')
 return q,s,R

def near(p,q):
 j=bisect_left(q,p)
 options=[x for x in q[max(0,j-1):min(len(q),j+1)] if 0<x<1]
 if not options: raise AssertionError('no interior node')
 return min(options,key=lambda x:(p-x)**2/(x*(1-x)))

def value(coeff,r,n):
 return sum((a*r**k*(1-r)**(n-k) for k,a in coeff.items()),F(0))

for M in reg['other_mesh_sizes']: mesh(M)
rows=[]
for item in reg['rounds']:
 n=item['n'];M=item['M'];before=checks
 q,s,R=mesh(M);C=1/(1-n*R)
 check(n*R==F(1,512) and C==F(512,511),'common certified amplitude factor')
 min_ratio=F(1)
 for k in range(n+1):
  p=F(k,n)
  if k in (0,n):
   check(p in q,'endpoint maximum included');continue
  x=near(p,q)
  ratio=(x/p)**k*((1-x)/(1-p))**(n-k)
  check((p-x)**2/(x*(1-x))<=R,'peak covered')
  check(ratio*C>=1,'single monomial exact maximum enclosure')
  min_ratio=min(min_ratio,ratio)
 for p in (F(1,3),F(2,5),F(1,2),F(3,5),F(2,3)):
  mean=n*p;lo=max(0,int(mean)-isqrt(n));hi=min(n,int(mean)+isqrt(n)+1)
  weights={lo:(hi-mean)/(hi-lo),hi:(mean-lo)/(hi-lo)}
  check(sum(weights.values())==1 and sum(k*v for k,v in weights.items())==mean,'prescribed posterior mean')
  coeff={k:w/(p**k*(1-p)**(n-k)) for k,w in weights.items()}
  check(value(coeff,p,n)==1,'stationary response normalization')
  derivative=sum((w*(F(k)/p-F(n-k)/(1-p)) for k,w in weights.items()),F(0))
  check(derivative==0,'exact stationary point')
  x=near(p,q)
  check(C*value(coeff,x,n)>=1,'positive two-mask stationary enclosure')
 # Exact positive coefficient subdivision, including a narrow future subinterval.
 coeff={n//4:F(1),n//2:F(7),3*n//4:F(3)}
 for a,b in ((F(1,5),F(4,5)),(F(2,5),F(3,5)),(F(7,10),F(71,100))):
  out=[F(0)]*(n+1)
  for k,w in coeff.items():
   for i in range(k+1):
    left=w*comb(k,i)*a**(k-i)*b**i
    for h in range(n-k+1):
     out[i+h]+=left*comb(n-k,h)*(1-a)**(n-k-h)*(1-b)**h
  for v in out: check(v>=0,'subdivision coefficient positivity')
  outdict=dict(enumerate(out))
  for v in (F(0),F(1,7),F(1,2),F(6,7),F(1)):
   check(value(coeff,a+(b-a)*v,n)==value(outdict,v,n),'exact affine subdivision identity')
 rows.append(dict(n=n,M=M,checks=checks-before,R=str(R),factor=str(C),min_single_peak_ratio_display=float(min_ratio)))
result=dict(status='PASS_EXACT_FRACTION',checks=checks,rounds=rows,registration_sha256=hashlib.sha256(REG.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Scalar nonnegative Bernstein continuum certificate only; no original spatial weak estimate')
(P/'positive_bernstein_rational_mesh_guard_results_20261007.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result))
