"""Numerical Fourier cap, before its asymptotic proof; not interval certified."""
import json,math
from pathlib import Path
from scipy.special import i0e,gammaln
from scipy.integrate import quad
p=Path(__file__).resolve().parent
rows=[]
for n in (16,128,10000,1000000):
 for z in (.125,.5):
  for h in (.15,.01,.001):
   k=int(n*z); tau=.2; alpha=1+tau*(z-.5)
   b=abs(tau)*h/(1-abs(alpha)*h); eta=b+b*b/4
   V=n*z*(1-z)*(1-eta)
   logmass=float(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*math.log(z)+(n-k)*math.log1p(-z))
   rows.append({'n':n,'k':k,'h':h,'V_lower':V,'eta':eta,'pointwise_cap':math.exp(math.log(float(i0e(V)))-logmass)})
checks=[]
for V in (.1,1.,10.,100.):
 q=quad(lambda t:math.exp(-2*V*math.sin(t/2)**2),-math.pi,math.pi,epsabs=1e-13)[0]/(2*math.pi)
 checks.append({'V':V,'integral':q,'i0e':float(i0e(V)),'error':abs(q-i0e(V))})
r={'status':'floating exploration only; exact integral is proved separately','rows':rows,'quadrature_checks':checks}
(p/'CAP_RESULTS.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'cases':len(rows),'max_quadrature_error':max(x['error'] for x in checks),'last':rows[-1]},indent=2))
