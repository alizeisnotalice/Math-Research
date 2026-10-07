from pathlib import Path
from fractions import Fraction as F
import json,hashlib
import numpy as np
P=Path(__file__).resolve().parent
out=P/'ordered_continuous_bounded_entropy_results_20261007.json'
assert not out.exists()
reg=P/'ordered_continuous_bounded_entropy_registration_20261007.json'
assert reg.exists()
checks=0
for n in (8,32,128):
 for lam in (F(1,2**n),F(1),F(2**n)):
  for U in (F(0),F(1,4),F(1),F(2),F(3),F(2**n)):
   for V in (F(0),F(1,8),F(1),F(2),F(2**n)):
    u,v=lam*U,lam*V
    left=(u/(lam+u)-v/(lam+v))*(u-v)
    right=lam*(u-v)**2/((lam+u)*(lam+v))
    assert left==right and right>=0;checks+=2
  for d in (F(1,8),F(1,2),F(1),F(2)):
   for U in (1+d,2+d,F(2**n)):
    for V in (F(0),F(1,2),F(1)):
     Q=(U-V)**2/((1+U)*(1+V))
     assert Q>=d*d/(2*(2+d));checks+=1
     assert U<=2*(1+d)*(2+d)/d**2*Q;checks+=1
rounds=[]
for n,c,a0 in ((8,1.,.5),(32,.75,.9),(128,.5,.99)):
 counts=np.bincount(1+np.arange(n)%3,minlength=4)[1:]
 b=np.arange(1,4,dtype=float)**2
 b=b/np.log1p(b)-1
 theta=c*b/(1+c*b)
 def state(z):
  z=np.asarray(z);r=-np.expm1(-z)
  a=a0*np.exp(np.sum(counts[:,None]*np.log1p(-theta[:,None]*r.ravel()),axis=0))
  d=c*np.exp(-z.ravel())
  ell=np.sum(counts[:,None]*(d[None,:]*b[:,None])/(1+d[None,:]*b[:,None]),axis=0)
  return a,ell
 def entropy(a,lam):
  b0=1+lam;dd=np.sqrt(b0*b0-a*a)
  return 1-lam*np.log((b0+dd)/(2*lam))
 rows=[]
 for lam in (1/16,1.,16.):
  vals=[];Z=4/n
  for N in (512,1024,2048):
   z=np.linspace(0,Z,N+1);a,ell=state(z)
   dd=np.sqrt((1+lam)**2-a*a)
   density=lam*a*a*ell/(dd*(1+lam+dd))
   integral=(Z/N)/3*(density[0]+density[-1]+4*density[1:-1:2].sum()+2*density[2:-1:2].sum())
   exact=float(entropy(a0,lam)-entropy(a[-1],lam))
   vals.append({'N':N,'integral':float(integral),'entropy_drop':exact,'absolute_error':abs(float(integral)-exact)})
  assert vals[-1]['absolute_error']<2e-10;checks+=1
  phaseerr=0.;generatorerr=0.
  for z in np.linspace(0,Z,7):
   a,ell=state([z]);a=float(a[0]);ell=float(ell[0])
   rr=-np.expm1(-z)
   derivative_r=-a*np.sum(counts*theta/(1-rr*theta))
   generatorerr=max(generatorerr,abs(derivative_r*np.exp(-z)+a*ell))
   for Nphase in (1024,4096):
    cs=np.cos(2*np.pi*np.arange(Nphase)/Nphase)
    v=1+a*cs
    E=float(np.mean(v-lam*np.log1p(v/lam)))
    D=float(a*ell*np.mean(v/(lam+v)*cs))
    dd=np.sqrt((1+lam)**2-a*a)
    target=lam*a*a*ell/(dd*(1+lam+dd))
    phaseerr=max(phaseerr,abs(E-float(entropy(a,lam))),abs(D-target))
  assert phaseerr<1e-10 and generatorerr<1e-10;checks+=2
  rows.append({'lambda':lam,'time_integrals':vals,'max_phase_quadrature_error':phaseerr,'max_original_generator_error':generatorerr})
 rounds.append({'n':n,'c':c,'amplitude':a0,'frequency_counts':counts.tolist(),'rows':rows})
r={'status':'PASS_EXACT_SCALAR_AND_FLOAT_ORIGINAL_FLOW','checks':checks,'rounds':rounds,'scope':'Scalar checks exact Fraction; periodic original-kernel flow diagnostics float not intervals, not actual geom, no dimension fit','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest()}
out.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'status':r['status'],'checks':checks,'max_finest_integral_error':max(v['time_integrals'][-1]['absolute_error'] for q in rounds for v in q['rows']),'max_phase_error':max(v['max_phase_quadrature_error'] for q in rounds for v in q['rows'])}))
