import json,time
from decimal import Decimal as Dec, getcontext
from pathlib import Path
from satellite import eval_one
getcontext().prec=80
n=3;lam=Dec(1);u=Dec(1)/40;kap=Dec(2);beta=2.5
opt=u/(1+2*u);p=lam*u**(n+2)
def f(d):return (2*d)**n+kap*p-(2*u*(1-2*d))**n
lo=Dec(0);hi=opt
assert f(lo)<0 and f(hi)>0
for _ in range(300):
 mid=(lo+hi)/2
 if f(mid)>0:hi=mid
 else:lo=mid
d=(lo+hi)/2;c=(opt-d)/u**3
res=eval_one(n,float(lam),beta,float(c),float(u),77123)
res['delta_exact_80digit']=str(d);res['budget_residual_80digit']=str(f(d));res['c_from_exact_budget']=str(c);res['scope']='floating Sobol numerical check; exact budget root 80-digit decimal; no sampling error bound'
(Path(__file__).resolve().parent/'satellite').mkdir(exist_ok=True)
(Path(__file__).resolve().parent/'satellite'/'kappa_calibration.json').write_text(json.dumps(res,indent=2));print(json.dumps({k:res[k] for k in ('gamma_estimate','main_missed_scaled','delta_exact_80digit','budget_residual_80digit','satellite_center_margin')}))
