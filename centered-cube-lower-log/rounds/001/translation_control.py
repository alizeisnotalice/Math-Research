"""Numerical checks of periodic averaging identity, before analytic estimates.
Run with NumPy; no statistical lower-bound certificate is claimed.
"""
import numpy as np, json
from fractions import Fraction as F
from pathlib import Path
p=Path(__file__).resolve().parent
rng=np.random.default_rng(2026102803)
x=rng.random(20000); L=rng.uniform(10/9,11/9,len(x))
def primitive(z):
    return z+.4*np.sin(2*np.pi*z)/(2*np.pi)+.1*np.sin(4*np.pi*z)/(4*np.pi)
def avg(x,L):return (primitive(x+L/2)-primitive(x-L/2))/L
rows=[]
for k in [1,2,3,8,32]:
    lhs=avg(x,L+k); rhs=1+L/(L+k)*(avg(x+k/2,L)-1)
    rows.append({'k':k,'max_abs_identity_error':float(np.max(abs(lhs-rhs)))})
pointwise=[]
for eps in [.1,.01,.001]:
 for k in [2,8,32,1024]:
    a=(11/9)/(11/9+k); d=np.array([-eps,eps])
    pointwise.append(float(np.max(np.log1p(a*d)/a-d)))
checks={'identity':rows,'max_log_mgf_pointwise_difference':max(pointwise),
'k1_decay':str(F(5,8)-F(45,32)*F(11,20)**2),
'k2_sum_decay':str(F(9,11)-F(1,2)),
'polynomial_decay_after_n_grid':str(16*F(511,2560)-1)}
checks['extra_events']={}
for n in [32,128,512,2048]:
 a=np.load(p/f'second_iso_n{n}_J16.npy');b=np.load(p/f'second_iso_n{n}_J64.npy')
 checks['extra_events'][str(n)]={'positive_gains':int(np.sum(b>a+1e-10)),'max_log_gain':float(np.max(b-a))}
(p/'translation_control.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
