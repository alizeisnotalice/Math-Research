from pathlib import Path
from fractions import Fraction as F
import numpy as np,json
P=Path(__file__).resolve().parent
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
p=.5;v=64.;rng=np.random.default_rng(2026102916);N=20000;x=rng.random(N);L=rng.uniform(10/9,11/9,N);c=np.sqrt(2*b/3)
# Positive background is certified by the existing D bound (20/3)H16<23.
assert F(20,3)*sum((F(1,j) for j in range(1,17)),F(0))<23
bg_lower=F(1,2)-F(23,2*K);assert bg_lower>0

def BB(length):return (np.sin(np.pi*length[:,None]*m)*np.cos(2*np.pi*x[:,None]*m))@c

def atomic(length):
 up=np.mod(x+length/2,1);lo=np.mod(x-length/2,1)
 H=np.searchsorted(atoms,up,side='right')-K*up-np.searchsorted(atoms,lo,side='right')+K*lo
 return 1+p*(H+BB(length))/(K*length)
rows=[]
for w in [4,16,64,256,1024]:
 delta=1/(w*K*v)
 smooth=1+p*(prim(x+L/2)-prim(x-L/2)+BB(L))/(K*L)
 lower=(L-delta)/L*atomic(L-delta);upper=(L+delta)/L*atomic(L+delta)
 violations=max(float(np.max(lower-smooth)),float(np.max(smooth-upper)));assert violations<1e-12
 rows.append(dict(w=w,delta=delta,N=N,max_lower_violation=float(np.max(lower-smooth)),max_upper_violation=float(np.max(smooth-upper)),changed_from_atomic=int(np.count_nonzero(abs(smooth-atomic(L))>1e-12))))
scale=[]
for n in [10**8,10**16,10**40]:
 degree=1
 while 9*.25*n/(29*(2**(degree+1)-1)**2)>=16*np.log(float(n)):degree+=1
 KK=2**degree-1
 vv=F(9)*F(1,4)*n/(29*KK*KK);dstar=1/(n*KK*vv);exponent=F(n)*dstar/F(10,9)
 scale.append(dict(n=str(n),K=KK,v=str(vv),n_delta_star=str(F(n)*dstar),padding_log_factor=str(exponent),padding_factor_numeric=float(np.exp(float(exponent)))))
for row in scale:
 exponent=F(row['padding_log_factor'])
 row['two_transfer_log_factor']=str(2*exponent)
 row['two_transfer_factor_numeric']=float(np.exp(float(2*exponent)))
# Exact single-atom checks, including contact at inner/outer closed endpoints.
endpoint_checks=[]
ll=F(4);dd=F(1)
for atom in [F(-3),F(-5,2),F(-2),F(-3,2),F(0),F(3,2),F(2),F(5,2),F(3)]:
 overlap=max(F(0),min(ll/2,atom+dd/2)-max(-ll/2,atom-dd/2))/dd
 lower=int(abs(atom)<=(ll-dd)/2);upper=int(abs(atom)<=(ll+dd)/2)
 assert lower<=overlap<=upper
 endpoint_checks.append(dict(atom=str(atom),lower=lower,overlap=str(overlap),upper=upper))
out=dict(background_lower=str(bg_lower),rows=rows,scaling=scale,exact_endpoint_checks=endpoint_checks)
(P/'padding_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
