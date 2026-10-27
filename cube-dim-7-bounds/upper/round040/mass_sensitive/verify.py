"""One exact finite common-f interface check, not a high-dimensional computation."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parent
a=F(16); eta=F(1,10); tau=F(1); lam=F(1)
Rs=[F(1),F(1,3),F(1,9)]
ys=[F(3,4)*rr for rr in Rs]
rs=[rr/a for rr in Rs]
Q=[(yy-rr,yy+rr) for yy,rr in zip(ys,rs)]
S=[(yy-RR,yy+RR) for yy,RR in zip(ys,Rs)]
assert all(lo<0<hi for lo,hi in S)
for i in range(3):
 for j in range(i):
  assert Q[i][1]<Q[j][0] or Q[j][1]<Q[i][0]
assert all(S[0][0]<=lo<hi<=S[0][1] for lo,hi in S)
V=[2*RR for RR in Rs]
v=[2*rr for rr in rs]
M=[(1+eta)*VV for VV in V]
p=[vv/mm for vv,mm in zip(v,M)]
U=sum(v); total_f=(1+eta)*V[0]
assert U==F(13,72)
assert total_f==F(11,5)
assert sum(p)==F(15,88)==(F(3)/a)/(1+eta)
assert sum(F(1)/a for _ in Rs)<tau
assert U/total_f==F(65,792)
assert (1+eta)/(1+1/a)>1 # all small cubes lie in true strict level set.
# Exact positive-volume incidence cells of the stars.
cuts=sorted(set(e for ss in S for e in ss))
integral=F(0); cell_loads=[]
for lo,hi in zip(cuts,cuts[1:]):
 mid=(lo+hi)/2
 load=sum(pp for pp,(sl,sh) in zip(p,S) if sl<mid<sh)
 integral+=(1+eta)*load*(hi-lo)
 cell_loads.append(load)
assert integral==U
assert max(cell_loads)==sum(p)
# The entire admissible budget interval brackets p. Check a hybrid point too.
theta=[F(0),F(1,2),F(1)]
h=[(1+eta*th)*VV for th,VV in zip(theta,V)]
ph=[vv/hh for vv,hh in zip(v,h)]
assert all(F(1)/(a*(1+eta))<=pp<=F(1)/a for pp in ph)
assert sum(p)<=sum(ph)<=F(3)/a
# Exact disjoint planting ratio; no floating point limit is used.
eps=F(1,5); fb=F(7); ub=F(3)
ratio=(ub+eps*U)/(fb+eps*total_f)
weighted=(fb*(ub/fb)+eps*total_f*(U/total_f))/(fb+eps*total_f)
assert ratio==weighted
receipt={
 'status':'PASS: one exact finite common-f/CF interface and planting identity',
 'dimension':1,'a':str(a),'eta':str(eta),'selected_volume':str(U),
 'input_mass':str(total_f),'CF_endpoint':str(sum(p)),
 'true_selected_ratio':str(U/total_f),
 'checks':['disjoint actual small increments','common open star intersection',
           'one common strict L1 witness','all allowed weights remain eligible',
           'D empty exact CF maximum','f-weighted equality integral fG=lambda U',
           'admissible-budget interval bracket','exact planting weighted ratio'],
 'scope':'interface checks only; high-dimensional saturation is analytic using R037',
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                  for name in ('PROOF.md','verify.py')},
 'no_production_search':True,'no_MATLAB':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
