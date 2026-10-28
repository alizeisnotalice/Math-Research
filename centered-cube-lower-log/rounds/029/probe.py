from pathlib import Path
from fractions import Fraction as F
from math import lcm
from bisect import bisect_left
import numpy as np,json
P=Path(__file__).resolve().parent
# Reuse actual source definition only; all new integrals below use their own partition.
src=(P.parent/'017/probe.py').read_text();exec(src.split('rows=[]')[0],globals());P=Path(__file__).resolve().parent
K=1023;v=64;L=F(29,25);base_atoms=atoms.copy()
# U_j are exact binary midpoint labels: 8*U_j is an odd integer.
offset=(8*U).astype(int);assert np.all(offset==8*U)
atomq=[F(8*j+int(offset[j]),8*K) for j in range(K)]
def min_gap(c1,c2):
 den=lcm(*[z.denominator for z in c1],*[z.denominator for z in c2]);x1=[int(z*den) for z in c1];x2=sorted(int(z*den) for z in c2);xx=[z-den for z in x2]+x2+[z+den for z in x2]
 gap=min(min(abs(z-xx[bisect_left(xx,z)-1]),abs(z-xx[bisect_left(xx,z)])) for z in x1)
 return F(gap,den)
rows=[]
for w in [4,16,64,256,1024]:
 dq=F(1,w*K*v);delta=float(dq)
 def prim(t):
  y=np.mod(t,1);lo=np.searchsorted(base_atoms,y-delta/2);hi=np.searchsorted(base_atoms,y+delta/2)
  return lo+np.where(hi>lo,(y-base_atoms[np.minimum(lo,K-1)]+delta/2)/delta,0)-K*y
 for tau in [F(1,16),F(1,4),F(1),F(4)]:
  step=tau/(K*v);J=L+step;t=step/2
  c1=[(aa-L/2-t/2)%1 for aa in atomq];c2=[(aa+L/2+t/2)%1 for aa in atomq]
  crossgap=min_gap(c1,c2);assert crossgap>=t+dq
  spacing=min([atomq[j+1]-atomq[j] for j in range(K-1)]+[atomq[0]+1-atomq[-1]]);assert spacing>=t+dq
  shape=t-dq/3 if t>=dq else t*t/dq-t**3/(3*dq*dq)
  exact=2*K*shape-(K*step)**2;limit=K*step-(K*step)**2;assert 0<=exact<=limit
  edges=np.unique(np.concatenate([np.array([0.,1.])]+[np.mod(base_atoms+edge+shift,1) for ll in [float(L),float(J)] for edge in [-delta/2,delta/2] for shift in [-ll/2,ll/2]]));length=np.diff(edges);x=(edges[:-1]+edges[1:])/2
  def H(ll):return prim(x+ll/2)-prim(x-ll/2)
  def slope(ll):
   out=np.zeros(len(x))
   for sign in [-1,1]:
    y=np.mod(x+sign*ll/2,1);count=np.searchsorted(base_atoms,y+delta/2)-np.searchsorted(base_atoms,y-delta/2);out+=sign*count/delta
   return out
  value=H(float(J))-H(float(L));der=slope(float(J))-slope(float(L));numeric=np.sum(length*(value*value+der*der*length*length/12));err=abs(numeric-float(exact));assert err<1e-8
  rows.append(dict(w=w,tau=str(tau),exact_increment_squared=str(exact),limit_increment_squared=str(limit),ratio_to_limit=float(exact/limit),piecewise_integral=float(numeric),absolute_error=err,cross_gap=str(crossgap),required_gap=str(t+dq),derivative_energy_normalized=str(F(w,2)-F(1,v))))
(P/'results.json').write_text(json.dumps(rows,indent=2));print([(r['w'],r['ratio_to_limit'],r['derivative_energy_normalized']) for r in rows if r['tau']=='1'])
