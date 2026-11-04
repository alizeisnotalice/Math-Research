"""Exact rational replay; SciPy solver output is not trusted by this checker."""
from fractions import Fraction as F
from pathlib import Path
import json
from witness_capacity import atoms

def verify(row):
 al=F(row['alpha']);ell=F(row['ell']);a=F(row['axis_radius']);C=[(a,0),(-a,0),(0,a),(0,-a),(a/2,a/2),(-a/2,a/2),(a/2,-a/2),(-a/2,-a/2)];boxes=[(x-ell/2,x+ell/2,y-ell/2,y+ell/2) for x,y in C];c=atoms(boxes);n=8
 for i in range(n):assert sum(v for s,v in c.items() if s>>i&1)==ell*ell
 for i in range(n):
  for j in range(i):
   I=max(F(0),ell-abs(C[i][0]-C[j][0]))*max(F(0),ell-abs(C[i][1]-C[j][1]));assert sum(v for s,v in c.items() if s>>i&1 and s>>j&1)==I
 lam=list(map(F,row['lambda']));z=F(row['z']);assert min(lam)>=0 and sum(lam)==1 and z>=0
 val=z+sum(v*max(F(0),sum(lam[i] for i in range(n) if s>>i&1)-z) for s,v in c.items())-al*ell
 assert val==F(row['dual_certificate'])
 m={int(s):F(v) for s,v in row['primal_masses'].items()};assert all(0<=v<=c.get(s,0) for s,v in m.items()) and sum(m.values())<=1
 slack=min(sum(v for s,v in m.items() if s>>i&1)-al*ell for i in range(n));assert slack==F(row['primal_strict_slack'])
 assert (row['status']=='exact-infeasible' and val<=0) or (row['status']=='exact-feasible' and slack>0)
 if row['status']=='exact-feasible':
  radii=[(sum(v for s,v in m.items() if s>>i&1)-al*ell)/(2*ell) for i in range(n)];LB=[(x-r,x+r,y-r,y+r) for (x,y),r in zip(C,radii)];assert sum(atoms(LB).values())==F(row['certified_halo_lower_bound'])
if __name__=='__main__':
 root=Path(__file__).resolve().parents[1]/'results';data=json.loads((root/'witness_capacity.json').read_text())
 for r in data['rows']:verify(r)
 print({'exact_rational_certificates_verified':len(data['rows']),'counts':data['counts'],'scope':'fixed finite rational squares, strict feasibility and neighbourhood volumes'})
