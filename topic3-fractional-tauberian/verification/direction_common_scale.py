"""Exact layer verification and finite common-scale capacity oracle.
Numerical optimization is regression only; proof is in the paper.
"""
from fractions import Fraction as F
from pathlib import Path
import json, platform, time, sys
import numpy as np
import scipy
import mpmath as mp
from witness_capacity import atoms
from eight_seeded_search import capacity
from eight_witness import Interval as I
ROOT=Path(__file__).resolve().parents[1]/'results'

def layers(ell,r,b):
 a=1-b
 return [(ell+2*b*r)*(ell+4*a*r-2*b*r),ell*ell+4*a*ell*r+(-8*a*b+12*b*b)*r*r,ell*ell+4*b*ell*r-12*b*b*r*r,ell*ell+4*b*ell*r+(-16*a*b+4*b*b)*r*r,(ell-2*a*r)*(ell+2*a*r-4*b*r),(ell-2*a*r)*(ell+2*a*r-4*b*r),(ell-2*a*r)**2,(ell-2*a*r)**2]

def total(ell,r,b):
 v=layers(ell,r,b)
 return min(k+sum(v[k:]) for k in range(9))

def oracle(alpha,r,b):
 # Vk(ell) is quadratic and increasing for ell>2ar.
 ls=[alpha,1/alpha]
 for k in range(8):
  c=layers(mp.mpf(0),r,b)[k];z=layers(mp.mpf(1),r,b)[k]-c-1
  root=(-z+mp.sqrt(z*z-4*(c-1)))/2
  if alpha<root<1/alpha:ls.append(root)
 vals=[total(e,r,b)-8*alpha*e for e in ls]
 j=max(range(len(ls)),key=lambda j:vals[j])
 return vals[j],ls[j]

def run():
 start=time.monotonic();exact=0;lp_err=0
 for b in [F(j,20) for j in range(1,10)]:
  a=1-b
  for r in [F(1,100),F(1,20),F(1,5)]:
   for ell in [F(9,10),F(1),F(11,10)]:
    centers=[(sx*a*r,sy*b*r) for sx in [-1,1] for sy in [-1,1]]+[(sx*b*r,sy*a*r) for sx in [-1,1] for sy in [-1,1]]
    boxes=[(x-ell/2,x+ell/2,y-ell/2,y+ell/2) for x,y in centers];cap=atoms(boxes)
    got=[sum(v for m,v in cap.items() if m.bit_count()>=k) for k in range(1,9)]
    assert got==layers(ell,r,b);exact+=1
    alpha=F(4,5);num=capacity(np.array(centers,float),np.repeat(float(ell),8),float(alpha));err=abs(num-float(total(ell,r,b)/8-alpha*ell));lp_err=max(lp_err,err);assert err<2e-12
 mp.mp.dps=80;rows=[]
 for bs in ['.05','.125','.25','.26393202250021030359','.375','.45']:
  b=mp.mpf(bs);a=1-b;gamma=(1-4*b*b)/2
  for ts in ['.001','.0001']:
   t=mp.mpf(ts);al=1-t;r0=(1-al*al)/(al+mp.sqrt(al*al-4*a*b*(1-al*al)))
   for d in [mp.mpf('-.01'),mp.mpf('.01')]:
    sig=gamma+d;r=r0*(1-sig*t);margin,e=oracle(al,r,b);assert (margin>0)==(d>0)
    e4=-2*b*r+mp.sqrt(1+16*a*b*r*r);e5=2*b*r+mp.sqrt(1+4*(a-b)**2*r*r)
    two=max(total(ee,r,b)-8*al*ee for ee in [e4,e5]);assert abs(two-margin)<mp.mpf('1e-70')
    grid=max(total(al+(1/al-al)*j/1000,r,b)-8*al*(al+(1/al-al)*j/1000) for j in range(1001));assert grid<=margin+mp.mpf('1e-70')
    rows.append(dict(b=bs,t=ts,sigma=mp.nstr(sig,45),margin=mp.nstr(margin,45),best_ell=mp.nstr(e,45),grid_gap=mp.nstr(margin-grid,20)))
 cert=[]
 for ts in ['1/1000','1/10000']:
  t=I(F(ts));al=1-t;b=F(3,8);a=1-b;d=a-b
  r=t+2*b*t*t+F(39,64)*t*t*t;L=(1+4*d*d*r*r).sqrt();z=L-2*d*r;e=L+2*b*r;vv=layers(e,r,b)
  active=1-d*r*z-al*e
  rect=t/(1-2*b*t);radial=r-rect
  for value in [active,radial,e-al,1-al*e,e-2*a*r,vv[3]-1,1-vv[6]]:assert value.lo>0
  cert.append(dict(t=ts,b='3/8',radius_formula='t+(3/4)t^2+(39/64)t^3',margin=active.encode(),rectangle_excess=radial.encode(),ell=e.encode(),r=r.encode(),L=L.encode(),z=z.encode(),sqrt_grid=80))
 from direction_common_coefficients import check
 formal=check()
 out=dict(formal_series_crosscheck=formal,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,precision=80,exact_layers=exact,LP_crosschecks=exact,max_LP_error=lp_err,oracle_samples=rows,exact_interval_certificates=cert,seconds=time.monotonic()-start,scope='same scale for all eight witnesses; no arbitrary independent-scale upper bound')
 (ROOT/'direction_common_scale_python.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k not in ['oracle_samples','exact_interval_certificates']})
if __name__=='__main__':run()
