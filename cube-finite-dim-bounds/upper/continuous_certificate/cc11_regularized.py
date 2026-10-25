"""Arb quadrature for H_eps; handles every projection degeneracy without branching."""
from pathlib import Path
from fractions import Fraction as F
import runpy,hashlib,json,time
from flint import arb,acb,ctx
ctx.prec=128
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/cb11_frame_green_locked.py')
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='ffbacd4f1a466f8f20d11dc59c458be2980b45be964751b8ea1a980cbbf1d849'
ns=runpy.run_path(str(SOURCE))
def q(v):
 v=F(v);return arb(v.numerator)/v.denominator
def integral_box(y,z,halfwidth='0',eps='1/50',tol='1e-8'):
 Y=q(y)+arb(0,q(halfwidth));Z=q(z)+arb(0,q(halfwidth));X=(arb(1),Y,Z);data=[]
 for cone,ks,w in ns['_prepared'](F(1,3)):
  bs=[sum((r[i]*X[i] for i in range(3)),arb(0))/q(k) for r,k in zip(cone,ks)]
  data.append((bs,q(w)))
 calls=0
 def f(u,analytic):
  nonlocal calls
  calls+=1;res=acb(0)
  for bs,w in data:
   A=acb(1);B=acb(1)
   for b in bs:A*=1-acb(0,b)*u;B*=1+acb(0,b)*u
   res+=w*u*(1/A+1/B)/2
  return res/(8*arb.pi()**3)
 val=acb.integral(f,0,1/q(eps),abs_tol=q(tol),rel_tol=q(tol),eval_limit=20000,depth_limit=20)
 return val,calls
if __name__=='__main__':
 start=time.time();data=[]
 for y,z,h in [('1','2272387/5000000','0'),('1','0','0'),('1','1','0'),('1/2','1/2','0'),('0','0','0')]:
  v,c=integral_box(y,z,h);row={'y':y,'z':z,'halfwidth':h,'enclosure':str(v),'calls':c,'above_target':bool(v.real>q('13/2500'))};data.append(row);print(row,flush=True)
 Path(__file__).with_name('regularized_checks.json').write_text(json.dumps({'epsilon':'1/50','precision_bits':ctx.prec,'quadrature_absolute_and_relative_tolerance':'1e-8','rows':data,'seconds':time.time()-start,'scope':'finite points only; not a minimum certificate'},indent=2)+'\n')

def point_divided_difference(y,z,eps='1/50'):
 """Exact-rational point classification; no box or transcendental-node promise."""
 X=(F(1),F(y),F(z));E=q(eps);total=arb(0)
 def f(a):return (E*E+q(a)**2).log()/2
 def df(a):return q(a)/(E*E+q(a)**2)
 for cone,ks,w in ns['_prepared'](F(1,3)):
  bs=[sum((F(r[i])*X[i] for i in range(3)),F(0))/k for r,k in zip(cone,ks)]
  if len(set(bs))==1:
   b=q(bs[0]);dd=(E*E-b*b)/(2*(E*E+b*b)**2)
  elif len(set(bs))==2:
   a=next(b for b in bs if bs.count(b)==2);b=next(b for b in bs if b!=a)
   dd=(f(b)-f(a)-df(a)*q(b-a))/q(b-a)**2
  else:
   dd=sum((f(b)/q((b-bs[(j+1)%3])*(b-bs[(j+2)%3])) for j,b in enumerate(bs)),arb(0))
  total+=q(w)*dd
 return total/(8*arb.pi()**3)
