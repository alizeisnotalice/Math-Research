"""Actual periodic source macroclass posterior; resolved small sources + exact-phase high-d pilot."""
from pathlib import Path
from fractions import Fraction as F
import argparse,os,sys,time,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
INPUT_DIR=Path(os.environ.get('CUBE_INPUT_DIR',str(HERE/'inputs')))
ORACLE_DIR=Path(os.environ.get('CUBE_ORACLE_DIR',str(HERE.parent/'round046_oracle')))
sys.path.insert(0,str(ORACLE_DIR))
from source_oracle import SourceOracle,NPZLabels,GFLabels,mp,mpf
GATES={'half_0':[(F(0),F(1,2))],'half_centered':[(F(0),F(1,4)),(F(3,4),F(1))]}
LD=np.longdouble;PI=np.arccos(LD(-1))

def ld(x):
 if isinstance(x,F):return LD(x.numerator)/LD(x.denominator)
 return LD(x)

def primitive_exact(oracle,t):
 t=F(t)
 fourier=sum((c/2*mp.sin(2*mp.pi*mpf((m*t)%1)) for m,c in enumerate(oracle.coeff,1)),mp.mpf(0))
 return mpf(t)+mpf(oracle.p/oracle.K)*(mpf(oracle.primitive(t))+fourier)

def gated_exact(oracle,t,intervals):
 t=F(t);period=t.numerator//t.denominator;r=t-period
 pi=sum((primitive_exact(oracle,b)-primitive_exact(oracle,a) for a,b in intervals),mp.mpf(0))
 val=period*pi
 for a,b in intervals:
  if r>a:val+=primitive_exact(oracle,min(r,b))-primitive_exact(oracle,a)
 return val

class SmallSource:
 def __init__(self,oracle):
  self.o=oracle;self.K=oracle.K;self.B=oracle.B;self.labels=oracle.labels.labels
  self.h=ld(oracle.h);self.width=ld(oracle.delta);self.p=ld(oracle.p)
  self.co=np.array([LD(str(c)) for c in oracle.coeff]);self.atoms=(np.arange(self.K,dtype=LD)+(self.labels.astype(LD)+LD('.5'))/self.B)/self.K
 def primitive(self,t):
  t=np.asarray(t,dtype=LD);z=(t%1)*self.K;j=np.floor(z).astype(np.int64);r=z-j
  u=(self.labels[j%self.K].astype(LD)+LD('.5'))/self.B
  pp=np.clip((r-u+self.h)/(2*self.h),0,1)-r
  low=np.zeros_like(t)
  for m,c in enumerate(self.co,1):low+=c/2*np.sin(2*PI*m*(t%1))
  return t+self.p/self.K*(pp+low)
 def gate_constants(self,intervals):
  return [(ld(a),ld(b),self.primitive(ld(a)).item(),self.primitive(ld(b)).item()) for a,b in intervals]
 def gate_prefix(self,t,Ft,constants):
  per=np.floor(t);r=t-per;Fr=Ft-per;pi=sum(Fb-Fa for a,b,Fa,Fb in constants);res=per*pi
  for a,b,Fa,Fb in constants:
   res+=np.where(r>=b,Fb-Fa,np.where(r>a,Fr-Fa,LD(0)))
  return res
 def integrate(self,L,order):
  L=ld(L);cuts=[LD(0),LD(1)]
  for sign in [-1,1]:
   cuts.extend(((self.atoms[:,None]+np.array([-1,1],dtype=LD)*self.width/2+sign*L/2)%1).ravel())
   cuts.extend((np.array([0,.25,.5,.75],dtype=LD)+sign*L/2)%1)
  cuts=np.unique(np.array(cuts,dtype=LD));gx,gw=np.polynomial.legendre.leggauss(order)
  x=((cuts[:-1,None]+cuts[1:,None])/2+(cuts[1:]-cuts[:-1])[:,None]*gx/2).ravel()
  w=((cuts[1:]-cuts[:-1])[:,None]*gw/2).ravel();hi=x+L/2;lo=x-L/2
  Fhi=self.primitive(hi);Flo=self.primitive(lo);Y=(Fhi-Flo)/L
  assert Y.min()>0
  Qw=w*Y;mass=Qw.sum();logY=np.log1p(Y-1);my=(Qw*logY).sum()/mass;vy=(Qw*(logY-my)**2).sum()/mass
  rows=[]
  for name,intervals in GATES.items():
   constants=self.gate_constants(intervals);pi=sum(Fb-Fa for a,b,Fa,Fb in constants)
   N=(self.gate_prefix(hi,Fhi,constants)-self.gate_prefix(lo,Flo,constants))/L
   post=N/Y;pm=(Qw*post).sum()/mass
   vp=(Qw*(post-pi)**2).sum()/mass;cov=(Qw*(post-pi)*(logY-my)).sum()/mass
   corr=cov/np.sqrt(vp*vy)
   rows.append(dict(gate=name,pi=float(pi),posterior_Qmean=float(pm),posterior_min=float(post.min()),posterior_max=float(post.max()),
      Var_Q_p=float(vp),Var_Q_logY=float(vy),r_squared=float(vp/(pi*(1-pi))),Corr_Q_p_logY=float(corr),
      Corr_Q_class_logY=float(cov/np.sqrt(pi*(1-pi)*vy)),Cov_Q_p_logY=float(cov),Q_mass=float(mass),quadrature_order=order,nodes=len(x)))
  return rows

def small_suite(which='pilot'):
 inputs=INPUT_DIR
 configs=[(10,F(1,4),32)] if which=='pilot' else [(10,F(1,4),32),(12,F(99,100),40)]
 result=[]
 for d,p,M in configs:
  lab=NPZLabels(inputs/f'source_d{d}.npz',3);o=SourceOracle(lab,F(1,64*lab.K*1024),p=p,M=M);ss=SmallSource(o)
  lengths=[F(7,6)] if which=='pilot' else [F(10,9),F(7,6),F(11,9)]
  for L in lengths:
   before=time.perf_counter();byorder=[ss.integrate(L,k) for k in [4,8]]
   comparisons=[]
   for a,b in zip(*byorder):
    diffs={k:abs(a[k]-b[k]) for k in ['pi','r_squared','Corr_Q_p_logY','Corr_Q_class_logY','Q_mass','posterior_Qmean']}
    comparisons.append(dict(**b,order4_vs8_absolute_differences=diffs))
   result.append(dict(d=d,K=lab.K,B=lab.B,p=str(p),M=M,v=1024,w=64,delta=str(o.delta),L=str(L),
       source_triangle_lower=float(1-p-F(20)*p/(3*lab.K)*sum((F(1,m) for m in range(1,M+1)),F(0))),
       rows=comparisons,seconds=time.perf_counter()-before))
 return result

def exact_geometry_reference():
 rows=[]
 for L in [F(10,9),F(7,6),F(11,9)]:
  def pre(t):
   n=t.numerator//t.denominator;r=t-n
   return F(n,2)+min(r,F(1,2))
  def pp(x):return (pre(x+L/2)-pre(x-L/2))/L
  cuts=sorted({F(0),F(1)}|{(a+sign*L/2)%1 for a in [F(0),F(1,2)] for sign in [-1,1]})
  variance=F(0);center=F(0)
  for a,b in zip(cuts,cuts[1:]):
   u,z=pp(a)-F(1,2),pp(b)-F(1,2)
   variance+=(b-a)*(u*u+u*z+z*z)/3;center+=(b-a)*(u+z)/2
  ell=L-1;closed=ell*ell/(L*L)*(1-F(4,3)*ell)
  assert center==0 and 4*variance==closed
  rows.append(dict(L=str(L),uniform_half_r_squared_exact=str(closed),numeric=float(closed)))
 K=2**47-1;p=F(99,100);V=F(40,9);Lmin=F(10,9)
 D=p/K*(1+V/2);kap=p/(K*Lmin)*(2+V);ep=(4*D/Lmin+kap)/(1-kap)
 er=(kap+2*ep+3*D+4*D*D)/(F(1,4)-4*D*D)
 assert F(76,3267)+er<F(1,36)
 return dict(rows=rows,high_d47_uniform_in_L_absolute_r_squared_error_bound=str(er),
             error_bound_numeric=float(er),high_d47_r_strictly_less_than='1/6',
             source_CDF_discrepancy_bound=str(D),Y_discrepancy_bound=str(kap),posterior_discrepancy_bound=str(ep))

def high_pilot(manifest,N=32):
 params=json.loads(Path(manifest).read_text())['parameters'];lab=GFLabels.from_manifest(manifest)
 p=F(params['p']);n=int(params['n']);v=9*p*p*n/(29*lab.K**2);delta=1/(int(params['w'])*lab.K*v)
 o=SourceOracle(lab,delta,p=p,M=int(params['M']));rng=np.random.default_rng(510001)
 xs=[o.point(((2*i+1)*lab.K)//(2*N),F(int(rng.integers(1,2**20)),2**20)) for i in range(N)]
 constants={name:[(a,b,primitive_exact(o,a),primitive_exact(o,b)) for a,b in intervals] for name,intervals in GATES.items()}
 def gatepart(t,Ft,cons):
  per=t.numerator//t.denominator;r=t-per;Fr=Ft-per;pi=sum((Fb-Fa for a,b,Fa,Fb in cons),mp.mpf(0));out=per*pi
  for a,b,Fa,Fb in cons:
   if r>=b:out+=Fb-Fa
   elif r>a:out+=Fr-Fa
  return out
 rows=[]
 for L in [F(10,9),F(7,6),F(11,9)]:
  points=[]
  for x in xs:
   hi=x+L/2;lo=x-L/2;Fhi=primitive_exact(o,hi);Flo=primitive_exact(o,lo)
   Y=(Fhi-Flo)/mpf(L);post={name:(gatepart(hi,Fhi,cons)-gatepart(lo,Flo,cons))/mpf(L)/Y for name,cons in constants.items()}
   points.append((Y,mp.log(Y),post))
  mass=sum((z[0] for z in points),mp.mpf(0))/N
  ml=sum((Y*g for Y,g,pp in points),mp.mpf(0))/(N*mass)
  vg=sum((Y*(g-ml)**2 for Y,g,pp in points),mp.mpf(0))/(N*mass)
  for name,cons in constants.items():
   pi=sum((Fb-Fa for a,b,Fa,Fb in cons),mp.mpf(0));pm=sum((Y*pp[name] for Y,g,pp in points),mp.mpf(0))/(N*mass)
   vp=sum((Y*(pp[name]-pm)**2 for Y,g,pp in points),mp.mpf(0))/(N*mass)
   cv=sum((Y*(pp[name]-pm)*(g-ml) for Y,g,pp in points),mp.mpf(0))/(N*mass)
   rows.append(dict(L=str(L),gate=name,pi=mp.nstr(pi,50),pilot_mean_post=mp.nstr(pm,50),pilot_r_squared=mp.nstr(vp/(pi*(1-pi)),30),
    pilot_Corr_p_logY=mp.nstr(cv/mp.sqrt(vp*vg),30),pilot_Corr_class_logY=mp.nstr(cv/mp.sqrt(pi*(1-pi)*vg),30),
    pilot_Var_logY=mp.nstr(vg,30),pilot_query_mass=mp.nstr(mass,40)))
 return dict(d=lab.d,b=lab.b,K=lab.K,p=str(p),M=o.M,v=str(v),delta=str(delta),
             query_design=f'{N} stratified integer coarse cells plus exact rational local phases, seed510001',
             all_geometry_Fraction=True,precision_decimal_digits=mp.mp.dps,rows=rows,
             scope='Few-query descriptive pilot, not a quadrature/MC error certificate. All tiny-width comparisons use Fraction before70-digit arithmetic.')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--high-pilot',action='store_true');ap.add_argument('--manifest',default=str(HERE.parent/'round046_gf/MANIFEST.json'));args=ap.parse_args();start=time.perf_counter()
 if args.high_pilot:
  result=high_pilot(args.manifest);out=dict(high_pilot=result,elapsed_seconds=time.perf_counter()-start,
      manifest_sha256=hashlib.sha256(Path(args.manifest).read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
  (HERE/'HIGH_PILOT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out));return
 rows=small_suite('pilot' if args.pilot else 'full')
 out=dict(scope='Actual positive periodic sources. Resolved longdouble piecewise Gauss quadrature; no interval certificate or high-dimensional tail transfer.',rows=rows,exact_geometry=exact_geometry_reference(),elapsed_seconds=time.perf_counter()-start,
          arithmetic=dict(longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant)+1,longdouble_epsilon=float(np.finfo(np.longdouble).eps),mpmath_decimal_digits=mp.mp.dps),
          inputs_dir=str(INPUT_DIR.resolve()),input_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in INPUT_DIR.glob('source_d*.npz')},
          oracle_source_sha256=hashlib.sha256((ORACLE_DIR/'source_oracle.py').read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 path=HERE/('PILOT.json' if args.pilot else 'RESULTS.json');path.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out))
if __name__=='__main__':main()
