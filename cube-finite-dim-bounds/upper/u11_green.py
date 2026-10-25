"""4D body Cauchy frame: exact cone geometry and Green interfaces.
Generic exact G is a rational coefficient times pi^-3.
Zero nodes: use the ALL-CONE positive-time cutoff, never guessed sign(0).
Numerical outputs do not certify floating roundoff.
"""
from fractions import Fraction as Q
from itertools import permutations,product
from math import prod,pi
import numpy as np
import mpmath as mp

V=tuple((1,)+s for s in product((-1,1),repeat=3))
BASE_RAYS={'A':(1,0,0,0),'B':(1,1,0,0),'C':(1,1,1,0),
           'D':(1,1,1,1),'E':(2,1,1,0),'F':(3,1,1,1)}
BASE_CONES=((('A','B','E','F'),1),(('E','B','C','D'),1),(('E','B','D','F'),2))
RAYS=[];CONES=[];_ids={}
for perm in permutations(range(4)):
 for signs in product((-1,1),repeat=4):
  for labels,det in BASE_CONES:
   ids=[]
   for lab in labels:
    r=tuple(signs[j]*BASE_RAYS[lab][perm[j]] for j in range(4))
    if r not in _ids:_ids[r]=len(RAYS);RAYS.append(r)
    ids.append(_ids[r])
   CONES.append((tuple(ids),det))
RAYS=tuple(RAYS);CONES=tuple(CONES)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def rates(tau=Q(0)):
 tau=Q(tau)
 if tau<0:raise ValueError('nonnegative face weight required')
 return tuple(Q(sum(abs(dot(v,r)) for v in V))+tau*sum(abs(r[i]+r[j])+abs(r[i]-r[j]) for i in range(4) for j in range(i+1,4)) for r in RAYS)
def point_nodes(x,tau=Q(0)):
 x=tuple(map(Q,x));q=rates(tau)
 return x,q,tuple(dot(x,r)/qj for r,qj in zip(RAYS,q))
def dd_sign(nodes):
 x=sorted(map(Q,nodes))
 if any(z==0 for z in x):raise ValueError('zero node: use full cutoff interface')
 vals=[Q(1 if z>0 else -1) for z in x]
 for k in range(1,len(x)):
  vals=[Q(0) if x[i+k]==x[i] else (vals[i+1]-vals[i])/(x[i+k]-x[i]) for i in range(len(x)-k)]
 return vals[0]
def green_pi3_coefficient(x,tau=Q(0)):
 """Exact Fraction C with G(x)=C/pi^3, nonzero nodes; confluent repeats legal."""
 x,q,b=point_nodes(x,tau);ss=Q(0)
 for ids,det in CONES:
  ss+=Q(det)/prod(q[i] for i in ids)*dd_sign([b[i] for i in ids])
 return -ss/32

def generic_residues(x,tau=Q(0)):
 """Exact real-p(t) residues: p_t=(2pi)^-4 sum_{b>0} Cb/(t^2+b^2)."""
 x,q,b=point_nodes(x,tau);rr={}
 for ids,det in CONES:
  z=[b[i] for i in ids]
  if len(set(z))<4 or 0 in z:raise ValueError('strict generic residues required')
  alpha=Q(det)/prod(q[i] for i in ids)
  for j,bj in enumerate(z):rr[bj]=rr.get(bj,Q(0))+alpha/prod(bj-z[k] for k in range(4) if k!=j)
 assert all(rr.get(-bj,Q(0))==-val for bj,val in rr.items())
 return {bj:-2*rr[bj]*bj for bj in rr if bj>0 and rr[bj]}

def _mpq(z):z=Q(z);return mp.mpf(z.numerator)/z.denominator

def density_mp_from_residues(t,res):
 t=mp.mpf(t)
 return mp.fsum(_mpq(co)/(t*t+_mpq(b)**2) for b,co in res.items())/(16*mp.pi**4)

def green_time_quad(x,dps=50):
 with mp.workdps(dps):
  res=generic_residues(x);bs=[_mpq(b) for b in res]
  cuts=sorted(set([mp.mpf(0),min(bs),mp.mpf(1),max(bs)]))+[mp.inf]
  val=mp.quad(lambda t:density_mp_from_residues(t,res),cuts)
  return +val

def dd_cutoff_log(nodes,epsilon):
 """Complex confluent DD3 of log(epsilon-i*b), epsilon>0."""
 x=sorted(map(Q,nodes));eps=_mpq(epsilon) if isinstance(epsilon,(Q,int)) else mp.mpf(epsilon)
 if eps<=0:raise ValueError('positive cutoff required')
 vals=[mp.log(eps-1j*_mpq(b)) for b in x]
 for k in range(1,len(x)):
  vals=[-(1j)**k/(k*(eps-1j*_mpq(x[i]))**k) if x[i+k]==x[i]
        else (vals[i+1]-vals[i])/_mpq(x[i+k]-x[i]) for i in range(len(x)-k)]
 return vals[0]
def green_cutoff_mp(x,epsilon,dps=80,tau=Q(0)):
 """Integral_epsilon^infty p_t dt; valid for ALL nodes. Float roundoff not enclosed."""
 with mp.workdps(dps):
  x,q,b=point_nodes(x,tau)
  terms=[]
  for ids,det in CONES:
   terms.append(_mpq(Q(det)/prod(q[i] for i in ids))*mp.re(-1j*dd_cutoff_log([b[i] for i in ids],epsilon)))
  return +(mp.fsum(terms)/(16*mp.pi**4))
def body_line_distance_squared(x):
 x=tuple(map(Q,x));return min(dot(x,x)-dot(x,v)**2/4 for v in V)
def body_cutoff_time_error_bound(x,epsilon):
 """Analytic omitted-time bound A_d*epsilon, n4 body; not a roundoff bound."""
 d2=body_line_distance_squared(x)
 if d2==0:return mp.inf
 return _mpq(Q(3670016)*Q(epsilon)/(d2*d2))/mp.pi**4

_ri=np.array(RAYS,dtype=np.longdouble);_ci=np.array([i for i,d in CONES]);_det=np.array([d for i,d in CONES],dtype=np.longdouble)
def density_float(t,x,tau=0):
 if t<=0:raise ValueError('positive time required')
 q=np.array([float(z) for z in rates(Q(str(tau)))],dtype=np.longdouble)
 ell=_ri@np.array(x,dtype=np.longdouble)
 den=(np.longdouble(t)*q-1j*ell).astype(np.clongdouble)
 val=np.sum(_det/np.prod(den[_ci],axis=1),dtype=np.clongdouble)
 return float(val.real/(16*np.longdouble(pi)**4))
def green_float_generic(x,tau=0,relative_guard=1e-12):
 q=np.array([float(z) for z in rates(Q(str(tau)))])
 b=np.asarray(_ri,dtype=float)@np.asarray(x,dtype=float)/q
 z=b[_ci];scale=max(1.,float(np.max(np.abs(z))))
 if np.min(np.abs(z))<=relative_guard*scale:raise ValueError('near zero nodes: cutoff or exact method needed')
 dif=z[:,:,None]-z[:,None,:];dif[:,range(4),range(4)]=1
 if np.min(np.abs(dif))<=relative_guard*scale:raise ValueError('near repeated nodes: confluent method needed')
 dd=np.sum(np.sign(z)/np.prod(dif,axis=2),axis=1)
 same=np.all(z>0,axis=1)|np.all(z<0,axis=1);dd[same]=0
 return float(-np.sum(np.asarray(_det,dtype=float)/np.prod(q[_ci],axis=1)*dd)/(32*pi**3))

# Exact ALL-NODE residue interface: zeros/repeats are combined before integration.
from collections import Counter
from math import comb

def confluent_residues(x,tau=Q(0)):
 """Return rational R_(b,k), with total integrand sum i^(k-4)R/(t-i*b)^k."""
 x,q,b=point_nodes(x,tau);out={}
 for ids,det in CONES:
  gr=Counter(b[i] for i in ids);alpha=Q(det)/prod(q[i] for i in ids)
  for bj,m in gr.items():
   poly=[Q(1)]+[Q(0)]*(m-1)
   for bc,mc in gr.items():
    if bc==bj:continue
    d=bj-bc
    part=[Q((-1)**ell*comb(mc+ell-1,ell))/d**(mc+ell) for ell in range(m)]
    poly=[sum(poly[k]*part[j-k] for k in range(j+1)) for j in range(m)]
   for k in range(1,m+1):
    key=(bj,k);out[key]=out.get(key,Q(0))+alpha*poly[m-k]
 return {z:v for z,v in out.items() if v}

def frame_line_distance_squared(x,tau=Q(0)):
 x=tuple(map(Q,x));vv=list(V)
 if Q(tau)>0:
  for i in range(4):
   for j in range(i+1,4):
    for sg in (-1,1):
     v=[0]*4;v[i]=1;v[j]=sg;vv.append(tuple(v))
 return min(dot(x,x)-dot(x,v)**2/dot(v,v) for v in vv)

def green_allnodes_pi3_coefficient(x,tau=Q(0)):
 """Exact Fraction C with G=C/pi^3 at EVERY rational point off active frame lines.
 Detect infinities first. Zero poles must cancel in the FULL sum, never per cone.
 """
 if frame_line_distance_squared(x,tau)==0:raise ValueError('Green infinite on active frame line (including origin)')
 rr=confluent_residues(x,tau)
 if any(b==0 for (b,k),val in rr.items()):raise ArithmeticError('uncancelled zero pole off active frame lines')
 return -sum(val*(1 if b>0 else -1) for (b,k),val in rr.items() if k==1)/32

def density_mp_from_confluent(t,rr):
 t=mp.mpf(t)
 val=mp.fsum(mp.re((1j)**(k-4)*_mpq(co)/(t-1j*_mpq(b))**k) for (b,k),co in rr.items())
 return val/(16*mp.pi**4)
