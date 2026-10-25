"""3D body+face Cauchy frame, exact 96-cone log representation.

L_tau=-sum_body |v.D|-tau sum_face |w.D|, tau>=0, unnormalized rays.
Finite Green values are analytic log expressions. Zero-knot divergences are
combined exactly before taking epsilon->0; no cone/negative term is removed.
Numerics are not interval certificates. Only n=3.
"""
from fractions import Fraction as F
from itertools import product,permutations
from functools import lru_cache
import math
import numpy as np

def _q(v):
    if isinstance(v,F):return v
    if isinstance(v,(int,np.integer)):return F(int(v))
    return F(str(v))

def rate_coeff(r):
    a,b,c=sorted(map(abs,r),reverse=True)
    return max(4*a,2*(a+b+c)),2*(2*a+b)

def _cones():
    bases=[((1,0,0),(1,1,0),(2,1,1)),((1,1,0),(2,1,1),(1,1,1))]
    out=[]
    for p in permutations(range(3)):
        for sg in product((-1,1),repeat=3):
            for base in bases:out.append(tuple(tuple(sg[i]*r[p[i]] for i in range(3)) for r in base))
    assert len(out)==96 and len({tuple(sorted(c)) for c in out})==96
    return tuple(out)
CONES=_cones()

@lru_cache(maxsize=32)
def _prepared(tau):
    tau=_q(tau)
    if tau<0:raise ValueError('tau must be nonnegative')
    seen=set();out=[]
    for cone in CONES:
        key=tuple(sorted(cone));opp=tuple(sorted(tuple(-v for v in r) for r in cone));pair=min(key,opp)
        if pair in seen:continue
        seen.add(pair);rates=tuple(F(a)+tau*b for a,b in map(rate_coeff,cone));weight=F(2)/math.prod(rates)
        out.append((cone,rates,weight)) # opposite cone is included by factor 2
    assert len(out)==48
    return tuple(out)

def _add(d,k,v):
    k=abs(k);d[k]=d.get(k,F(0))+v

def green_expression(x,tau):
    """Exact rational coefficients of constant+sum c_b log|b|+D log epsilon.

    Coordinates/tau may be Fraction or decimal strings. Float inputs are
    interpreted as their decimal spelling. Returned D<0 means true +infinity;
    D=0 gives a finite expression. Origin is returned separately.
    """
    x=tuple(map(_q,x));tau=_q(tau)
    if len(x)!=3:raise ValueError('x must be a 3-vector')
    scale=max(map(abs,x))
    if not scale:return {'origin':True,'scale':F(0),'log_divergence':None}
    x=tuple(v/scale for v in x);constant=F(0);D=F(0);logs={};cache={}
    for cone,rates,W in _prepared(tau):
        nodes=[]
        for r,k in zip(cone,rates):
            key=(r,k)
            if key not in cache:cache[key]=sum(F(v)*xx for v,xx in zip(r,x))/k
            nodes.append(cache[key])
        zeros=nodes.count(0)
        if zeros==2:
            b=next(v for v in nodes if v);D-=W/(b*b);_add(logs,b,W/(b*b))
        elif zeros==1:
            b,c=[v for v in nodes if v];D+=W/(b*c)
            if b==c:_add(logs,b,-W/(b*b));constant+=W/(b*b)
            else:_add(logs,b,W/(b*(b-c)));_add(logs,c,W/(c*(c-b)))
        elif zeros==3:raise ArithmeticError('nonzero x orthogonal to an independent cone')
        elif len(set(nodes))==1:constant-=W/(2*nodes[0]**2)
        elif len(set(nodes))==2:
            a=next(v for v in nodes if nodes.count(v)==2);b=next(v for v in nodes if v!=a)
            _add(logs,b,W/(b-a)**2);_add(logs,a,-W/(b-a)**2);constant-=W/(a*(b-a))
        else:
            for j,b in enumerate(nodes):_add(logs,b,W/math.prod(b-nodes[k] for k in range(3) if k!=j))
    logs={a:c for a,c in logs.items() if c};assert sum(logs.values(),F(0))+D==0
    if D>0:raise ArithmeticError('positive log-epsilon coefficient contradicts true density positivity')
    return {'origin':False,'scale':scale,'constant':constant,'log_coefficients':logs,'log_divergence':D}

def _mpq(q,mp):return mp.mpf(q.numerator)/q.denominator

def green(x,tau, *, method='auto',dps=80):
    """True G_tau(x). auto returns float; mp returns an mpf at fixed dps.

    Exact nonzero D log(epsilon) divergence is preserved as +infinity.
    Finite expressions retain every signed cone contribution.
    """
    e=green_expression(x,tau)
    if e['origin'] or e['log_divergence']<0:return math.inf
    logs=e['log_coefficients'];reference=max(logs) if logs else F(1)
    if method not in ('auto','float','mp'):raise ValueError('method must be auto, float or mp')
    terms=[float(e['constant'])]
    for a,c in logs.items():
        ratio=a/reference
        lg=math.log1p(float((a-reference)/reference)) if ratio>F(1,2) else math.log(float(ratio))
        terms.append(float(c)*lg)
    total=math.fsum(terms);condition=sum(abs(v) for v in terms)/max(abs(total),1e-300)
    if method=='float' or (method=='auto' and total>0 and condition<1e5):return total/(8*math.pi**3*float(e['scale'])**2)
    import mpmath as mp
    with mp.workdps(dps):
        val=_mpq(e['constant'],mp)+mp.fsum(_mpq(c,mp)*mp.log(_mpq(a/reference,mp)) for a,c in logs.items())
        val/=8*mp.pi**3*_mpq(e['scale'],mp)**2
        if val<=0:raise ArithmeticError('finite exact log expression evaluated nonpositive; raise precision explicitly')
        return val if method=='mp' else float(val)

def green_face(y,z,tau,**kwargs):return green((1,y,z),tau,**kwargs)

def green_array(x,tau,**kwargs):
    x=np.asarray(x,dtype=float)
    if x.shape[-1]!=3:raise ValueError('last dimension must be 3')
    return np.array([green(v,tau,**kwargs) for v in x.reshape(-1,3)]).reshape(x.shape[:-1])

@lru_cache(maxsize=1024)
def _heat_data(x,tau):
    x=tuple(map(_q,x));rows=[];weights=[]
    for cone,rates,w in _prepared(tau):
        rows.append([float(sum(F(v)*xx for v,xx in zip(r,x))/k) for r,k in zip(cone,rates)]);weights.append(float(w))
    return np.array(rows),np.array(weights)

def heat_kernel(t,x,tau):
    if t<=0:raise ValueError('t must be positive')
    b,w=_heat_data(tuple(map(_q,x)),_q(tau))
    return float(np.real(np.sum(w/np.prod(float(t)-1j*b,axis=1)))/(8*math.pi**3))

def peak_coefficient(tau):
    """Exact rational C with p_t(0)=C/(pi^3 t^3)."""
    tau=_q(tau)
    if tau<0:raise ValueError('tau must be nonnegative')
    return F(5)/(8*(1+tau)*(2+3*tau)*(4+5*tau))

def green_time_tail_upper(T,tau):return float(peak_coefficient(tau))/(2*math.pi**3*float(T)**2)
