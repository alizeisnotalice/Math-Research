#!/usr/bin/env python3
"""Fixed exact/512-bit truncated-Green identity checks. No expectation sampling."""
from pathlib import Path
from itertools import product
import hashlib,json,platform,time
import sympy as sp
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
SOURCES=['AGENTS.md','cube-dim-7-bounds/upper/round001/gaussian/proof.md',
 'cube-dim-7-bounds/upper/round001/theory/proof.md',
 'cube-dim-7-bounds/upper/round002/hadamard_moment/completeproof.md',
 'cube-dim-7-bounds/upper/round002/g_probe/report.md']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def snapshot():return {f:sha(ROOT/f) for f in SOURCES}
def ar(n,d=1):return arb(fmpq(n,d))
def enc(x):return {'arb':str(x),'lower_rational':str(x.lower().fmpq()),'upper_rational':str(x.upper().fmpq())}
start=time.perf_counter(); before=snapshot()
V=sp.Matrix.hstack(*(sp.Matrix((1,)+s) for s in product((-1,1),repeat=6)))
assert V*V.T==64*sp.eye(7)
assert 64*(V*V.T)==4096*sp.eye(7)
z=sp.symbols('z',positive=True); b=sp.Rational(65,2); c=b+3
F=1-(1-z)**b*(1+b*z+b*(b+1)*z*z/2)
B=2/(b*(b+1)*(b+2))
assert sp.simplify(sp.diff(F,z)-z*z*(1-z)**(b-1)/B)==0
assert F.subs(z,0)==0 and F.subs(z,1)==1
assert sp.limit(F/z**3,z,0)==b*(b+1)*(b+2)/6
assert sp.Rational(7,2)-3==sp.Rational(1,2)
assert sp.Integer(32)+sp.Rational(7,2)==sp.Rational(71,2)
assert sp.Integer(32)+sp.Rational(1,2)==b
# Exact prefactor and Q=0 consistency use Gamma recurrence.
assert sp.simplify(sp.gamma(c)*B/2-sp.gamma(b))==0
assert sp.simplify(sp.gamma(b)*b*(b+1)*(b+2)/6-sp.gamma(c)/6)==0
# Uniform determinant envelope at P uniform; degree and frame-leverage exponents.
assert (sp.Integer(4096)**7)==sp.Integer(64)**14
assert sp.Rational(7,64)/14==sp.Rational(1,128)
assert sp.Rational(39,64)*14==sp.Rational(273,32)
ctx.prec=512; pi=arb.pi(); bb=ar(65,2); cc=ar(71,2); eps=ar(1,10)
K0=cc.gamma()/(pi**ar(7,2)*ar(32).gamma())
KD=bb.gamma()/(pi**ar(7,2)*ar(32).gamma())
assert (K0*2/(bb*(bb+1)*(bb+2))/2-KD).contains(0)
checks=[]
for sn,sd in [(0,1),(1,64),(1,1),(64,1)]:
    s=ar(sn,sd); zz=s/(1+s)
    direct=zz.beta_lower(3,bb,regularized=True)
    closed=1-(1-zz)**bb*(1+bb*zz+bb*(bb+1)*zz*zz/2)
    assert (direct-closed).contains(0)
    if sn:
        profile=direct/s**3
        assert profile>0
    else:
        profile=bb*(bb+1)*(bb+2)/6
    checks.append({'s':f'{sn}/{sd}','I_beta':enc(direct),'closed_form':enc(closed),'dimensionless_Qminus3_profile':enc(profile),'difference_contains_zero':True})
kap=(ar(72).sqrt()**-1)*(ar(72,73)**ar(73,2))
theta=(ar(73,75)**ar(75,2)) # (c+1)/(c+2) = (73/2)/(75/2)
B0=K0/(6*ar(64)**7*eps**6)
Lx=K0*cc*kap/(28*ar(64)**7*eps**7)
Hx=K0*(4*cc*theta+2*cc)/(512*ar(64)**7*eps**8)
fourier={str(k):enc(ar(int(sp.gamma(6+k)))/(120*pi**4*ar(8)**(7+k)*eps**(6+k))) for k in range(3)}
assert snapshot()==before
out={'status':'PASS','claim':'exact truncated-Green Beta formula and deterministic error interfaces; no g or C7 certificate',
 'environment':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'arb_bits':ctx.prec},
 'formulas':{'b':'65/2','c':'71/2','K0':'Gamma(71/2)/(pi^(7/2)*Gamma(32))','KD':'Gamma(65/2)/(pi^(7/2)*Gamma(32))',
 'H_epsilon':'KD*E_P det(A)^(-1/2)*Q^(-3)*I_{Q/(Q+epsilon^2)}(3,65/2)',
 'Q_zero_profile':'b*(b+1)*(b+2)/(6*epsilon^6)',
 'P_bulk_first_chart_derivative':'(273/32)*B_epsilon/delta^2',
 'Fourier_order_k':'Gamma(6+k)/(120*pi^4*8^(7+k)*epsilon^(6+k))'},
 'scalar_beta_checks':checks,'fixed_epsilon':'1/10',
 'uniform_bound_expressions':{'K0':enc(K0),'KD':enc(KD),'shape_amplitude_B':enc(B0),'shape_gradient_L':enc(Lx),'shape_Hessian_operator_H':enc(Hx),'Fourier_H_derivatives_orders012':fourier},
 'sources':before,'registration_sha256':sha(HERE/'registration.md'),'check_sha256':sha(Path(__file__)),'elapsed_seconds':time.perf_counter()-start,
 'new_MC_samples':0,'g_certificate':False,'C7_certificate':False}
(HERE/'sources.json').write_text(json.dumps({'root':str(ROOT),'sha256':before},indent=2)+'\n')
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'epsilon':'1/10','b':'65/2','beta_cases':4,'B':str(B0),'L':str(Lx),'Hessian':str(Hx),'elapsed_seconds':out['elapsed_seconds'],'new_MC_samples':0},indent=2))
