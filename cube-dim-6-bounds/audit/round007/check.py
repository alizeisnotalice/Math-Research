#!/usr/bin/env python3
"""F007: three bounded analytic checks, no body32 kernel evaluation."""
import hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
import sympy as s
from flint import arb,ctx
ctx.prec=256
HERE=Path(__file__).resolve().parent
q=F(6,5)
def mul(a,b,N):
    c=[F(0) for _ in range(N+1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=N:c[i+j]+=x*y
    return c
def constant(N):
    a=[F(0)]+[F(math.comb(63+k,k),32**k) for k in range(1,N+1)]
    power=[F(1)]+[F(0)]*N
    binom=F(1);total=F(0)
    for j in range(1,N+1):
        power=mul(power,a,N)
        binom*= (q-j+1)/j
        total+=abs(binom)*power[N]
    return math.factorial(N)*total
D6=constant(6);D12=constant(12)
assert D12==F(2354965600826387008780486926466707,8589934592000000000000)
c6=F(1,46080)
def ar(x):return arb(x.numerator)/x.denominator
eps=arb(1)/1000;B=arb(1)/16000
A=arb(243)**(-arb(1)/6)*arb.pi()**(-arb(11)/6)
a=arb(12)-arb(7)/6
M12upper=A*eps**(-a)/a
M12lower=M12upper*(1-arb(2)**(-a))/(1+4096*eps**2)
pair_at_h_eps_lower=ar(c6*c6*D12)*eps**12*M12lower/B
pair_at_h_eps_upper=ar(c6*c6*D12)*eps**12*M12upper/B
all_pairs_lower=10*pair_at_h_eps_lower
h_over_eps_necessary=all_pairs_lower**(-arb(1)/12)
assert all_pairs_lower>29
assert h_over_eps_necessary<arb(76)/100

# Test2: positive smooth F with zero integrated one-axis errors but nonzero tensor error.
x,y=s.symbols('x y',real=True)
nodes=[(5-s.sqrt(15))/10,s.Rational(1,2),(5+s.sqrt(15))/10]
weights=[s.Rational(5,18),s.Rational(4,9),s.Rational(5,18)]
f=x**6-s.Rational(1,7)
ef=s.simplify(sum(w*f.subs(x,v) for w,v in zip(weights,nodes)))
assert ef==-s.Rational(1,2800)
Fxy=1+(x**6-s.Rational(1,7))*(y**6-s.Rational(1,7))
II=s.integrate(Fxy,(x,0,1),(y,0,1))
QQ=s.simplify(sum(wx*wy*Fxy.subs({x:nx,y:ny}) for nx,wx in zip(nodes,weights) for ny,wy in zip(nodes,weights)))
assert II==1 and QQ-II==s.Rational(1,7840000)
assert s.integrate(f,(x,0,1))==0

# Test3: exact single-factor extremal ratio and new sparse-gauge scalar factors.
u=s.Rational(10,3);z=s.symbols('z',real=True)
M=(s.sqrt(u*u+4)+u)**2/4
quadratic=s.expand(M*(1+(z+u)**2)-(1+z*z))
assert s.simplify(s.discriminant(quadratic,z))==0
shift_rows=[]
for r in (5,20,100):
    phi=(12*(arb(r)/12).asinh()).exp()
    old=arb(r).exp()
    assert phi<old
    shift_rows.append(dict(r_over_t=r,new_factor=phi.str(50),old_over_new=(old/phi).str(50)))

inputs=[HERE.parents[1]/'upper/round007/mixed_time_error/report.md',HERE/'sparse_shift_bound.md']
# Test3 continuation: independent scalar replay of A's saved best local box.
# No fan evaluation or new kernel point is performed.
apath=HERE.parents[1]/'upper/round007/g_time_moments/results.json'
adata=json.loads(apath.read_text())
hs=[arb(row['H_arb384']) for row in adata['time_nodes']]
taus=[F(row['tau']) for row in adata['time_nodes']]
grad=[arb(v) for v in next(row for row in adata['gradients'] if row['time_index']==4)['gradient_x2_to_x6']]
rad=F(717,50000)
centers=[F(v,adata['denominator']) for v in adata['center_numerators'][1:]]
intervals=[(max(-rad,-v),min(rad,1-v)) for v in centers]
linear=arb(0);normsq=F(0);volume=F(120)
for g,(lo,hi) in zip(grad,intervals):
    endpoints=[(g*ar(lo)).lower(),(g*ar(hi)).lower()]
    linear+=endpoints[0] if endpoints[0]<endpoints[1] else endpoints[1]
    normsq+=max(abs(lo),abs(hi))**2
    volume*=hi-lo
U2=arb(0)
for k in range(4,7):
    mass_upper=(hs[k]-hs[k+1]).upper()
    U2+=ar(taus[k])**-2*(12*(ar(rad/taus[k])/12).asinh()).exp()*mass_upper
U2+=ar(taus[7])**-2*(12*(ar(rad/taus[7])/12).asinh()).exp()*hs[7].upper()
local_lower=hs[4].lower()+linear-ar(normsq)/32*U2.upper()
assert local_lower>arb(16666817)/10**12>arb(1)/60000
coordinate_intervals=[(v+lo,v+hi) for v,(lo,hi) in zip(centers,intervals)]
assert all(coordinate_intervals[i+1][1]<coordinate_intervals[i][0] for i in range(4))
assert volume==F(614593713570664440310497,488281253417968750000000000000)
inputs.append(apath)
out=dict(status='PASS_THREE_ANALYTIC_TESTS_NO_GLOBAL_C6_BOUND',kernel_evaluations=0,
 test1=dict(D6=str(D6),D12=str(D12),D12_over_D6_squared=str(D12/D6**2),
   pair_budget_at_h_eps_lower=pair_at_h_eps_lower.str(50),pair_budget_at_h_eps_upper=pair_at_h_eps_upper.str(50),
   ten_pair_budget_at_h_eps_lower=all_pairs_lower.str(50),
   necessary_h_over_eps_upper=h_over_eps_necessary.str(50),
   limitation='Allocated scalar envelope budget only; no lower bound on actual derivative error or model cost.'),
 test2=dict(F='1+(x^6-1/7)*(y^6-1/7)',positive_lower='43/49',II=str(II),QQ=str(QQ),
  E1I2='0',I1E2='0',E1E2='1/7840000',omitting_mixed_fails=True),
 test3=dict(single_factor_discriminant='0',sparse_gauge_support_bound=6,
  majorant='exp(12*asinh(r/(12*t)))',values=shift_rows,
  A_saved_box_replay=dict(epsilon='2/125',radius=str(rad),
    independent_G_lower=local_lower.str(70),strict_rational_G_lower='16666817/1000000000000',
    independent_J2_upper=U2.str(70),permutation_images_disjoint_interiors=True,
    exact_covered_positive_face_volume=str(volume),
    global_coverage=False,kernel_replayed=False)),
 inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
