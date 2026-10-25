"""One conservative positive-rate matrix + one Cauchy semigroup analytic family.
No continuum obstacle solver, kernel cost search, or new weak-type bound.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
import sympy as s
import mpmath as mp

N=5;t=Q(1,4);kap=Q(1)
rate=[[Q(0) for _ in range(N)] for _ in range(N)]
for i in range(N):
    for j in range(N):
        if i!=j:
            d=min(abs(i-j),N-abs(i-j));rate[i][j]=Q(1,d*d)
    rate[i][i]=-sum(rate[i])
P=[[Q(int(i==j))+t*rate[i][j] for j in range(N)] for i in range(N)]
assert all(sum(row)==1 for row in P)
assert all(P[i][j]>0 and P[i][j]==P[j][i] for i in range(N) for j in range(N))
mv=lambda v:[sum(P[i][j]*v[j] for j in range(N)) for i in range(N)]
def update(v,f):
    pv=mv(v)
    return [max(Q(0),pv[i]+t*(f[i]-kap)) for i in range(N)]
f=[Q(2)]+[Q(0)]*4;u=[Q(2,5)]+[Q(0)]*4
assert update(u,f)==u
pu=mv(u);sigma=[f[i]+(pu[i]-u[i])/t for i in range(N)]
assert sigma==[Q(1),Q(2,5),Q(1,10),Q(1,10),Q(2,5)]
assert all(0<=z<=kap for z in sigma) and sum(sigma)==sum(f)
energy=sum(u[i]*(u[i]-pu[i])/t for i in range(N))
assert energy+kap*sum(u)==sum(f[i]*u[i] for i in range(N))
v=[Q(0)]*N
for _ in range(40):v=update(v,f)
error=max(abs(v[i]-u[i]) for i in range(N))
assert error==Q(2,5)*Q(3,8)**40
superf=[Q(6)]+[Q(0)]*4;w=[Q(0)]*N
for _ in range(20):
    nxt=update(w,superf);assert sum(nxt)-sum(w)>=t*(sum(superf)-N*kap);w=nxt
assert sum(w)>=Q(5)

z,M,ss,T=s.symbols('z M ss T',positive=True)
U=M*s.exp(-ss*z)*(1-s.exp(-T*z))/z
assert s.simplify(-z*U-M*(s.exp(-(ss+T)*z)-s.exp(-ss*z)))==0
assert s.limit(U,z,0)==M*T
mp.mp.dps=40;analytic=[]
for n in (3,4):
    r=s.symbols('r',positive=True)
    peak=s.integrate(M/(s.pi**n*(ss+r)**n),(r,0,T))
    expected=M/(s.pi**n*(n-1))*(ss**(1-n)-(ss+T)**(1-n))
    assert s.simplify(peak-expected)==0
    assert Q(100)/(Q(22,7)**n)>1 and Q(100,12**n)<1
    potentials=[]
    for dt in (mp.mpf('.5'),mp.mpf('.25'),mp.mpf('.125')):
        vv=100/mp.pi**n*dt**(1-n)*mp.zeta(n,1/dt)
        bound=dt*100/mp.pi**n+100/(mp.pi**n*(n-1))
        assert vv<=bound
        potentials.append({'t':str(dt),'discrete_potential_peak':mp.nstr(vv,30),
                           'right_sum_bound':mp.nstr(bound,30)})
    analytic.append({'n':n,'input':'f=100*p_1, T=3, kappa=1',
                     'uT_L1_exact':'300','f_peak':mp.nstr(100/mp.pi**n,25),
                     'sigma_peak':mp.nstr(100/(4*mp.pi)**n,25),
                     'uT_peak':mp.nstr(100/(mp.pi**n*(n-1))*(1-4**(1-n)),25),
                     'full_potential_peak':mp.nstr(100/(mp.pi**n*(n-1)),25),
                     'full_potential_L1':'infinite by Tonelli',
                     'saturation':'fails: uT>0 everywhere, sigma not kappa',
                     'discrete_potential_peak_checks':potentials})

data={'batch':'U08','status':'finite diagnostics PASS, not a continuum numerical proof or new upper bound',
      'matrix':{'N':N,'t':str(t),'P_exact':[[str(x) for x in row] for row in P],
                'fixed_u_exact':list(map(str,u)),'sigma_exact':list(map(str,sigma)),
                'mass_exact':str(sum(sigma)),'energy_exact':str(energy),
                '40_step_error_exact':str(error),'40_step_error_display':float(error),
                'supercritical_20step_mass_exact':str(sum(w)),
                'supercritical_20step_mass_display':float(sum(w)),
                'supercritical_obstruction':'input mass6>total cap5, so no finite mass-conserving capped fixed point'},
      'analytic_family':analytic,'new_upper_bound':False,
      'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('u08_results.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
