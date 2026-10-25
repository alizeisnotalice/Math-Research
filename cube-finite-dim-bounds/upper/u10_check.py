"""Fixed U10 support/mass/tail arithmetic, not a candidate/LP or source integral."""
from pathlib import Path
import hashlib,json,itertools
import sympy as s
import mpmath as mp
q=s.Rational
centers=[s.Matrix([sg*q(3,4) if j==i else 0 for j in range(3)]) for i in range(3) for sg in (-1,1)]
d2=sorted(set(s.simplify((a-b).dot(a-b)) for a,b in itertools.combinations(centers,2)))
assert d2==[q(9,8),q(9,4)] and all(x>1 for x in d2)
a=q(5,4);R=q(1,2)
assert q(3,4)+R==a and q(3,2)-a==q(1,4)
r=s.symbols('r',positive=True)
M0=s.integrate(4*s.pi*(1-r*r)**2*r*r,(r,0,1));M1=6*R**3*M0
assert M0==32*s.pi/105 and M1==8*s.pi/35 and M1/M0==q(3,4)
assert s.simplify(5/(16*s.pi**3*a*a)-1/(8*s.pi**3))==3/(40*s.pi**3)
assert s.simplify(5/(16*s.pi**3*q(3,2)**2)-1/(8*s.pi**3))==1/(72*s.pi**3)
V=[s.Matrix(v) for v in [(1,1,1),(1,1,-1),(1,-1,1),(1,-1,-1)]]
dets=[abs(s.Matrix.hstack(*[V[i] for i in ids]).det()) for ids in itertools.combinations(range(4),3)]
assert dets==[4]*4
B=s.Matrix.hstack(*V[:3]);Bi=B.inv();extra=Bi*V[3]
assert all(abs(z)==1 for z in extra)
assert all(sum(Bi[i,j]**2 for j in range(3))<=1 for i in range(3))
tau,u=s.symbols('tau u',nonnegative=True)
y=V[0]+tau*V[1]
dist2=[s.expand(y.dot(y)-y.dot(v)**2/3) for v in V]
assert all(all(z>=0 for z in s.Poly(s.expand((z-q(8,3)).subs(tau,u+1)),u).all_coeffs()) for z in dist2)
assert q(8,3)-q(3,2)**2==q(5,12)>0
# min(A t,C t^-3) time integral and fixed frame constants.
A,C,d=s.symbols('A C d',positive=True)
t=s.symbols('t',positive=True);t0=(C/A)**q(1,4)
upper=s.integrate(A*t,(t,0,t0))+s.integrate(C*t**-3,(t,t0,s.oo))
assert s.simplify(upper-s.sqrt(A*C))==0
Cframe=1/(4*s.pi**3);delta=d/(8*s.sqrt(3));Aframe=6*Cframe/delta**4
Gupper=s.simplify(s.sqrt(Aframe*Cframe))
assert s.simplify(Gupper-48*s.sqrt(6)/(s.pi**3*d*d))==0
assert 2304*q(5,2)/q(5,3)==3456 and 3456<60**2
W,cc=s.symbols('W c',positive=True)
LB=s.simplify(4*(2*a/s.sqrt(3))*4*W/(s.pi*t*t))
assert s.simplify(LB-40*W/(s.pi*s.sqrt(3)*t*t))==0
pay=s.factor(cc/(4*s.pi*t)-LB)
assert s.simplify(pay-(s.sqrt(3)*cc*t-160*W)/(4*s.pi*s.sqrt(3)*t*t))==0
mp.mp.dps=50
c=(4*mp.log(2)-1)/(24*mp.pi**3)
syntheticW=mp.mpf('0.002');Ts=max(mp.mpf(60),160*syntheticW/(mp.sqrt(3)*c));Rs=mp.sqrt(3)*Ts+mp.mpf('1.25')
data={'batch':'U10','status':'fixed exact arithmetic PASS, no candidate or LP/source integration',
      'support':{'centers_distance_squared':list(map(str,d2)),'outer_radius':str(a),'Omega_buffer_lower':'1/4','b1_pointwise_upper':'1'},
      'mass':{'M0':str(M0),'M1':str(M1),'M1_over_M0':'3/4','Lb_small_ball_scaling_factor':'2'},
      'line_basis':{'triple_determinants':list(map(str,dets)),'remaining_basis_coordinates':list(map(str,extra))},
      'tail':{'segment_distance_squared':list(map(str,dist2)),'G_off_line_upper':str(Gupper),
              'uniform_time_threshold_sufficient':'60','joint_source_slack_factor':str(pay),
              'T_formula':'max(60,160*(w0+w1)/(sqrt(3)*c))','R_formula':'sqrt(3)*T+5/4'},
      'synthetic_arithmetic_only':{'w0':'1/1000','w1':'1/1000','c_display':mp.nstr(c,25),'T_display':mp.nstr(Ts,25),'R_display':mp.nstr(Rs,25),'candidate_feasibility':'NOT tested or claimed'},
      'new_upper_bound':False,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('u10_results.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(data,indent=2,ensure_ascii=False))
