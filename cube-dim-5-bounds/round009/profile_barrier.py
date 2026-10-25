"""New exact all-rho certificate-family obstruction; not an upper bound on G."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import sympy as s
import json,hashlib,time,os
P=Path(__file__).resolve().parents[1];began=time.monotonic();t=s.symbols('t',positive=True);z=s.symbols('z',positive=True);C=s.symbols('C',positive=True);rho=s.symbols('rho',real=True)
x=s.Matrix([1,1,1,0,0]);V=[s.Matrix((1,)+sig) for sig in product((-1,1),repeat=4)]
D=s.zeros(5);counts={9:0,1:0}
for v in V:
 h=int((v.dot(x))**2);assert h in counts;counts[h]+=1;D+=(t if h==9 else 1)*(v*v.T)
assert counts=={9:4,1:12}
detD=s.factor(D.det());bD=s.factor((x.T*D.adjugate()*x)[0]);assert s.expand(detD-16384*(3*t+1)*(t+3)**2)==0;assert s.expand(bD-12288*(t+3)**2)==0
assert s.simplify(D*x-(12*t+4)*x)==s.zeros(5,1)
I1=s.integrate(s.exp(-C*s.sqrt(z)),(z,0,s.oo));I32=s.integrate(z**s.Rational(3,2)*s.exp(-C*s.sqrt(z)),(z,0,s.oo));assert s.simplify(I1-2/C**2)==0 and s.simplify(I32-48/C**5)==0
Zinv=8*s.sqrt(s.pi);U2=s.simplify(Zinv**2*(detD/(16*C**5))*I1);VV=s.simplify(Zinv*3/s.gamma(s.Rational(9,2))*bD/C**5);assert s.simplify(U2-8*s.pi*detD/C**7)==0;assert s.simplify(VV-s.Rational(128,35)*bD/C**5)==0
c=1/(32*s.pi**3);budget_squared=s.factor((c*c*U2**3/VV**4).subs(C,4/t+12));claimed_squared=(s.Rational(1225,32768)**2*2/s.pi**3)*t*(3*t+1)**2/(144**2*(t+3)**2);assert s.simplify(budget_squared-claimed_squared)==0
# rho<1 gives t²<11/7<81/49=(9/7)²; all denominators positive.
assert s.simplify(s.Rational(11,7)-(12-rho)/(12-5*rho)-48*(rho-1)/(7*(5*rho-12)))==0
assert F(11,7)<F(9,7)**2 and F(9,7)<F(8,7)**2
assert s.simplify(s.diff((3*t+1)/(t+3),t)-8/(t+3)**2)==0
assert F(2,3)<F(5,6)**2
point_upper=F(6125,589824)*F(17,1890);assert point_upper<F(1,10000)
patch_factor=F(1009,991)**9*F(1000,999)**4;assert patch_factor<F(6,5);patch_upper=F(1,10000)*patch_factor;assert patch_upper<F(3,25000)<F(15198,10**8)
paths=[Path(__file__),P/'round009/profile_barrier_registration.md',Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/astra_sqrt030/report.md'),Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper-audit-cont/cube-dim-5-bounds/upper_alpha1_029/report.md')]
out={'status':'PASS_exact_continuous_profile_family_obstruction_pending_independent_review','family':'alpha1 Z=a0*(1-rho*(W_x+2)/12)^(-3/2), all rho<1,a0>=0 including spatial choices','point':[1,1,1,0,0],'direction_counts':counts,'D':[[str(a) for a in row] for row in D.tolist()],'detD':str(detD),'bD':str(bD),'point_budget_upper_rational':str(point_upper),'strict_point_cap':'1/10000','eps_l2':'1/1000','patch_factor':str(patch_factor),'strict_face_patch_cap':'3/25000','explicit_positive_measure_face_patch':'x1=1; x2,x3 in [3999/4000,1]; x4,x5 in [0,1/4000]','limitation':'Upper cap concerns the Young-derived lower-function budget, not actual G; actual G(x0)=infinity. Full family uniform cap extends to positive-measure patch by relative lambda and Loewner estimates in accompanying proof. No new C5 bound.','checks':{'symbolic_16_direction_D_and_cofactor':True,'Mellin_constants':True,'budget_squared_identity':True,'uniform_parameter_comparisons_exact':True,'positive_measure_patch_margin_exact':True},'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began,'sympy':s.__version__}}
(P/'round009/profile_barrier.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','point_budget_upper_rational','strict_point_cap','strict_face_patch_cap','runtime']},indent=2))
