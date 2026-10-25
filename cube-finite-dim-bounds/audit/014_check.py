#!/usr/bin/env python3
"""Independent exact Poisson/Green monotonicity and angular-stratum identities."""
from pathlib import Path
import hashlib,json
import sympy as s
R,r=s.symbols('R r',positive=True);z=s.symbols('z',real=True)
D=R**2+r**2-2*R*r*z
assert s.expand(D-((R-r)**2+2*R*r*(1-z)))==0
N3=(R+2*r)*D-3*r*(R+r)*(r-R*z)
F3=(R-r)*(R**2+3*R*r+r**2+R*r*z)
assert s.expand(N3-F3)==0
assert s.expand((R**2+3*R*r+r**2+R*r*z)-(R+r)**2-R*r*(1+z))==0
q3=r*(R+r)/D**s.Rational(3,2)
assert s.simplify(s.diff(q3,r)-F3/D**s.Rational(5,2))==0
q4=r*r/D**2;F4=2*r*(R**2-r**2)
assert s.simplify(s.diff(q4,r)-F4/D**3)==0
# Sphere stationary strata: equal multiplicities give a continuum, constant value.
a,T,d,m=s.symbols('a T d m',real=True);b=T-a;e=-3*d*T/2
g=s.expand(m*(d*(a**3+b**3)+e*(a*a+b*b)))
assert s.simplify(g+m*d*T**3/2)==0
rows=[]
B4,B6=s.symbols('B4 B6',real=True)
for n in (3,4):
 u=s.symbols('u0:'+str(n),positive=True);S=sum(u);p=s.Rational(n-2,2)
 alpha=s.Rational(15,n+8);beta=s.Rational(30,(n+8)*(n+4))
 h4=sum(x*x for x in u)-s.Rational(3,n+2)*S*S
 h6=sum(x**3 for x in u)-alpha*S*sum(x*x for x in u)+beta*S**3
 f=S**(-p)-B4*h4-B6*h6
 expected=(u[0]-u[1])*(-2*B4-3*B6*(u[0]+u[1])+2*alpha*B6*S)
 assert s.simplify(s.diff(f,u[0])-s.diff(f,u[1])-expected)==0
 rows.append({'n':n,'free_coordinate_difference':str(expected),'scope':'at most two distinct free positive squared coordinates for B6 != 0; degeneracies still need value/root treatment'})
out={'status':'exact Poisson and stratum identities PASS','engine':'SymPy '+s.__version__,'legal_parameter_domain':'n3/4, c>0, h harmonic on ball and continuous closure, F=c*r^(2-n)-h strictly negative on full sphere','distance_square':str(D),'n3_q':'r*(R+r)/D^(3/2)','n3_derivative_numerator':str(F3),'n3_numerator_lower':'(R-r)*(R+r)^2 > 0 for 0<r<R, -1<=z<=1','n4_q':'r^2/D^2','n4_derivative_numerator':str(F4),'equal_multiplicity_angular_value':str(s.factor(g)),'cube_rows':rows,'scope_limit':'no cost, no asserted cube corner minimum; unguarded directional multiroot examples not ruled out','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('014_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
