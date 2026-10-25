"""Conditional homogeneous mass algebra; does not certify a Green kernel."""
import json, hashlib
from pathlib import Path
import sympy as s
p=Path(__file__).resolve().parents[1]
z,c,g,J=s.symbols('z c g J', positive=True)
rows=[]
for n in (3,4):
 q=s.Rational(1,n-1); radial=s.integrate(g-c*z**(n-1),(z,0,(g/c)**q))
 expected=s.Rational(n-1,n)*g**(1+q)*c**(-q)
 assert s.simplify(radial-expected)==0
 mass=s.Rational(n-1,n)*J*c**(-q)/(2**n*(g-c))
 assert s.simplify(s.diff(mass,c).subs(c,g/n))==0
 # J=2*n*2**(n-1)*g**(1+q)*E for identical symmetric faces.
 E=s.symbols('E', positive=True)
 best=s.simplify(mass.subs(c,g/n).subs(J,n*2**n*g**(1+q)*E))
 assert s.simplify(best-n**(1+q)*E)==0
 logder=s.simplify(s.diff(mass,c)/mass)
 assert s.simplify(logder-(n*c-g)/((n-1)*c*(g-c)))==0
 rows.append({'n':n,'radial_integral':str(expected),'log_derivative':str(logder),'optimal_c':'g/'+str(n),'cost':str(best),'required_face_average':str(2**n/n**(1+q)),'required_face_average_numeric':str(s.N(2**n/n**(1+q),30))})
out={'status':'exact conditional algebra PASS; no actual G/g/J or upper bound','assumptions':['G positive homogeneous degree 1-n','g positive full cube infimum','J finite and equals face integral','full sign/permutation symmetry for single-face average','kernel distribution comparison separately required'],'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(p/'results/round009_mass_interface.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
