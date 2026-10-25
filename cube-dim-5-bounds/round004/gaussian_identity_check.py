import sympy as s,json,hashlib
from pathlib import Path
z,a=s.symbols('z a',positive=True)
# z=a/(2t²), t=(a/(2z))^(1/2)
t=s.sqrt(a/(2*z))
trans=s.simplify(t**(-5)*s.exp(-z)*(-s.diff(t,z)))
value=s.integrate(trans,(z,0,s.oo))
assert s.simplify(value-2/a**2)==0
out={'status':'symbolic_change_of_variable_check_only','transformed_integrand':str(trans),'integral':str(value),'sympy':s.__version__,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
