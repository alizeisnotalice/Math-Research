"""Finite new-identity check; does not rerun dimension-one baselines."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import sympy as s
x=s.symbols('delta',positive=True)
tube=s.expand(sum(s.binomial(4,k)*(2*x)**(4-k)*(k*(2*x)**(k-1)-(k-1)*(2*x)**k) for k in range(1,5)))
assert tube==256*x**3-272*x**4
V=F(14455,191102976)
cp=F(3,2)*V
cf=4096*cp
local=3*256*14**4*F(125,6)
assert cp==F(14455,127401984) and cf==F(14455,31104)
assert local==614656000
out={'status':'PASS_finite_algebra_only','tube_polynomial':str(tube),'Cperp_times_pi4':str(cp),'Cface_times_pi5':str(cf),'local_d2_density_coefficient_times_pi5':str(local),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
raw=json.dumps(out,indent=2)+'\n'
Path(__file__).with_suffix('.json').write_text(raw)
print(raw)
