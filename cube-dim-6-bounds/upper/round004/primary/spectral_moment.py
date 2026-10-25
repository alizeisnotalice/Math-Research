from fractions import Fraction as F
from pathlib import Path
import json,math,hashlib
from flint import arb,ctx
ctx.prec=256
root=Path(__file__).resolve().parents[3];src=root/'upper/round002/fan/exact_cells_results.json'
d=json.loads(src.read_text());s0=F(0);s2=F(0)
for c in d['cells']:
 r=[[F(x) for x in v] for v in c['rays_a0_equals_one']];q=list(map(F,c['rates']));alpha=F(c['determinant'])/math.prod(q)
 z=[[r[j][i]/q[j] for j in range(6)] for i in range(6)]
 mean_norm2=sum(sum(v)**2+sum(t*t for t in v) for v in z)
 s0+=alpha;s2+=alpha*mean_norm2
E=s2/s0;A=lambda f:arb(f.numerator)/f.denominator
C0=A(s0*46080)/(2*arb.pi())**6
# cos>=1-u²/2; p_t(x)>=C0 t^-6(1-E/(2t²)) on unit cube.
# integrate from T=sqrt(E/2) to infinity gives positive lower g.
T=(A(E)/2).sqrt();g=C0*(1/(5*T**5)-A(E)/(14*T**7))
out={'status':'exact_spectral_second_moment_cube_lower_bound','E_norm2':str(E),'C0':str(C0),'T':str(T),'uniform_g_lower':str(g),'scope':'valid uniform G lower bound, not an I upper bound or C6 bound','fan_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
