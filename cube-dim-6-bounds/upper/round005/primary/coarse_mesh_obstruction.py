from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import json,hashlib
ctx.prec=256
r=Path(__file__).resolve().parents[3]
mesh=r/'upper/round005/cover/simplex_cover_audit.json';m=json.loads(mesh.read_text());D=F(m['leaf_exact_variance_min']);assert D==F(25,192)
# The spectral second-moment lower bound comes from integration over t>=T.
spectral=json.loads((r/'upper/round004/primary/spectral_moment.json').read_text());assert arb(spectral['uniform_g_lower'])>arb(6)/10**8;assert arb(spectral['T'])>arb(1)/100
# Thus H_e >=L on the cube for all e<=1/100.
L=arb(6)/10**8;U=arb(23)/10**6
point=r/'upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json';p=json.loads(point.read_text())
# Only the verified single-point upper bound on essential g is needed.
value=next(v for k,v in p.items() if k.startswith('G_arb'))
assert arb(value)<U
K=(arb(6)/7)*(arb(1)/243)**(arb(1)/6)*arb.pi()**(-arb(11)/6)
e=arb(1)/100;limit=arb(16)/3*U
coarse_H_lower=L*(arb(D.numerator)/D.denominator/(32*e**2)).exp()
tail_at_e=K*e**(arb(7)/6)
assert coarse_H_lower>limit and tail_at_e>limit
out={'status':'PASS_fixed_coarse_mesh_budget_cannot_beat64_for_any_epsilon','epsilon_split':'1/100','allowed_norm_budget_at_g_upper':str(limit),'coarse_H_budget_lower_for_epsilon_le_split':str(coarse_H_lower),'tail_budget_lower_for_epsilon_ge_split':str(tail_at_e),'scope':'excludes this 120-leaf variance-maximum upper-bound method only, not actual I, model cost, or C6','mesh_sha256':hashlib.sha256(mesh.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
