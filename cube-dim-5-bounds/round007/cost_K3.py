from pathlib import Path
from fractions import Fraction as Q
from flint import arb,ctx
import json,hashlib
p=Path(__file__).resolve().parents[1];ctx.prec=256
f=p/'astra_upper007/verification.json';v=json.loads(f.read_text());g=Q(v['global_g_exact_rational']);g5=max(Q(x['global_g_exact_rational']) for x in json.loads((p/'astra_upper005/verification.json').read_text())['candidates'])
assert g>Q(v['global_g_simple_rational'])>g5
Ipath=p/'round006/joint_cost.json';I=arb(json.loads(Ipath.read_text())['I_ball']);ga=arb(g.numerator)/g.denominator;cost=(arb(10)/3)**(arb(10)/7)*I/ga**(arb(10)/7);assert cost>32
out={'status':'new_strict_auxiliary_g_K3_no_C5_breakthrough','g_exact_rational':str(g),'g_simple_rational':v['global_g_simple_rational'],'new_minus_005_g_exact_rational':str(g-g5),'cost_ball':cost.str(70,more=True),'C5_upper':'32','source_hashes':{str(x.relative_to(p)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [f,Ipath,p/'audit_g014/report.md']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(p/'round007/cost_K3.json').write_text(json.dumps(out,indent=2)+'\n');print('model cost',cost)
