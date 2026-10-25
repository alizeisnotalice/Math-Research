from pathlib import Path
from flint import arb,ctx
import json,hashlib
p=Path(__file__).resolve().parents[1];ctx.prec=256
files=['upper_g_round010/tilted_gaussian_bound.json','upper/round011_shape_integral_results.json','round006/joint_cost.json'];a,b,c=[json.loads((p/f).read_text()) for f in files]
oldg=arb(a['green_bound_arb']);oldI=arb(b['arb_256bit']['complete_I_upper']);newg=arb(c['g_ball']);newI=arb(c['I_ball']);dg=newg-oldg;di=newI-oldI
assert dg>0 and di<0
oldcost=(arb(10)/3)**(arb(10)/7)*oldI/oldg**(arb(10)/7);newcost=arb(c['auxiliary_cost_ball']);delta=newcost-oldcost;assert delta<0 and newcost>32
out={'status':'strict_complete_auxiliary_certificate_improvements_no_C5_breakthrough','new_minus_old_g':dg.str(65,more=True),'new_minus_old_I':di.str(65,more=True),'preceding_full_formula_cost':oldcost.str(65,more=True),'new_minus_old_cost':delta.str(65,more=True),'source_hashes':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in files},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(p/'round006/compare_aux.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
