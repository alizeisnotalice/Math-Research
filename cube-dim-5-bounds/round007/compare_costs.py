from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import json,hashlib,time
p=Path(__file__).resolve().parents[1];start=time.monotonic();ctx.prec=256
oldpath=p/'round006/joint_cost.json';old=json.loads(oldpath.read_text());I=arb(old['I_ball']);oldg=arb(old['g_ball']);oldU=arb(old['auxiliary_cost_ball']);rows=[]
for batch in ['astra_upper004','astra_upper005']:
 f=p/batch/'verification.json';d=json.loads(f.read_text());c=max(d['candidates'],key=lambda x:F(x['global_g_exact_rational']));g=F(c['global_g_exact_rational']);ga=arb(g.numerator)/g.denominator;U=(arb(10)/3)**(arb(10)/7)*I/ga**(arb(10)/7)
 assert ga>oldg and U>32 and U<oldU
 rows.append({'batch':batch,'g_exact_rational':str(g),'g_simple_rational':c['global_g_simple_rational'],'new_minus_round006_g':(ga-oldg).str(65,more=True),'complete_model_cost':U.str(65,more=True),'new_minus_round006_cost':(U-oldU).str(65,more=True),'input_sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
out={'status':'auxiliary_g_improvements_still_dominated_by_C5_upper32','I_input_sha256':hashlib.sha256(oldpath.read_bytes()).hexdigest(),'rows':rows,'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(p/'round007/compare_costs.json').write_text(json.dumps(out,indent=2)+'\n');print([(x['batch'],x['complete_model_cost']) for x in rows])
