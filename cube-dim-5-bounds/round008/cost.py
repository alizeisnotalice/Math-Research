from pathlib import Path
from fractions import Fraction as Q
from flint import arb,ctx
import json,hashlib,time,os
p=Path(__file__).resolve().parents[1];ctx.prec=256;began=time.monotonic()
f=p/'astra_upper009/verification.json';v=json.loads(f.read_text());row=next(x for x in v['candidates'] if x['eta']=='27/250');g=Q(row['global_g_exact_rational']);gprev=Q(json.loads((p/'astra_upper007/verification.json').read_text())['global_g_exact_rational'])
assert g>Q(2821,50000000)>gprev
Ipath=p/'round006/joint_cost.json';I=arb(json.loads(Ipath.read_text())['I_ball']);ga=arb(g.numerator)/g.denominator;cost=(arb(10)/3)**(arb(10)/7)*I/ga**(arb(10)/7);assert cost>32
out={'status':'complete_auxiliary_g_qualified_no_C5_improvement','exact_g':str(g),'simple_g':'2821/50000000','exact_new_minus_previous_g':str(g-gprev),'I_source':str(Ipath.relative_to(p)),'cost_ball':cost.str(70,more=True),'C5_upper':'32','C5_lower_certificate':'lower_round003/a2_equals_a3_results.json','source_sha256':{str(x.relative_to(p)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [f,Ipath,p/'audit_g015/report.md']},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began,'arb_bits':ctx.prec}}
(p/'round008/cost.json').write_text(json.dumps(out,indent=2)+'\n');print('new model cost:',cost)
