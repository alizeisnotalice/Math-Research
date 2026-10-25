from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import json,hashlib,time,os
p=Path(__file__).resolve().parents[1];ctx.prec=256;began=time.monotonic()
def ar(x):
 x=F(x);return arb(x.numerator)/x.denominator
files=[p/'round008/hermite_linear.json',p/'upper_I024/results.json',p/'audit_g016/report.md',p/'audit_I026/report.md',p/'round007/cost_K3.json',p/'round006/joint_cost.json']
gdata=json.loads(files[0].read_text());Idata=json.loads(files[1].read_text());g=F(gdata['exact_g']);lo=F(Idata['I_optimal_budget_interval']['lower_fraction']);hi=F(Idata['I_optimal_budget_interval']['upper_fraction']);simpleI=F('1620961/125000000000');simpleg=F(gdata['simple_g'])
assert hi<simpleI and g>simpleg
q=ar(F(10,7));factor=ar(F(10,3))**q;cost=factor*ar(hi)/ar(g)**q;safe=factor*ar(simpleI)/ar(simpleg)**q
old=json.loads(files[4].read_text());oldg=F(old['g_exact_rational']);oldI=arb(json.loads(files[5].read_text())['I_ball'])
assert g>oldg and ar(hi)<oldI and cost>32
out={'status':'complete_qualified_auxiliary_g_I_no_C5_improvement','g_exact':str(g),'g_simple':str(simpleg),'I_formula':'upper_I024/report.md; 2Gamma7/4*m(3/20)^5*Z^(-3/20), same orbit interface','I_formula_exact_interval':[str(lo),str(hi)],'I_simple_rational_upper':str(simpleI),'I_input_used':'exact_saved_upper_root_endpoint_of_complete_formula_not_display_decimal','joint_budget_upper_ball':cost.str(70,more=True),'conservative_rational_input_budget_ball':safe.str(70,more=True),'g_new_minus_round007_exact':str(g-oldg),'I_new_budget_upper_minus_round006_ball':(ar(hi)-oldI).str(70,more=True),'joint_budget_new_minus_round007_ball':(cost-arb(old['cost_ball'])).str(70,more=True),'g_needed_if_using_new_I_budget_for_cost32':(ar(F(10,3))*(ar(hi)/32)**ar(F(7,10))).str(70,more=True),'C5_upper':'32','C5_lower':'complete C_cert; lower_round003/a2_equals_a3_results.json','continuation_bound_differences':{'lower':'0','upper':'0'},'post_dual_rounds':0,'sources_sha256':{str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began,'arb_bits':ctx.prec}}
(p/'round008/joint_certified.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','joint_budget_upper_ball','conservative_rational_input_budget_ball','g_needed_if_using_new_I_budget_for_cost32','joint_budget_new_minus_round007_ball']},indent=2))
