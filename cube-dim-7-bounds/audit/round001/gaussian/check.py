#!/usr/bin/env python3
"""Independent 768-bit full-I expression audit; no Monte Carlo or author imports."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time,platform
import sympy as sp
import flint
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
NAMES=['upper/round001/gaussian/proof.md','upper/round001/gaussian/check.py','upper/round001/gaussian/check-results.json','upper/round001/gaussian/results.json','upper/round001/gaussian/sources-hash.json','src/check_exact_budget_round001.py','audit/round001/exact_budget_results.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inputs():return {x:sha(ROOT/x)for x in NAMES}
def sealed():return {p.name:sha(p)for p in sorted((HERE.parent/'theory').iterdir())if p.is_file()}
start=time.perf_counter(); before=inputs(); seal_before=sealed()
manifest=json.loads((HERE.parent/'theory'/'manifest.json').read_text())
for n,v in manifest['files'].items(): assert seal_before[n]==v['sha256']
provenance=json.loads((ROOT/NAMES[4]).read_text())
for name in ('proof.md','check.py','check-results.json','results.json'):
    assert sha(ROOT/'upper/round001/gaussian'/name)==provenance['files'][name]['sha256']
locked=F('74319631963157785256240159489393050185362460367580814299903508439387034491195/7588550360256754183279148073529370729071901715047420004889892225542594864082845696')
parent=json.loads((ROOT/NAMES[6]).read_text()); assert F(parent['I_certified_upper'])==locked
assert sha(ROOT/NAMES[5])==parent['checker_sha256']
# Exact analytical constants, independently arranged before numerical evaluation.
n=sp.Integer(7); q=sp.Rational(7,6)
C=sp.gamma(3)/(4*sp.pi**sp.Rational(7,2))
H=2*sp.pi**sp.Rational(7,2)/sp.gamma(sp.Rational(7,2))
assert sp.simplify(C-1/(2*sp.pi**sp.Rational(7,2)))==0
assert sp.simplify(H-16*sp.pi**3/15)==0
assert -sp.Rational(1,2)+1/(2*q)==-sp.Rational(1,14)
assert sp.simplify(C**q*H/896-1/(105*2**sp.Rational(25,6)*sp.pi**sp.Rational(13,12)))==0
assert 2*7*2**6==896
# Dirichlet beta factor cancels radial Gamma32/Gamma(65/2), leaving2/sqrt64=1/4.
assert sp.simplify(2/sp.sqrt(64))==sp.Rational(1,4)
assert 64*sp.Rational(1,128)==sp.Rational(1,2)
assert 32+sp.Rational(1,2)==sp.Rational(65,2)
ctx.prec=768
pi=arb.pi()
# Independent log-Gamma evaluation; boundary coefficient uses original C,H rather than simplified K.
logratio=(arb(65)/128).lgamma()-(arb(1)/2).lgamma()
mg=(64*logratio).exp()/4
ca=1/(2*pi**(arb(7)/2)); ha=16*pi**3/15
ival=ca**(arb(7)/6)*ha*mg**(arb(7)/6)/896
simple=mg**(arb(7)/6)/(105*arb(2)**(arb(25)/6)*pi**(arb(13)/12))
assert ival.overlaps(simple)
def exact_interval(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
mi=exact_interval(mg); ii=exact_interval(ival)
il,iu=map(F,ii)
assert 0<il<iu<locked
# Also check the legacy256 endpoint's own conditional comparison, not only a tighter replacement.
g=F(1,179000)
cost_new=F(7**7)*iu**6/g**7
cost_old=F(7**7)*locked**6/g**7
assert 0<cost_new<cost_old<128**6
assert F(parent['conditional_cost_sixth_power_exact'])==cost_old
assert F(parent['conditional_gap_below_128_sixth_power'])==128**6-cost_old
assert parent['global_g_target_proved'] is False and parent['new_C7_upper_bound'] is False
assert inputs()==before and sealed()==seal_before
out={'status':'PASS independent complete-I analytic supplement only','bits':768,'environment':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__},'fixed_expression':'(Gamma(65/128)/Gamma(1/2))^(224/3)/(4^(7/6)*105*2^(25/6)*pi^(13/12))','evaluation':'loggamma ratio; C^(7/6)*H*M^(7/6)/896 from unsimplified boundary formula','M_Gamma_outward_rational_interval':mi,'I_expression_outward_rational_interval':ii,'I768_upper_strictly_less_than_locked_I256_upper':True,'locked_I256_upper':str(locked),'locked_minus_I768_upper_exact':str(locked-iu),'conditional_g_target':str(g),'conditional_legacy_cost_sixth_power_exact':str(cost_old),'conditional_tighter_cost_sixth_power_exact':str(cost_new),'legacy_cost6_strictly_below_1286':True,'legacy_exact_gap':str(128**6-cost_old),'global_g_target_proved':False,'new_C7_upper_bound':False,'MC_status':'candidate diagnostics only; not read/replayed or certified here','new_random_samples':0,'new_Green_nodes':0,'source_hashes':before,'sealed_theory_hashes':seal_before,'source_before_after_equal':True,'sealed_theory_before_after_equal':True,'elapsed_seconds':time.perf_counter()-start}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ('status','bits','I768_upper_strictly_less_than_locked_I256_upper','legacy_cost6_strictly_below_1286','global_g_target_proved','new_C7_upper_bound','elapsed_seconds')},indent=2))
