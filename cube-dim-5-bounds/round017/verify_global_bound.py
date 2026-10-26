"""Portable exact bound check; analytical and point certificates stay separate."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent.parent
def read(n):return json.loads((P/n).read_bytes())
a=read('round017/primary_global_acceptance.json')
b=read('audit_global148/global_actual/results.json')
old=read('round016/primary_upper139_acceptance.json')
assert (b['roots'],b['passed'],b['splits'],b['pending'])==(52360,111201,58841,0)
assert b['epsilon_leaf_counts']=={'1/64':98642,'1/128':12559}
for n,k in [('leaf_lineage.jsonl.gz','leaf_lineage_sha256'),('root_coverage.jsonl.gz','root_coverage_sha256')]:
    assert hashlib.sha256((P/'audit_global148/global_actual'/n).read_bytes()).hexdigest()==b[k]
m=Q(b['minimum_saved_log_margin']);g=Q(11,100000)*(1+m)
assert m>0 and g==Q(a['qualified_g'])==Q(b['optional_rational_G_cube_floor_for_primary'])
assert a['minimum_witness']==b['minimum_margin_witnesses'][0]
I=Q(a['complete_I_upper']);assert I==Q(old['complete_I_upper'])
c4=3125*I**4/g**5;assert c4==Q(a['best_upper']['fourth_power'])<Q(old['best_upper']['fourth_power'])<32**4
lo,hi=map(Q,a['best_upper']['interval']);assert 0<lo<hi and lo**4<=c4<=hi**4
simple=Q(a['simple_rational_upper']);assert hi<=simple<32
assert a['general_L1_all_positive_common_radii'] and a['lower_changed'] is False
print('PASS complete-g source/coverage hashes and exact new C5 upper; no H/log/corner replay')
