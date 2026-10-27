#!/usr/bin/env python3
"""Exact comparison against complete old costs; no search or old experiment replay."""
import sys
if hasattr(sys,'set_int_max_str_digits'): sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,csv
p=Path(__file__).resolve().parent
project=p.parent
old=json.loads((p/'literature/internal_baselines.json').read_text())
base={x['dimension']:x for x in old['rows']}
rows=json.loads((p/'numeric/table.json').read_text())
old7=json.loads((project/'audit/round028/primary_result.json').read_text())
lo7,hi7=map(Q,old7['F'])
assert hi7<Q(base[5]['upper']['exact_interval'][0])
assert hi7<64
# Machin identity pi=16 atan(1/5)-4 atan(1/239), alternating-series bounds.
pi_lower=16*(Q(1,5)-Q(1,5)**3/3)-4*Q(1,239)
assert pi_lower>Q(157,50)
c4=Q(6244212653911047061276408,1944624818861484527587890625)
new={x['n']:x for x in rows}
assert Q(new[4]['F_upper'])**3 < Q(new[4]['U'])**3 < c4*pi_lower**12
for n in (5,6): assert Q(new[n]['F_upper'])<Q(new[n]['U'])<lo7
assert Q(1547906585017734375000000000,24707828464688747938727861)<8**2
assert Q(base[2]['upper']['enclosure'][1])<4
assert new[7]['t']==old7['exact_t'] and Q(new[7]['U'])==Q(old7['strict_C7_upper'])
assert new[8]['t']==base[8]['upper']['t']
assert Q(base[8]['upper']['strict_rational_upper'])<Q(new[8]['U'])
result=[]
for r in rows:
 n=r['n'];candidate=r['U']
 if n==2:
  before='359091/100000';after=before;status='retain stronger existing two-dimensional result'
 elif n==3:
  before='1979/250';after=before;status='retain stronger existing three-dimensional result'
 elif n==4:
  before='14371/1000';after=candidate;status='strict improvement vs complete old pi/root cost'
 elif n in (5,6):
  before='4910992105193/200000000000';after=candidate;status='strict improvement vs complete inherited seven-dimensional cost'
 elif n==7:
  before='4910992105193/200000000000';after=before;status='same exact parameter and bound; no new C7 bound'
 elif n==8:
  before=base[8]['upper']['strict_rational_upper'];after=before;status='same exact parameter; retain tighter old rational enclosure'
 else:
  before=str(2**n);after=candidate;status='strict improvement vs verified generic baseline only; not a best-known claim'
  assert Q(r['F_upper'])<2**n
 result.append({'n':n,'t':r['t'],'candidate_strict_upper':candidate,'old_effective_display_bound':before,'selected_usable_bound':after,'comparison':status,'lower_source':f'literature/internal_baselines.json#/rows/{n-2}/lower' if n<=8 else 'dimension monotonicity from dimension8; see literature/internal_baselines.json#/rows/6/lower','global_family_status':'endpoint infimum=2^n, unattained' if n<=3 else 'an interior global minimum exists; this candidate not certified global; uniqueness unproved','finite_a':'2^256; explicit rational tau in numeric/table.json; exact checker PASS','public_scope':'Tomar v4 (2026-10-22), upper 3.59123; public preprint; not rerun' if n==2 else 'generic 2^n verified; finite literature search only; no exhaustive strongest-bound claim','complete_candidate_certificate':f'numeric/table.json row n={n}'})
receipt={'status':'PASS exact old/new complete-cost comparisons','pi_lower':str(pi_lower),'machin_lower_gt_157_over_50':True,'n4_new_upper_cube_below_old_cost_cube_lower':True,'n5_n6_new_upper_below_old_full_cost_lower':True,'n7_same_exact_parameter':True,'n8_same_exact_parameter_old_enclosure_retained':True,'no_new_lower_computation':True,'snapshot':old['snapshot'],'source_sha256':{str(f.relative_to(project)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [p/'numeric/table.json',p/'literature/internal_baselines.json',project/'audit/round028/primary_result.json']},'rows':result}
(p/'dimension_comparison.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
with (p/'dimension_comparison.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=result[0].keys(),lineterminator='\n');w.writeheader();w.writerows(result)
print('PASS complete-cost comparison; n4,n5,n6 improve; n7 unchanged; n8 old retained; no lower work')
