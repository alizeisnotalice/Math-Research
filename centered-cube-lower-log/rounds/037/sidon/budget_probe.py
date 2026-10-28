"""Exact scalar gates and finite far-correlation diagnostic before proof."""
from pathlib import Path
from fractions import Fraction as F
import json,math
rows=[]
for b in [2,3,4,5,6,7,8,12,13,14,16]:
 B=2**b;K=B**3-1;D=(B-1)*(B-2)
 H=F(1,12)+F(1,6*B*B);E=F(1,48)+F(1,6*B*B)
 Q=math.sqrt(float(H*E))+float(E)
 C=Q+1/(12*B*B)+1/(2*K)
 rows.append({'b':b,'B':B,'K':str(K),'D':D,
  'primitive_count_sufficient_exact':K**3>(6*D)**4,
  'two_segment_candidate_cap':Q,
  'far_upper_limit_with_quantization_and_zero_state_only':(4*C+34/261)/(100/261),
  'asymptotic_far_upper_exact':'397/400','asymptotic_far_lower_exact':'-27/40'})
out={'scope':'Exact primitive sufficient inequality and ordinary-float source far-bound diagnostics only; no complete correlation or C_n theorem.', 'rows':rows}
Path(__file__).with_name('budget_results.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:print(r['B'],r['primitive_count_sufficient_exact'],r['far_upper_limit_with_quantization_and_zero_state_only'])
