from pathlib import Path
from fractions import Fraction as Q
import json,time
from rational_lib import F24,E24,grid,encode
p=Path(__file__).resolve().parents[2];r=json.loads((p/'upper/round024/parameter_numeric/result.json').read_text());start=time.perf_counter()
a=Q(r['candidate']);lo=Q(r['root_bracket_lower']);hi=Q(r['root_bracket_upper']);assert hi-lo==Q(1,2**79) and a==(lo+hi)/2
el=E24(lo);eh=E24(hi);assert el[1]<0<eh[0]
f=F24(a);f36=F24(Q(36));U=Q(r['candidate_grid_upper']);U36=Q(r['F36_grid'])
assert f[1]<U and grid(f[1])==U and f36[1]<U36 and grid(f36[1])==U36 and f[1]<f36[0]
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
old=dy(r['R023_old_interval']['lower']);assert U<old
print(json.dumps(dict(status='PASS primary independent Fraction log/sqrt certificates',source_candidate=str(a),root_bracket=[str(lo),str(hi)],endpoint_E=[encode(el),encode(eh)],F_candidate=encode(f),F36=encode(f36),strict_C7_upper=str(U),scope='candidate upper not family minimum-value enclosure; final root endpoints independently checked, no80step replay',wall_seconds=time.perf_counter()-start),indent=2))
