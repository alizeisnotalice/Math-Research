from pathlib import Path
from fractions import Fraction as Q
import json,time
from cap_lib import cap,base
p=Path(__file__).resolve().parents[2];r=json.loads((p/'upper/round027/cube_cap_scalar/result.json').read_text());start=time.perf_counter();rows=[]
for t in [Q(4),Q(5),Q(11,2),Q(6),Q(8)]:
 B,F,c,k,terms=cap(t);rows.append(dict(t=str(t),Bcap=base.encode(B),F=base.encode(F),rational_constant=str(c),log_coefficient=str(k),signed_terms=[str(v)for v in terms]))
w=min(rows,key=lambda z:Q(z['F'][1]));U=base.grid(Q(w['F'][1]));assert w['t']=='6' and str(U)==r['winner']['grid_upper']
assert all(Q(w['F'][1])<Q(z['F'][0])for z in rows if z is not w)
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
assert U<dy(r['R026_old_candidate_interval']['lower'])
print(json.dumps(dict(status='PASS independent exact signed-coefficient cap certificate',winner=w,strict_C7_upper=str(U),rows=rows,continuous_t_optimum_claim=False,wall_seconds=time.perf_counter()-start),indent=2))
