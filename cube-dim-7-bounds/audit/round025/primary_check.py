from pathlib import Path
from fractions import Fraction as Q
import importlib.util,json,time
p=Path(__file__).resolve().parents[2];spec=importlib.util.spec_from_file_location('rational_lib',p/'audit/round024/rational_lib.py');lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(lib)
r=json.loads((p/'upper/round025/buffer_scalar/result.json').read_text());start=time.perf_counter();rows=[]
for t in [Q(5),Q(11,2),Q(6)]:
 f=lib.F25(t);rows.append(dict(t=str(t),interval=lib.encode(f)))
w=min(rows,key=lambda z:Q(z['interval'][1]));U=lib.grid(Q(w['interval'][1]));assert w['t']=='11/2' and str(U)==r['winner']['grid_upper']
assert all(Q(w['interval'][1])<Q(z['interval'][0])for z in rows if z is not w)
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
assert U<dy(r['R024_candidate_interval']['lower'])
l,u=lib.log(Q(5));A=Q(6,5)**7;tau=Q(20,3);cal=((1+A/tau)*(tau+1+7*l),(1+A/tau)*(tau+1+7*u));assert cal[1]<Q(2233,75)<30
assert sum((Q(5,3)**j/Q(__import__('math').factorial(j))for j in range(5)),Q())==Q(10009,1944)>5
print(json.dumps(dict(status='PASS primary independent Fraction three-row and calibration certificate',winner=w,strict_C7_upper=str(U),rows=rows,calibration=lib.encode(cal),continuous_t_optimality_claimed=False,wall_seconds=time.perf_counter()-start),indent=2))
