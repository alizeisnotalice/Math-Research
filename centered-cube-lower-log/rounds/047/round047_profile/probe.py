from pathlib import Path
import json,math,hashlib
from scipy.special import log_ndtr
from scipy.optimize import minimize_scalar
H=Path(__file__).resolve().parent
rows=[]
for v in [8,16,32,64,128,256]:
 def logW(h):return h+float(log_ndtr(-(h+v/2)/math.sqrt(v)))
 r=minimize_scalar(lambda h:-logW(h),bounds=(-v,2*v),method='bounded',options={'xatol':1e-11})
 rows.append(dict(v=v,best_h=float(r.x),best_W=math.exp(logW(r.x)),sqrt_2piv_times_best_W=math.sqrt(2*math.pi*v)*math.exp(logW(r.x)),baseline_h=v/2-2,baseline_W=math.exp(logW(v/2-2)),low_h_minus_point1=math.exp(logW(-.1)),constant_input_weak_at_h_minus_point1=math.exp(-.1)))
out=dict(rows=rows,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Exact Gaussian tail formula evaluated numerically; no HL transfer or outward certificate.')
(H/'PROFILE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
