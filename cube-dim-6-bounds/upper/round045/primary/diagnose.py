from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time,os
P=Path(__file__).resolve().parent
reg=json.loads((P/'registration.json').read_text());s=Path(reg['input'])
assert hashlib.sha256(s.read_bytes()).hexdigest()==reg['sha256']
tic=time.process_time();d=json.loads(s.read_text())
def midpoint(v):return Q(int(v['lower_mantissa'])+int(v['upper_mantissa']),2)*Q(2)**int(v['common_binary_exponent'])
g=list(map(midpoint,d['gradient']));blocked=max(map(abs,g[1:4]));free=max(map(abs,g[4:6]));allfree=max(map(abs,g[1:6]));assert blocked>0 and free>0
out={'pid':os.getpid(),'arithmetic':'Fraction exact saved dyadic midpoints; display only float','source_sha256':reg['sha256'],'gradient_midpoints':list(map(str,g)),'gradient_displays':list(map(float,g)),'blocked_normalization_max':str(blocked),'two_free_normalization_max':str(free),'free_over_old_normalization_exact':str(free/allfree),'free_over_old_normalization_display':float(free/allfree),'native_calls':0,'new_candidates':0,'C6_change':False,'CPU_seconds':time.process_time()-tic}
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in ('gradient_midpoints',)},indent=2))
