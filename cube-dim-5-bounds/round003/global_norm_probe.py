import mpmath as m
from pathlib import Path
import json,hashlib,time
m.mp.dps=50
start=time.monotonic()
P=m.mpf(35509)/(764411904*m.pi**5)
F=m.mpf(14455)/(31104*m.pi**5)
a=m.sqrt(P/F)
def integrand(t):
 return min(F*t**-3,P*t**-5)**m.mpf('.2')*(t/(m.pi*(1+256*t*t)))**m.mpf('.8')
value=m.quad(integrand,[0,a])+m.quad(integrand,[a,1,m.inf])
g=m.mpf('0.000146450542661776629693871130288973764')
out={'status':'nonrigorous_budget_diagnostic_not_C5_bound','crossover':str(a),'global_Minkowski_norm_proxy':str(value),'best_sample_g_proxy':str(g),'cost_proxy':str((5*value/g)**m.mpf('1.25')),'seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
raw=json.dumps(out,indent=2)+'\n';Path(__file__).with_suffix('.json').write_text(raw);print(raw)
