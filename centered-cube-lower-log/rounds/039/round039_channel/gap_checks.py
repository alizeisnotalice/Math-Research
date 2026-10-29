"""Exact arithmetic interface between continuous and bin-query signal."""
from pathlib import Path
from fractions import Fraction as F
import json,math
HERE=Path(__file__).resolve().parent
B=F(1,2)*F(1,48960)+F(1,510)*2*F(1,6)+F(1,510)*2*F(110,9)*F(1,48960)
eta0=F(127,255);etaerr=F(1,1020)*F(245,81)
assert F(49,100)<eta0-etaerr<eta0+etaerr<F(51,100)
correction=B*B/F(2499,10000)
assert B==F(29879,44945280)
assert correction==F(22318866025,12620438518754304)
trace=json.loads((HERE/'trace_results.json').read_text())
pack=[]
for row in trace['binary_float']:
    n=row['n'];D=math.isqrt(36*n)+1;lower=F(3,4)*(n//D+1)
    assert row['E_H']>=float(lower)
    pack.append({'n':n,'D':D,'packing_lower':str(lower),'observed_E_H':row['E_H']})
out={'B_exact':str(B),'B_float':float(B),'eta_interval_exact':[str(eta0-etaerr),str(eta0+etaerr)],
     'continuous_minus_bin_signal_upper_exact':str(correction),'upper_float':float(correction),
     'qualification':'This exact gap budget does not certify saved ordinary-FP bin signal values.',
     'binary_packing_checks':pack}
p=HERE/'gap_results.json'
if p.exists():raise RuntimeError('Refusing replacement')
p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
