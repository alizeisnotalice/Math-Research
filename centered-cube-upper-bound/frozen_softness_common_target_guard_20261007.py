from fractions import Fraction as F
from pathlib import Path
import json, hashlib
p=Path(__file__).resolve().parent
# Invert log(1+v)/v as a formal rational series, independently of the displayed B expansion.
a=[F(1),F(-1,2),F(1,3)]
b=[F(1)]
for k in (1,2): b.append(-sum(a[j]*b[k-j] for j in range(1,k+1)))
assert b[1:]==[F(1,2),F(-1,12)]
checks=1
rows=[]
for n in (2,16,512):
    corrections=[]
    for c in (F(1,2),F(3,4),F(1)):
        psi1=b[1]; psi2=b[2]-c*b[1]**2
        m2=2*c*psi1; m4=-24*c*psi2
        assert m2==c and m4==2*c+6*c*c
        checks+=1
        # L x_i^4 = -6 x_i^2 - (2+6c); L |x|^2=-n.
        constant=-F(n)*(1-F(3,n+2))*m4/c
        correction=-F(n-1,n+2)*(2+6*c)
        assert constant-correction*n==0
        assert -6+F(3,n+2)*2*(n+2)==0
        checks+=2
        corrections.append(correction)
        rows.append(dict(n=n,c=str(c),psi2=str(psi2),jump_m4=str(m4),harmonic_correction=str(correction)))
    for i in range(3):
        for j in range(i+1,3):
            assert corrections[i]!=corrections[j]
            checks+=1
out=dict(status='PASS_FINITE_MOMENT_COMMON_TARGET_RIGIDITY_COMPONENTS_ONLY',checks=checks,rows=rows,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(p/'frozen_softness_common_target_guard_20261007_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status=out['status'],checks=checks)))
