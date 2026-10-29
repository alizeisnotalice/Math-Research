"""Independent exact counterexample and tilted-variance algebra checks."""
from pathlib import Path
from fractions import Fraction as F
import json,math,hashlib

HERE=Path(__file__).resolve().parent
OTHER=HERE.parent/'round039_posterior'
checks=0
for p in [F(3,10),F(1,2),F(7,10)]:
    for eta in [F(49,100),F(1,2),F(51,100)]:
        a=p*(1-eta)/((1-p)*eta)
        assert F(49,119)<=a<=F(119,49)
        for r in [F(j,100)for j in range(101)]:
            ratio=a/(1+(a-1)*r)**2
            assert ratio>=min(a,1/a)>=F(49,119)>F(2,5)
            u=(r-eta)/(eta*(1-eta))
            assert abs(u*(p-eta))<=F(3,7)
            checks+=1
n=256
pb=sum((F(math.comb(n,k)*9**k*7**(n-k),16**n)for k in range(147,n+1)),F(0))
prior=sum((F(math.comb(n,k),2**n)for k in range(147,n+1)),F(0))
gap=pb-16*prior-F(1,8)
assert gap>F(87,1000)
out={'scope':'Exact Fraction parameter algebra and independent counterexample; no new MC',
     'rational_parameter_cases':checks,'counterexample_gap':str(gap),'gap_float':float(gap),
     'audited_proof_SHA256':hashlib.sha256((OTHER/'PROOF.md').read_bytes()).hexdigest(),
     'proof_reading':'Exact tilt, Fourier cap, Binomial integer mode, Chebyshev lower cap, endpoint limits and all-set summation reviewed; no substantive flaw found.',
     'limits':'No expectation of exp(C*d^2/tau) or maximum-scale bound follows.'}
target=HERE/'posterior_audit.json'
if target.exists():raise RuntimeError('Refusing replacement')
target.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
