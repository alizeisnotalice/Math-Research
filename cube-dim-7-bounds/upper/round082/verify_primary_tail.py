from fractions import Fraction as F
from math import comb
from pathlib import Path
import json,hashlib
p=Path(__file__).resolve().parent
coeff=[F(-5,2),F(-3,2),F(-1,4),F(1,12),F(1,16)]
shift=[sum((coeff[j]*comb(j,k)*4**(j-k) for j in range(k,5)),F(0)) for k in range(5)]
assert shift==[F(53,6),F(33,2),F(27,4),F(13,12),F(1,16)] and all(v>0 for v in shift)
out={'status':'PASS','shifted_coefficients':[str(v) for v in shift],'registration_sha256':hashlib.sha256((p/'primary_tail_registration.json').read_bytes()).hexdigest(),'scope':'polynomial identity only; supports y*S16prime-(5/2)*S16>0 for y>=4 after retaining quartic and positive higher coefficients'}
(p/'primary_tail_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
