from pathlib import Path
from fractions import Fraction as Q
from flint import fmpq
import json,hashlib
from exact_reference import sign_dd
HERE=Path(__file__).resolve().parent
reg=json.loads((HERE/'registration.json').read_text());native=json.loads((HERE/'anchors_native.json').read_text())

def bounds(row):
    factor=Q(2)**int(row['common_binary_exponent'])
    return Q(int(row['lower_mantissa']))*factor,Q(int(row['upper_mantissa']))*factor

def enclosed(v,row):
    v=Q(int(v.numerator),int(v.denominator));lo,hi=bounds(row);assert lo<=v<=hi
    return {'numerator':str(v.numerator),'denominator':str(v.denominator)}

checks=[]
for i,(nodes,row)in enumerate(zip(reg['anchors'],native['anchors'])):
    z=[fmpq(n,d)for n,d in nodes]
    if any(not v for v in z):
        assert row=={'unsupported_zero_node':True};checks.append({'index':i,'zero_rejected':True});continue
    value=sign_dd(z);out={'index':i,'value_exact':enclosed(value,row['value']),'partials_exact':[]}
    for j in range(6):out['partials_exact'].append(enclosed(sign_dd(z+[z[j]]),row['partials'][j]))
    checks.append(out)
result={'status':'PASS_ALL6_CASES_FIVE_NONZERO_PLUS_ONE_ZERO_REJECTION','cases':checks,
        'native_source_sha256':hashlib.sha256((HERE/'generic_green.cpp').read_bytes()).hexdigest(),
        'reference_source_sha256':hashlib.sha256((HERE/'exact_reference.py').read_bytes()).hexdigest()}
(HERE/'anchor_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['status'])
