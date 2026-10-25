"""Recheck saved exact endpoints only. This script never calls the kernel."""
from pathlib import Path
from fractions import Fraction
import json
from flint import arb,ctx,fmpq
ctx.prec=1024
HERE=Path(__file__).resolve().parent

def endpoints(row):
    e=int(row['common_binary_exponent'])
    values=[Fraction(int(row[k]))*Fraction(2)**e for k in ['lower_mantissa','upper_mantissa']]
    return [arb(fmpq(v.numerator,v.denominator))for v in values]

av=json.loads((HERE/'anchor_validation.json').read_text())
native=json.loads((HERE/'anchor_native.json').read_text())
for n,r in zip(native['anchors'],av['anchors']):
    lo,hi=endpoints(n);rl,ru=endpoints(r['python640_exact_outward_dyadic_interval'])
    assert lo<=rl and hi>=ru
old=arb(json.loads((HERE/'validation_results.json').read_text())['old_certificate_display_read_verbatim'])
for fn in ['full_native.json','full_native_v2.json','full_native_v3.json']:
    r=json.loads((HERE/fn).read_text());assert r['full_group_terms']==967680
    lo,hi=endpoints(r['dyadic_interval']);assert lo>old.lower() and hi<old.upper()
assert json.loads((HERE/'full_native_v2.json').read_text())['dyadic_interval']==json.loads((HERE/'full_native_v3.json').read_text())['dyadic_interval']
print('PASS:6anchor intervals and3existing-point intervals;0kernel calls')
