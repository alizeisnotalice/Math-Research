"""Read saved full dyadics and B011 original ball. Zero kernel calls."""
from pathlib import Path
from fractions import Fraction as Q
import json
from flint import arb, fmpq, ctx
ctx.prec=640
here=Path(__file__).resolve().parent
x=json.loads((here/'native_result.json').read_text())
ref=json.loads((here.parent/'cutoff48/results.json').read_text())
row=next(v for v in ref['rows']if v['label']=='A007_reference')
assert x['reference_complete_arb384']==row['H_arb384']
r=arb(row['H_arb384']);n=x['native'];d=n['dyadic_interval']
factor=Q(2)**int(d['common_binary_exponent'])
lo=Q(int(d['lower_mantissa']))*factor;hi=Q(int(d['upper_mantissa']))*factor
b=lambda v:arb(fmpq(v.numerator,v.denominator))
assert r.lower()<b(lo)<=b(hi)<r.upper()
assert x['exit_code']==0 and n['full_group_terms']==967680
assert n['precision_bits']==512 and n['computed_DD_evaluations']==483840
print('PASS primary saved H48 cross-engine; zero kernel calls')
