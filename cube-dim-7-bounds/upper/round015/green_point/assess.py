#!/usr/bin/env python3
"""Exact threshold and fixed-budget tests from the saved outward G endpoint; no phases evaluated."""
import json
from fractions import Fraction
from pathlib import Path
from flint import arb,fmpq,ctx
ctx.prec=256;p=Path(__file__).resolve().parent;d=json.loads((p/'result.json').read_text());assert d['status']=='completed' and d['positive'] and d['width_pass']
def dy(v):return Fraction(int(v[0]))*Fraction(2)**int(v[1])
G=dy(d['full_G']['upper']);assert G>0
thresholds={str(k):dict(strict=G<Fraction(1,k),upper_gap=str(Fraction(1,k)-G)) for k in (219000,194000,179000)}
I0=Fraction(7723627498344,10**18);cost6=Fraction(7**7)*I0**6/G**7;passes=cost6>128**6
x=arb(fmpq(cost6.numerator,cost6.denominator));root=(x.log()/6).exp();out=dict(G_upper_exact=str(G),thresholds=thresholds,fixed_I0=str(I0),fixed_cost_power6_lower_exact=str(cost6),fixed_cost_power6_strictly_gt_128_power6=passes,fixed_cost_lower_display=root.str(35),scope='If strict, excludes <128 only from this fixed-I0 upper-bound interface using any legal g; not a C7 lower bound, and tighter I may change it.')
(p/'assessment.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('G_upper_exact','fixed_cost_power6_lower_exact')},indent=2))
