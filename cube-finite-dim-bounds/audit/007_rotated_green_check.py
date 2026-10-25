#!/usr/bin/env python3
"""Exact domain-comparison counterexample to coordinatewise Green monotonicity."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import sympy as s
r=s.sqrt(2);a=s.Rational(19,20)*r
upper=s.simplify(1/a-1/(3-a))
lower=1/r-s.Rational(6,11)
delta=s.simplify(lower-upper)
c0=F(600,1439)-F(6,11);c1=F(190,1439)-F(1,38)
assert c1>0 and c0+c1*F(7,5)>0
assert F(7,5)**2<2<F(3,2)**2
assert s.simplify(delta-(s.Rational(c0.numerator,c0.denominator)+s.Rational(c1.numerator,c1.denominator)*r))==0
out={'scope':'strict analytic domain-comparison bounds, not direct heat quadrature','dimension':3,'domain':'|u|<3/2, |v|<4, |z|<4; u=(x+y)/sqrt2, v=(x-y)/sqrt2','points':{'x':['19/20','19/20','0'],'y':['1','-1','0']},'coordinatewise_abs_x_le_abs_y':True,'upper_4pi_times_Gx':str(upper),'lower_4pi_times_Gy':str(lower),'difference_lower_exact':str(delta),'rational_lower_for_4pi_difference_using_sqrt2_gt_7p5':str(c0+c1*F(7,5)),'numerical_upper_Gx':float(upper/(4*s.pi)),'numerical_lower_Gy':float(lower/(4*s.pi)),'status':'G_D(y)>G_D(x) strictly proved; does not disprove minimum-at-some-corner','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('007_rotated_green_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
