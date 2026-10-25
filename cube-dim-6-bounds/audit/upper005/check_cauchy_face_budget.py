#!/usr/bin/env python3
"""Exact algebra audit for the round005 Cauchy-product face estimate."""
from __future__ import annotations
import hashlib,json
from itertools import product
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
N=sp.Integer(32); t,s=sp.symbols('t s', positive=True, real=True)
z=sp.symbols('z0:32', real=True)
sumz=sum(z); sumsq=sum(q*q for q in z)
product_gap=sp.expand(N*(N*t*t+sumsq)-(sumz*sumz+N*N*t*t))
pair_gap=sp.expand(sum((z[i]-z[j])**2 for i in range(32) for j in range(i+1,32)))
Cf=sp.factor(N*N/sp.Integer(12)**5)
K=sp.factor(sp.Rational(6,7)*Cf**sp.Rational(1,6)*sp.pi**sp.Rational(-sp.Rational(11,6)))
signs=list(product((-1,1),repeat=5))
one=(1,)*5
quotient=[tuple(a-b for a,b in zip(v,one,strict=True)) for v in signs if v!=one]
rank=sp.Matrix(quotient).rank()
equal_sum=sum(abs(sum(v)) for v in signs)
out={
 'status':'PASS_exact_cauchy_product_face_algebra',
 'product_inequality_gap_equals_pair_squares':sp.simplify(product_gap-pair_gap)==0,
 'quotient_rank':rank,
 'quotient_equal_coordinate_abs_sum':equal_sum,
 'Cf_without_pi6':str(Cf),
 'Knew':str(K),
 'l1_interpolation_exponents':'sup: (1/(243*pi^6))*t^-4; positive-face L1: t/pi; q=6/5 gives pi^(-11/6)*243^(-1/6)*t^(1/6)',
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
assert out['product_inequality_gap_equals_pair_squares'] and rank==5 and equal_sum==60 and Cf==sp.Rational(1,243)
(HERE/'check_cauchy_face_budget_results.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
