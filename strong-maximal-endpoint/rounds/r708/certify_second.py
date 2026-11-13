"""Directed enclosure for the exact general second-scale constant."""
import json,re
from fractions import Fraction
from pathlib import Path
from mpmath import iv
from physical import outward
HERE=Path(__file__).resolve().parent

def main():
    iv.dps=125
    w=json.loads((HERE/'second_constant.json').read_text())['values'][-1]
    mu_candidate=Fraction(re.sub(r'`[^*]*','',w['mu']).replace('*^','e'))
    scale=10**90;lower=(mu_candidate*scale).numerator//(mu_candidate*scale).denominator
    left=iv.mpf(lower)/scale;right=iv.mpf(lower+1)/scale
    F=lambda x:x-1-iv.exp(-x)
    assert F(left).b<0 and F(right).a>0
    mu=iv.mpf([left.a,right.b]);sigma=iv.sqrt(mu)
    out=dict(interval_precision=125,mu_interval=outward(mu,92),sigma_interval=outward(sigma,92),
        left_defining_function_interval=outward(F(left),120),right_defining_function_interval=outward(F(right),120),
        uniqueness='F(x)=x-1-exp(-x) has derivative1+exp(-x)>0; verified bracket identifies the unique positive root.',
        scope='sigma is the coefficient of sqrt(n log n) in mu*n-log C_n; no claim about the finer log n coefficient.')
    (HERE/'constant_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
