#!/usr/bin/env python3
"""Small new symbolic identities only; no scalar-budget or H evaluations."""
import sys
sys.dont_write_bytecode=True
import sympy as s
import json,os,hashlib
from pathlib import Path

def main():
    H=Path(__file__).resolve().parent;t,e,C3,C1,P5=s.symbols('t e C3 C1 P5',positive=True)
    B=C3*(e**-2-t**-2)/2+C1*s.log(t/e)+P5/(4*t**4)
    derivative=s.factor(s.diff(B,t)-(C1*t**4+C3*t**2-P5)/t**5);assert derivative==0
    primitive=s.log(1+256*t*t)/(512*s.pi);mass=s.factor(s.diff(primitive,t)-t/(s.pi*(1+256*t*t)));assert mass==0
    h,r,u=s.symbols('h r u',positive=True);q=s.Rational(5,4)
    remainder=s.simplify(q*(h**(q-1)*r+r**q/q)-(q*h**(q-1)*r+r**q));assert remainder==0
    P3=1/(1024*s.sqrt(2)*s.pi**2);radial=(2*s.pi)**-5*(8*s.pi**2/3)*s.factorial(4)/(8*s.sqrt(2))**5
    peak=s.simplify(radial-P3/(64*s.pi));assert peak==0
    out=dict(status='four_new_exact_symbolic_identities_passed',pid=os.getpid(),sympy=s.__version__,residuals=dict(B_derivative=str(derivative),trace_mass_primitive=str(mass),subadditive_integral_coefficients=str(remainder),P5_from057_general_formula=str(peak)),new_H=0,new_numeric_log=0,new_budget_cases=0,
        note='inequality and extended-value MCT require the written proof;symbolic identities alone are not that proof')
    (H/'symbolic_result.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
