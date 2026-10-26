#!/usr/bin/env python3
"""Exact symbolic derivative of one PF term; no actual H/log evaluation."""
import sys
sys.dont_write_bytecode=True
import sympy as s
import common as C

def main():
    t=s.symbols('t');A,B,det=s.symbols('A B det',nonzero=True);k=s.symbols('k:5');v=s.symbols('v:5')
    nodes=[k[j]+t*v[j] for j in range(5)];D=s.prod(nodes[0]-nodes[j] for j in range(1,5));w=A*A+B*B*nodes[0]**2
    term=-det*s.log(w)/(2*D);D0=D.subs(t,0);w0=w.subs(t,0)
    T=sum((v[0]-v[j])*s.prod(k[0]-k[m] for m in range(1,5) if m!=j) for j in range(1,5))
    formula=det*T*s.log(w0)/(2*D0**2)-det*B*B*v[0]*k[0]/(D0*w0)
    residual=s.factor(s.diff(term,t).subs(t,0)-formula);assert residual==0
    C.write(C.HERE/'symbolic_probe.json',dict(status='exact_onepole_fourphysicalderivative_formula_verified',sympy_version=s.__version__,residual=str(residual),new_H=0,new_gradient_values=0,new_numerical_log=0,route='symbolic identity only'))
    print('exact derivative residual: 0; no H or numerical logarithm evaluated')

if __name__=='__main__':main()
