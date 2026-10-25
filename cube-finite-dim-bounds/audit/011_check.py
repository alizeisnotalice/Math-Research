#!/usr/bin/env python3
"""Exact continuous cube-minimum identities for the harmonic quartic kernel."""
from pathlib import Path
import hashlib,json
import sympy as s
rows=[]
for n in (3,4):
    xs=s.symbols('x0:'+str(n),real=True);ds=s.symbols('d0:'+str(n),nonnegative=True)
    rr=sum(x*x for x in xs);h=sum(x**4 for x in xs)-s.Rational(3,n+2)*rr**2
    assert s.expand(sum(s.diff(h,x,2) for x in xs))==0
    D=sum(ds);hcorner=s.Rational(n)-s.Rational(3,n+2)*n*n
    delta=s.expand(sum((1-d)**2 for d in ds)-s.Rational(3,n+2)*(n-D)**2-hcorner)
    bound=s.Rational(n-1,n+2)*(4*D+D**2)
    slack=s.expand(bound-delta)
    assert s.expand(slack-2*sum(ds[i]*ds[j] for i in range(n) for j in range(i+1,n)))==0
    q=s.Rational(2*(n-1),n*(n+2));bcrit=(n-2)*(n+2)/(8*(n-1)*s.Integer(n)**s.Rational(n,2))
    C=(n-2)/(8*s.Integer(n)**s.Rational(n,2))
    assert s.simplify(bcrit*s.Rational(n-1,n+2)-C)==0
    rows.append({'n':n,'H4_corner':str(hcorner),'sphere_H4_min_coefficient':str(-q),'sphere_H4_max_coefficient':str(s.Rational(n-1,n+2)),'Bcrit':str(s.simplify(bcrit)),'quartic_deficit_upper_slack':str(slack),'continuous_cube_min':'n^((2-n)/2) + 2*n*(n-1)*B/(n+2)'})
t=s.symbols('t',positive=True);D3=3*(1-t*t)
rem3=s.factor(1/t-1-(4*D3+D3**2)/24)
assert s.simplify(rem3-(1-t)**3*(3*t*t+9*t+8)/(8*t))==0
d=s.symbols('D',nonnegative=True)
rem4=s.factor(1/(4-d)-s.Rational(1,4)-(4*d+d*d)/64)
assert s.simplify(rem4-d**3/(64*(4-d)))==0
B,R,n,q,A=s.symbols('B R n q A',positive=True)
# Stationary envelope substitution B*q*R^(n+2)=(n-2)/4.
env=R**(2-n)-A+B*q*R**4
stationary=s.simplify(env.subs(B,(n-2)/(4*q*R**(n+2)))-((n+2)*R**(2-n)/4-A))
assert stationary==0
out={'status':'exact continuous cube-minimum identities PASS','engine':'SymPy '+s.__version__,'scope':'n3/n4; B>=0 and B<=Bcrit, i.e. RB>=sqrt(n); no mass/cost computation','dimensions':rows,'n3_radial_remainder_times_sqrt3':str(rem3),'n3_remainder_domain':'t=sqrt((3-D)/3), 0<t<=1','n4_radial_remainder':str(rem4),'n4_remainder_domain':'0<=D<4','outer_stationary_envelope':'(n+2)/4*R^(2-n)-A','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('011_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
