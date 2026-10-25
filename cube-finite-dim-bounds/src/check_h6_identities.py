"""Exact identities only: H6 harmonicity and coordinate derivative difference.

No assertion of an extremum, radial-root completeness, or maximal bound.
"""
import json
import sympy as s

rows = []
for n in (3, 4):
    x = s.symbols('x:'+str(n))
    r2 = sum(z*z for z in x)
    a = s.Rational(15, n+8)
    b = s.Rational(30, (n+8)*(n+4))
    h6 = sum(z**6 for z in x)-a*r2*sum(z**4 for z in x)+b*r2**3
    lap = s.expand(sum(s.diff(h6, z, 2) for z in x))
    assert lap == 0
    y = s.symbols('y:'+str(n))
    total = sum(y)
    B4, B6 = s.symbols('B4 B6')
    h4y = sum(z*z for z in y)-s.Rational(3,n+2)*total**2
    h6y = sum(z**3 for z in y)-a*total*sum(z*z for z in y)+b*total**3
    f = -B4*h4y-B6*h6y
    difference = s.factor(s.diff(f,y[0])-s.diff(f,y[1]))
    expected = (y[0]-y[1])*(-2*B4-3*B6*(y[0]+y[1])+2*a*B6*total)
    assert s.expand(difference-expected) == 0
    rows.append({'n':n,'a':str(a),'b':str(b),'laplacian':str(lap),
                 'cube_derivative_difference':str(difference)})
print(json.dumps({'status':'exact identities PASS','sympy':s.__version__,
                  'rows':rows,'scope':'identities only; no new upper bound'},indent=2))
