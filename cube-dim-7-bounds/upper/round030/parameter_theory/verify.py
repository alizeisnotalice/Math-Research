"""Finite algebra checks and high precision root candidates; no completeness claim."""
from pathlib import Path
import hashlib
import json
import mpmath as mp
import sympy as sp

mp.mp.dps = 90
ROOT = Path(__file__).resolve().parent
u, h = sp.symbols('u h')
checks = []
for n in (2, 4, 7, 64):
    q = sp.Poly((1-u/sp.Integer(n))**(n-1), u)
    def moment(k):
        p = sp.Poly(q.as_expr()*u**k, u).integrate()
        return sp.factor(p.eval(n)-p.eval(0))
    assert moment(0) == 1
    assert moment(1) == sp.Rational(n, n+1)
    assert moment(2) == sp.Rational(2*n*n, (n+1)*(n+2))
    v = sp.symbols('v')
    mass = sp.integrate(n*v**(n-1), (v, h, 1))
    first = sp.integrate(n*n*(1-v)*v**(n-1), (v, h, 1))
    assert sp.expand(mass-(1-h**n)) == 0
    assert sp.expand(first-(sp.Rational(n,n+1)-n*h**n+sp.Rational(n*n,n+1)*h**(n+1))) == 0
    checks.append(f'exact moments and truncated moments n={n}')
assert 16*sp.Rational(3,4)**16 < 1
assert sp.Rational(1,4)*(1+sp.Rational(1,128))*sp.Rational(33,32)**2 < 1
L = mp.mpf(18)
assert 3*mp.log(L)/L < mp.mpf('0.5')
assert L*L/120 > 1+2/L
assert 3*mp.log(L) > 2
checks.append('displayed-bound thresholds at L=18; analytic monotonicity in proof')

def parts(n, t):
    x = n/t
    H = (n-x)/2
    hh = (1+1/t)/2
    I = mp.quad(lambda z: (1-z/n)**(n-1)/(z+x/2), [0,H])
    J = mp.quad(lambda z: (1-z/n)**(n-2)/(z+x/2), [0,H])
    A = (1+1/t)**n
    B = 1-hh**n+n*mp.log(t)-n*I/2
    D = n/(2*t)+n*hh**(n-1)/(2*t*t)+n*(n-1)*J/(4*t*t)
    g = D/mp.sqrt(B)-n/t**2*(1+1/t)**(n/2-1)
    return g, (mp.sqrt(A)+mp.sqrt(B))**2

def bisect(f, a, b):
    fa, fb = f(a), f(b)
    assert fa*fb < 0
    for _ in range(200):
        c = (a+b)/2
        fc = f(c)
        if fa*fc < 0:
            b = c
        else:
            a, fa = c, fc
    return (a+b)/2

candidates = []
for nn in (4,7,64):
    n = mp.mpf(nn)
    scale = n*2**nn
    z = bisect(lambda z: parts(n,1+z/scale)[0], mp.mpf('0.25'),mp.mpf(64))
    t = bisect(lambda t: parts(n,t)[0],mp.mpf(2),4*n)
    candidates.append({'n':nn,'z_near_endpoint':mp.nstr(z,45),
                       't_bulk':mp.nstr(t,45),'F_bulk':mp.nstr(parts(n,t)[1],45),
                       'status':'high precision candidate only; no root-completeness certificate'})
receipt = {'mpmath':mp.__version__,'sympy':sp.__version__,'decimal_precision':mp.mp.dps,
           'exact_checks':checks,'candidates':candidates,
           'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                            for name in ('proof.md','verify.py')}}
(ROOT/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
