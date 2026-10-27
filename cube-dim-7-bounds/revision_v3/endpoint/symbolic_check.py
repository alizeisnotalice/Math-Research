"""Exact checks for endpoint/proof.md. Requires SymPy; no numerical search."""
import sympy as s

t, x = s.symbols('t x', positive=True)
L = s.log(t)


def B(n):
    c = 1 + 1/(2*t)
    integral = c**(n-1)*L + sum(
        s.binomial(n-1, j)*c**(n-1-j)*(-s.Rational(1, 2))**j
        *(1-t**(-j))/j for j in range(1, n)
    )
    return s.expand(1-((1+1/t)/2)**n+n*L-s.Rational(n, 2)*integral)


B2 = s.Rational(5, 4)-x-x*x/4+(1-x/2)*(-s.log(x))
B3 = (35-18*x-15*x*x-2*x**3+(24-24*x-6*x*x)*(-s.log(x)))/16
assert s.simplify(B(2).subs(t, 1/x)-B2) == 0
assert s.simplify(B(3).subs(t, 1/x)-B3) == 0
d = 1-x
assert s.simplify(B2-d*d-(d*(5*x+1)/4+(1-x/2)*(-s.log(x)))) == 0
P = -37+126*x-87*x*x-2*x**3
Q = 24-24*x-6*x*x
assert s.simplify(16*(B3-s.Rational(9, 2)*d*d)-P-Q*(-s.log(x))) == 0
assert s.simplify(P+Q*2*d/(1+x)-d*(2*x**3+79*x*x+4*x+11)/(1+x)) == 0
assert s.simplify(P+Q*d/x-d*(2*x**3+83*x*x-61*x+24)/x) == 0
assert 61**2-4*83*24 == -4247

b42 = s.simplify(B(4).subs(t, 2))
assert s.simplify(b42-(1841+72*s.log(2))/768) == 0
assert s.Rational(1913, 768) < s.Rational(49, 16)
delta5 = 2**s.Rational(5, 2)-s.Rational(3, 2)**s.Rational(5, 2)
assert s.simplify(delta5**2-(s.Rational(1267, 32)-18*s.sqrt(3))) == 0
assert s.Rational(7, 4)**2 > 3
assert s.Rational(1267, 32)-18*s.Rational(7, 4) == s.Rational(259, 32)
assert s.Rational(259, 32) > 6

for n in range(1, 9):
    assert s.simplify(B(n).subs(t, 1)) == 0
    assert s.simplify(s.diff(B(n), t).subs(t, 1)) == n

print('PASS: exact B2/B3/B4 formulas, logarithm-bound factorizations,')
print('      discriminant, n>=4 witness certificates, B(1), and Bprime(1).')
print('No grid search or floating-point certification was used.')
