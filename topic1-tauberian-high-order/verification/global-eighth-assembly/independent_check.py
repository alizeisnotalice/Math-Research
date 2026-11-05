"""Independent exact arithmetic for R057 review; no geometric claims."""
from fractions import Fraction as F
from math import comb, factorial
import json


def series(n, length):
    def scalar(c):
        return [F(c)] + [F(0)] * length

    def add(*values):
        return [sum(value[j] for value in values) for j in range(length + 1)]

    def scale(value, c):
        return [c * x for x in value]

    def shift(value, k):
        return [F(0)] * k + value[:length + 1 - k]

    def mul(p, q):
        return [sum(p[k] * q[j-k] for k in range(j+1)) for j in range(length+1)]

    def reciprocal(value):
        out = scalar(1 / value[0])
        for j in range(1, length + 1):
            out[j] = -sum(value[k] * out[j-k] for k in range(1, j+1)) / value[0]
        return out

    def binomial(exponent):
        out = scalar(1)
        for j in range(1, length + 1):
            out[j] = out[j-1] * (exponent-j+1) / j
        return out

    A = binomial(n)
    a = shift(scalar(1), n)
    K = F((n-1)*n**3, 54)
    c8 = -F(13,18)-F(n*(n-2),3)+F(400,243)
    Q = scale(add(shift(binomial(-4),6),
                  scale(shift(binomial(-5),7),F(10,9)),
                  scale(shift(binomial(-2),8),c8)),K)
    if n == 2:
        Q = add(Q, shift(scalar(F(8,27)),8))
    actual = mul(A, reciprocal(add(A,scale(a,-1),scale(mul(a,Q),-1))))
    Jhat = F((n-1)*n**3*(2609+1971*n+81*n*n),26244)
    if n == 2:
        Jhat += F(16,27)
    expected = []
    for j in range(9):
        baseline = sum((-1)**(j-(q-1)*n)*comb(n+j-1,j-(q-1)*n)
                       for q in range(1,2+j//n))
        correction = K if j == 6 else (-(n+F(26,9))*K if j == 7 else (Jhat if j == 8 else F(0)))
        expected.append(F(baseline)+correction)
    residuals = [actual[n+j]-expected[j] for j in range(9)]
    assert residuals == [0]*9
    return {'n':n,'residuals':list(map(str,residuals)),
            'Jhat':str(Jhat),'baseline_j8':str(expected[8]-Jhat),
            'ordinary_beta8':str(expected[8]*2**(n+8))}


# The generic bracket has coefficients in n of degrees 0,1,2.
generic_bracket = [F(40,9)-F(13,18)+F(400,243),
                   F(2,3)+F(26,9)+F(1,2),-F(1,3)+F(1,2)]
assert generic_bracket == [F(2609,486),F(1971,486),F(81,486)]

rows = [series(n,n+8) for n in range(2,13)]
sigma5 = []
for q in (F(0),F(1,100),F(1,2),F(1)):
    product = F(1)
    for r in range(1,5):
        product *= r-q
    sigma5.append({'q':str(q),'coefficient':str(product/factorial(5))})

result = {'arithmetic':'Python Fraction exact',
          'generic_bracket_coefficients':list(map(str,generic_bracket)),
          'dimensions_2_through_12':rows,
          'sigma5':sigma5,
          'uniform_sigma6_bound':str(F(factorial(5)*64,factorial(6))),
          'geometry_certified_by_script':False,
          'all_assertions_passed':True}
print(json.dumps(result,ensure_ascii=False,indent=2))
