"""Independent finite algebra checks; no imported author code, no geometry claims."""
from fractions import Fraction as F
from math import comb, factorial
from hashlib import sha256
from pathlib import Path
import json


def mul(a, b, degree):
    out = [F(0)] * (degree + 1)
    for i, x in enumerate(a[:degree + 1]):
        for j, y in enumerate(b[:degree + 1 - i]):
            out[i + j] += x * y
    return out


def logarithm(a, degree):
    t = a[:degree + 1] + [F(0)] * max(0, degree + 1 - len(a))
    assert t[0] == 1
    t[0] = F(0)
    power = [F(1)] + [F(0)] * degree
    out = [F(0)] * (degree + 1)
    for h in range(1, degree + 1):
        power = mul(power, t, degree)
        weight = F((-1) ** (h + 1), h)
        for k in range(degree + 1):
            out[k] += weight * power[k]
    return out


def exponential(t, degree):
    assert t[0] == 0
    out = [F(1)] + [F(0)] * degree
    for k in range(1, degree + 1):
        out[k] = sum(F(i) * t[i] * out[k - i] for i in range(1, k + 1)) / k
    return out


def direct_power(a, p, degree):
    return exponential([p * x for x in logarithm(a, degree)], degree)


def c_product(q, degree):
    out = [F(1)]
    for k in range(1, degree + 1):
        value = F(1, factorial(k + 1))
        for l in range(1, k + 1):
            value *= l - q
        out.append(value)
    return out


def a_recurrence(c, p, degree):
    out = [F(1)]
    for k in range(1, degree + 1):
        out.append(sum(((p + 1) * i - k) * c[i] * out[k - i]
                       for i in range(1, k + 1)) / k)
    return out


q_values = (F(0), F(1, 3), F(7, 11), F(1))
power_checks = 0
for q in q_values:
    c = c_product(q, 8)
    assert c[0] == 1
    for k in range(1, 9):
        assert c[k] == (k - q) * c[k - 1] / (k + 1)
        assert 0 <= c[k] <= F(1, k + 1)
    for p in (F(-3, 2), F(0), F(1), F(7, 3)):
        rec = a_recurrence(c, p, 8)
        direct = direct_power(c, p, 8)
        assert rec == direct
        power_checks += 9

composition_cases = []
coefficient_assertions = 0
for n in (2, 3, 5, 7):
    for r in (0, 3, 6, 10, 15):
        for q in q_values:
            degree = n + r
            c = c_product(q, r // n)
            # Build A_q(w^n), then exp(p log A_q(w^n)) as a w-series.
            aw = [F(0)] * (degree + 1)
            for k, value in enumerate(c):
                aw[n * k] = value
            for j in range(r + 1):
                p = F(n + j, n)
                direct = direct_power(aw, p, degree)
                rec = a_recurrence(c, p, r // n)
                for J in range(r + 1):
                    coefficient = direct[J - j] if J >= j else F(0)
                    index, residue = divmod(J - j, n)
                    predicted = rec[index] if J >= j and residue == 0 else F(0)
                    assert coefficient == predicted, (n, r, q, j, J)
                    if q == 1:
                        assert coefficient == int(J == j)
                    coefficient_assertions += 1
            composition_cases.append({'n': n, 'r': r, 'q': str(q), 'pass': True})

baseline_checks = 0
baseline_rows = []
for n in range(2, 10):
    degree = 24
    ell = [F(1)]
    for j in range(1, degree + 1):
        ell.append(-sum(comb(n, k) * ell[j-k]
                        for k in range(1, min(j, n - 1) + 1)))
    # Independently sum each geometric-series term v^(nq)(1+v)^(-nq).
    direct = [F(0)] * (degree + 1)
    for count in range(1, degree // n + 2):
        offset = (count - 1) * n
        for h in range(degree - offset + 1):
            direct[offset + h] += (-1) ** h * comb(n * count + h - 1, h)
    for j in range(degree + 1):
        closed = sum((-1) ** (j-(count-1)*n) * comb(n+j-1, j-(count-1)*n)
                     for count in range(1, j // n + 2))
        assert ell[j] == direct[j] == closed, (n, j)
        baseline_checks += 1
    baseline_rows.append({'n': n, 'through_j': degree, 'pass': True})
assert all(F((-2) ** j) == (-sum(comb(2, k) * F((-2) ** (j-k))
                                  for k in range(1, min(j, 1) + 1)))
           for j in range(1, 25))
assert sum((-1) ** (5-(k-1)*3) * comb(7, 5-(k-1)*3) for k in (1, 2)) == 0

result = {
    'runtime': 'Python standard library Fraction; self-written log/exp formal series',
    'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    'power_coefficient_checks': power_checks,
    'composition_cases': len(composition_cases),
    'composition_basis_assertions': coefficient_assertions,
    'composition_scope': {'n': [2, 3, 5, 7], 'r': [0, 3, 6, 10, 15],
                          'q': [str(q) for q in q_values]},
    'baseline_checks': baseline_checks,
    'baseline_scope': baseline_rows,
    'm1_identity': True,
    'q0_auxiliary_only': True,
    'r0_included': True,
    'n3_j5_baseline_zero': True,
    'all_pass': True,
    'scope_limit': 'Finite exact algebra only; not a proof of true C_n assumptions or uniform remainders.'
}
print(json.dumps(result, ensure_ascii=False, indent=2))
