"""Exact rational certificate interface for the rank>=4 remainder only.

No density/kernel evaluation; no frame is executed on import. The accompanying
checker calls this only for its one preregistered 6x8 toy frame, never body32.
All returned matrix/scalar majorants carry an explicit common factor pi^-6.
"""
from fractions import Fraction as Q
from itertools import combinations


def matrix(rows):
    return [[Q(x) for x in row] for row in rows]


def rank(rows):
    a = matrix(rows)
    if not a or not a[0]:
        return 0
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][c]
        a[r] = [x / scale for x in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][c]
            a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def selected(a, cols):
    return [[row[j] for j in cols] for row in a]


def determinant_and_inverse(rows):
    a = matrix(rows)
    n = len(a)
    assert n and all(len(row) == n for row in a)
    inverse = [[Q(int(i == j)) for j in range(n)] for i in range(n)]
    determinant = Q(1)
    for c in range(n):
        pivot = next((i for i in range(c, n) if a[i][c]), None)
        if pivot is None:
            return Q(0), None
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            inverse[c], inverse[pivot] = inverse[pivot], inverse[c]
            determinant = -determinant
        scale = a[c][c]
        determinant *= scale
        a[c] = [x / scale for x in a[c]]
        inverse[c] = [x / scale for x in inverse[c]]
        for i in range(n):
            if i != c:
                scale = a[i][c]
                a[i] = [x - scale * y for x, y in zip(a[i], a[c])]
                inverse[i] = [x - scale * y for x, y in zip(inverse[i], inverse[c])]
    return determinant, inverse


def independent_four_bases(rows):
    """Lexicographic J4 and first lexicographic two-column completion K(J).

    N is ordered as J then K. Its inverse is fixed before any differentiation,
    and is shared by every pattern assigned to J. No pattern is enumerated.
    """
    a = matrix(rows)
    assert len(a) == 6 and len({len(row) for row in a}) == 1
    m = len(a[0])
    assert m >= 6 and rank(a) == 6
    records = []
    for j in combinations(range(m), 4):
        if rank(selected(a, j)) != 4:
            continue
        remainder = [i for i in range(m) if i not in j]
        for k in combinations(remainder, 2):
            n = selected(a, j + k)
            det, inverse = determinant_and_inverse(n)
            if det:
                records.append(dict(J=j, K=k, N=n, determinant=det, inverse=inverse))
                break
        else:
            raise AssertionError('Independent four columns failed to extend in a spanning frame')
    return records


def scalar_constants(radius, cutoff_time):
    r, t = Q(radius), Q(cutoff_time)
    assert r > 0 and 0 < t <= 1
    a1, a2 = 4 / r, 20 / r**2
    b1, b2 = 1 + 4 / r, 2 + 8 / r + 12 / r**2
    dtail = (a2 + 3 * a1**2) * t**3 / 3 + a1 * b1 * t**2
    dcompletion = (b2 + b1**2) * t + 2 * a1 * b1 * t**2
    return dict(R=r, T=t, a1=a1, a2=a2, b1=b1, b2=b2,
                d_tail=dtail, d_completion=dcompletion,
                h0=Q(1), h1=max(a1, b1), h2=max(a2, b2, max(a1, b1)**2))


def majorant_from_bases(records, radius, cutoff_time):
    """Return pi^6 M4, pi^6 tr(M4), and pi^6 K_k for k=0,1,2.

    This helper assumes records comprise every independent J exactly once;
    certificate() supplies that exhaustive list without enumerating patterns.
    pi^6 K_k bounds all unit-direction k-linear derivatives by t^(2-k).
    """
    c = scalar_constants(radius, cutoff_time)
    out = [[Q(0) for _ in range(6)] for _ in range(6)]
    ks = [Q(0), Q(0), Q(0)]
    for record in records:
        inverse = record['inverse']
        scale = 1 / (c['R']**8 * abs(record['determinant']))
        weights = [c['d_tail']] * 4 + [c['d_completion']] * 2
        for i in range(6):
            for j in range(6):
                out[i][j] += scale * sum(weights[k] * inverse[k][i] * inverse[k][j]
                                         for k in range(6))
        # A deliberately rational upper bound on ||N^-1||_(2 -> 1).
        lam = sum(abs(x) for row in inverse for x in row)
        for k in range(3):
            ks[k] += scale * c['h' + str(k)] * lam**k
    return dict(pi6_M4=out, pi6_trace_M4=sum(out[i][i] for i in range(6)),
                pi6_K012=ks, constants=c)


def certificate(rows, radius, cutoff_time):
    bases = independent_four_bases(rows)
    return dict(scope='rank>=4 cutoff remainder only; no global Green/Cn claim',
                independent_four_count=len(bases), basis_records=bases,
                **majorant_from_bases(bases, radius, cutoff_time))


if __name__ == '__main__':
    raise SystemExit('Library interface only. Use a separately registered finite caller.')
