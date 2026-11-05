"""Small independent exact curve check, using (delta,t) Fraction coefficients."""
from fractions import Fraction as F
from math import comb
import json


def plus(*items):
    result = {}
    for item in items:
        for key, value in item.items():
            result[key] = result.get(key, F(0)) + value
    return {key: value for key, value in result.items() if value}


def scale(item, factor):
    return {key: value * factor for key, value in item.items() if value * factor}


def product(a, b, order=None):
    result = {}
    for (da, ta), va in a.items():
        for (db, tb), vb in b.items():
            if order is None or da + db <= order:
                key = (da + db, ta + tb)
                result[key] = result.get(key, F(0)) + va * vb
    return {key: value for key, value in result.items() if value}


def derivative(a):
    return {(d, t - 1): value * t for (d, t), value in a.items() if t}


def primitive(a):
    return {(d, t + 1): value / (t + 1) for (d, t), value in a.items()}


def endpoint(a, right):
    out = {}
    for (d, t), value in a.items():
        if right or t == 0:
            out[(d, 0)] = out.get((d, 0), F(0)) + value
    return {key: value for key, value in out.items() if value}


def reverse(a):
    out = {}
    for (d, t), value in a.items():
        for k in range(t + 1):
            out[(d, k)] = out.get((d, k), F(0)) + value * comb(t, k) * (-1) ** k
    return {key: value for key, value in out.items() if value}


one = {(0, 0): F(1)}
delta = {(1, 0): F(1)}
initial = [({(0, 1): F(1)}, {(1, 1): F(1)}, {(1, 2): F(1, 2)})]


def project(a, order):
    return {key: value for key, value in a.items() if key[0] <= order}


def advance(curves, order=None):
    total = endpoint(curves[-1][2], True)
    output = []
    for high in (False, True):
        old_pieces = reversed(curves) if high else curves
        for x, y, z in old_pieces:
            if high:
                x, y, z = reverse(x), reverse(y), reverse(z)
                q = plus(one, scale(x, -1), y)
                m = plus(scale(total, 2), scale(z, -1), scale(product(y, y, order), F(1, 2)))
            else:
                q = y
                m = plus(scale(product(y, y, order), F(1, 2)),
                         product(y, plus(one, scale(x, -1)), order), z)
            output.append((q, product(delta, plus(q, scale(m, -1)), order)))
    prefix = {}
    result = []
    for x, y in output:
        z = plus(prefix, primitive(product(y, derivative(x), order)))
        if result:
            assert endpoint(result[-1][0], True) == endpoint(x, False)
            assert endpoint(result[-1][1], True) == endpoint(y, False)
            assert endpoint(result[-1][2], True) == endpoint(z, False)
        result.append((x, y, z))
        prefix = endpoint(z, True)
        assert derivative(z) == product(y, derivative(x), order)
    assert endpoint(result[0][0], False) == {}
    assert endpoint(result[-1][0], True) == one
    return result


def total_area(curves):
    return scale(endpoint(curves[-1][2], True), 2)


def moment(curves, order):
    out = {}
    for x, y, z in curves:
        density = plus(product(plus(one, x), y, order), scale(product(y, y, order), -1))
        out = plus(out, scale(endpoint(primitive(product(density, derivative(x), order)), True), 2))
    return out


def evaluate(a, d, t):
    return sum(value * d ** i * t ** j for (i, j), value in a.items())


# Check the first generated graph against direct two-piece min-cap integration.
curves = advance(initial)
point_checks = 0
d = F(1, 13)
for x, y, z in curves:
    for t in (F(0), F(1, 3), F(2, 3), F(1)):
        q = evaluate(x, d, t)
        if q <= d:
            m0 = q - (1-d) * q*q / (2*d)
        else:
            m0 = d * (1-2*d+2*q-q*q) / (2*(1-d))
        assert evaluate(y, d, t) == d * (q-m0)
        point_checks += 1

# Compare exact curves with projection through j=2, including every primitive offset.
exact = initial
truncated = initial
for layer in range(3):
    assert len(exact) == len(truncated) == 2 ** layer
    assert [tuple(project(p, 6) for p in piece) for piece in exact] == truncated
    if layer < 2:
        exact = advance(exact)
        truncated = advance(truncated, 6)

# Small independent truncated run through R=6,j=5, with no piece pruning.
order = 6
curves = initial
rows = []
previous_moment = None
for layer in range(6):
    area = total_area(curves)
    cap_moment = moment(curves, order)
    if previous_moment is not None:
        assert area == plus(delta, scale(product(delta, previous_moment, order), -1))
    rows.append({'layer': layer, 'pieces': len(curves),
                 'area': [str(area.get((i, 0), F(0))) for i in range(order + 1)],
                 'formally_flat_pieces': sum(not derivative(x) for x, y, z in curves),
                 'flat_pieces_with_nonzero_primitive_prefix':
                     sum(not derivative(x) and bool(endpoint(z, False)) for x, y, z in curves)})
    previous_moment = cap_moment
    if layer < 5:
        curves = advance(curves, order)
expected = list(map(F, [0, 1, '-5/3', '13/4', '-33/5', '1657/120', '-1661/56']))
assert list(map(F, rows[-1]['area'])) == expected

print(json.dumps({'all_pass': True, 'exact_projection_layers': [0, 1, 2],
                  'first_graph_direct_min_cap_point_checks': point_checks,
                  'delta_for_direct_checks': str(d),
                  'formal_order': order, 'depth': 5, 'pruning': False, 'rows': rows,
                  'scope': 'Finite exact algebra and implementation checks only; general tail and closure proved analytically.'},
                 indent=2))
