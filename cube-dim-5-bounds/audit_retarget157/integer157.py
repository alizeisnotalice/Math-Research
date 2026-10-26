"""Independent153 integer coefficients and one-clear-denominator corners."""
from math import gcd
from rational157 import Q


def lcm(a, b):
    return a//gcd(a, b)*b


def roots():
    return [[a,b,c,d] for a in range(32) for b in range(a+1)
            for c in range(b+1) for d in range(c+1)]


def bounds(node, table):
    rid, path = node
    assert type(rid) is int and 0 <= rid < 52360 and not set(path)-{'0','1'}
    lo, den = list(table[rid]), [32]*4
    for j, bit in enumerate(path):
        k = j%4
        lo[k], den[k] = 2*lo[k]+int(bit), 2*den[k]
    return [(lo[k],lo[k]+1,den[k]) for k in range(4)]


def inside(point, node, table):
    assert Q(point[0]) == 1
    return all(Q(a,d) <= Q(x) <= Q(b,d)
               for (a,b,d),x in zip(bounds(node,table),point[1:]))


def coefficient(source):
    c = source['old_core']
    e, old, upper = c['epsilon'], Q(c['penalty_coefficient']), Q(source['old_H_interval'][1])
    assert upper > Q(0) and old > Q(0)
    mass = source['mass_qualification']
    if e == '1/64':
        assert old == 256 and mass is None
        A, lower = 256, Q(0)
    else:
        assert e == '1/128'
        if mass is None:
            assert c['curvature_rule_class'] == '075_full_frame' and old == 1024
        A, lower = 1024, Q(0) if mass is None else Q(mass['lower'])
    assert Q(0) <= lower < upper
    # Common denominator30: 7/10-8r/15=(21-16r)/30.
    split = Q(A*(21*upper.n*lower.d-16*lower.n*upper.d),30*upper.n*lower.d)
    uniform = Q(2*A,3)
    q = min(old,uniform,split)
    assert Q(0) < q <= uniform and q <= old
    return q, {'A':Q(A),'mass_ratio_lower':lower/upper,
               'effective_mass_ratio_lower':max(Q(1,16),lower/upper),
               'uniform_time_cap':uniform,'split_time_mass_cap':split,
               'old_q_upper':old,'new_q_upper':q}


def corners(core, q, leaf, log_upper, table):
    b = bounds(leaf['node'],table)
    for (a,z,d),lo,hi in zip(b,leaf['box_lower'],leaf['box_upper']):
        assert Q(a,d) == Q(lo) and Q(z,d) == Q(hi)
    point = [Q(x) for x in core['point'][1:]]
    gradient = [[Q(a),Q(z)] for a,z in core['loggradient_directed']]
    delta = Q(core['Delta110_lower'])
    W = 1
    for x in point:
        W = lcm(W,x.d)
    for a,z,d in b:
        W = lcm(W,d)
    anchor = [x.n*(W//x.d) for x in point]
    ib = [(a*(W//d),z*(W//d)) for a,z,d in b]
    T = lcm(delta.d,q.d*W*W)
    for a,z in gradient:
        assert a <= z
        T = lcm(T,a.d*W)
        T = lcm(T,z.d*W)
    d0 = delta.n*(T//delta.d)
    qi = q.n*(T//(q.d*W*W))
    gi = [(a.n*(T//(a.d*W)),z.n*(T//(z.d*W))) for a,z in gradient]
    best = None
    for mask in range(16):
        bits = [(mask>>(3-j))&1 for j in range(4)]
        v = [ib[j][bits[j]]-anchor[j] for j in range(4)]
        value = d0+sum(min(a*x,z*x) for (a,z),x in zip(gi,v))-qi*sum(x*x for x in v)
        if best is None or value < best[0]:
            best = value,bits
    margin110 = Q(best[0],T)
    assert margin110 >= Q(leaf['saved_g110_margin']) > Q(0)
    return margin110,margin110-log_upper,best[1]


def public_log(certificate):
    """Positive polynomial integral, one integer denominator for64 terms."""
    assert certificate['ratio'] == '23/22' and certificate['terms'] == 64
    denominators = [(2*k+1)*45**(2*k+1) for k in range(64)]
    T = 1
    for d in denominators:
        T = lcm(T,d)
    lower = Q(sum(2*(T//d) for d in denominators),T)
    # Integral remainder 2∫0^(1/45) s^128/(1-s²)ds.
    tail = Q(2*45*45,129*45**129*(45*45-1))
    assert Q(certificate['lower']) == lower
    assert Q(certificate['tail_upper']) == tail
    upper = Q(certificate['upper'])
    assert upper == lower+tail and upper > Q(0)
    return upper
