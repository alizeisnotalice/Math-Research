"""Integer outward arithmetic; no numerical or production library is imported."""
from fractions import Fraction as F
from math import lcm, prod
from itertools import product

P = 768
S = 1 << P
LOG_TERMS = 282
ATAN_TERMS = 216


def ceildiv(n, d):
    assert d > 0
    return -((-n) // d)


def atanh_log_bounds(n, d):
    """2atanh(n/d), 0<=n/d<=1/3, returned as integer multiples of2^-P."""
    assert 0 <= 3*n <= d
    low, high = S*n//d, ceildiv(S*n, d)
    sq_low, sq_high = low*low//S, ceildiv(high*high, S)
    power_low, power_high = low, high
    total_low = total_high = 0
    for j in range(LOG_TERMS):
        odd = 2*j + 1
        total_low += power_low//odd
        total_high += ceildiv(power_high, odd)
        power_low = power_low*sq_low//S
        power_high = ceildiv(power_high*sq_high, S)
    # 2sum_{j>=N}y^(2j+1)/(2j+1) <=9/[4(2N+1)3^(2N+1)].
    odd = 2*LOG_TERMS + 1
    tail = ceildiv(9*S, 4*odd*3**odd)
    return 2*total_low, 2*total_high + tail


LOG_TWO = atanh_log_bounds(1, 3)


def log_bounds(n, d=1):
    assert n > 0 and d > 0
    exponent = n.bit_length() - d.bit_length()
    if exponent >= 0:
        if n < d << exponent:
            exponent -= 1
    elif n << (-exponent) < d:
        exponent -= 1
    a, b = (n, d << exponent) if exponent >= 0 else (n << (-exponent), d)
    assert b <= a < 2*b
    lo, hi = atanh_log_bounds(a-b, a+b)
    if exponent >= 0:
        return lo + exponent*LOG_TWO[0], hi + exponent*LOG_TWO[1]
    return lo + exponent*LOG_TWO[1], hi + exponent*LOG_TWO[0]


def atan_inverse_bounds(q):
    low = high = 0
    power = q
    for j in range(ATAN_TERMS):
        d = (2*j+1)*power
        a, b = S//d, ceildiv(S, d)
        if j % 2:
            low -= b; high -= a
        else:
            low += a; high += b
        power *= q*q
    tail = ceildiv(S, (2*ATAN_TERMS+1)*power)
    if ATAN_TERMS % 2:
        low -= tail
    else:
        high += tail
    return low, high


A5, A239 = atan_inverse_bounds(5), atan_inverse_bounds(239)
PI = (16*A5[0]-4*A239[1], 16*A5[1]-4*A239[0])
assert 3*S < PI[0] < PI[1] < 4*S
# Exact tangent identity: (5+i)^4*(239-i)=114244*(1+i).
assert 476*239+480 == 480*239-476 == 114244


def add_product(acc, numerator, denominator, interval):
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    assert denominator > 0
    low, high = interval
    if numerator < 0:
        low, high = high, low
    acc[0] += numerator*low//denominator
    acc[1] += ceildiv(numerator*high, denominator)


def direct_values(point, fan):
    """All3504 cones: H and all4 physical derivatives via direct rational poles."""
    point = list(map(F, point));D = lcm(*(v.denominator for v in point))
    integer_point = [int(v*D) for v in point]
    vertices = fan['vertices'];cones = fan['simplices']
    assert fan['denominator'] == 144 and len(vertices) == 280 and len(cones) == 3504
    k = [sum(x*y for x,y in zip(v,integer_point)) for v in vertices]
    assert len(set(k)) == 280 and 0 not in k
    a = F(144*D,64);A,B = a.numerator,a.denominator
    squared = [A*A+B*B*x*x for x in k]
    logs = [log_bounds(x) for x in squared]
    slopes = [[D*x for x in vertex[1:]] for vertex in vertices]
    totals = [[0,0] for _ in range(5)]
    for cone in cones:
        ids,det = cone['vertices'],cone['det_numerator']
        assert len(ids) == len(set(ids)) == 5 and det > 0
        for vertex in ids:
            others = [m for m in ids if m != vertex]
            differences = [k[vertex]-k[m] for m in others]
            denominator = prod(differences)
            add_product(totals[0], -det, 2*denominator, logs[vertex])
            quotients = [denominator//d for d in differences]
            for axis in range(4):
                # sigma = sum((s_j-s_m)/(k_j-k_m)) = T/denominator.
                T = sum((slopes[vertex][axis]-slopes[m][axis])*q for m,q in zip(others,quotients))
                add_product(totals[axis+1], det*T, 2*denominator*denominator, logs[vertex])
                add_product(totals[axis+1], -det*B*B*slopes[vertex][axis]*k[vertex], denominator*squared[vertex], (S,S))
    factor_numerator = D**4*S**4
    factor_denominators = (4608*PI[0]**5, 4608*PI[1]**5)
    output=[]
    for low,high in totals:
        # All four products are safe for either sign of a derivative interval.
        endpoints = [F(factor_numerator*v,d) for v in (low,high) for d in factor_denominators]
        output.append((min(endpoints),max(endpoints)))
    assert output[0][0] > 0
    return output, {'common_denominator':str(D),'all_cones':len(cones),'distinct_projections':len(set(k)),
                    'minimum_integer_projection_gap':str(min(y-x for x,y in zip(sorted(k),sorted(k)[1:]))) }


def dyadic_outward(interval, bits=512):
    scale=1<<bits;low,high=interval
    a=(low*scale).__floor__();b=(high*scale).__ceil__()
    return str(F(a,scale)),str(F(b,scale))
