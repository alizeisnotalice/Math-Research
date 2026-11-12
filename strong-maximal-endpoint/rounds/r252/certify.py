from fractions import Fraction as Q
from functools import lru_cache

def series(x):
    z = (x-1)/(x+1)
    low = 2*sum(z**(2*j+1)/(2*j+1) for j in range(12))
    return low, low + 2*z**25/(25*(1-z*z))

ln2lo, ln2hi = series(Q(2))
scale = 10**12

@lru_cache(None)
def logs(x):
    shifts = 0
    while x < 1:
        x *= 2
        shifts += 1
    assert 1 <= x <= 2
    lo, hi = series(x)
    return (Q((scale*(lo-shifts*ln2hi))//1, scale),
            Q(-((-scale*(hi-shifts*ln2lo))//1), scale))

w, a, b, h = Q(1,2), Q(31,100), Q(19,25), Q(1,2000)
def first(t, idx):
    return t*(logs(b*t)[idx]-1) if t else Q(0)
def second(t, idx):
    return (1-t)*logs(a+b*t)[idx]

upper = -2*h*sum((1-i*h)*logs(a+b*i*h)[0] for i in range(1000))
peak = upper
for i in range(500):
    l, r = i*h, (i+1)*h
    increment = 2*h*(first(l,1)+second(r,1)-first(w-l,0)-second(w-r,0))
    peak = max(peak,upper+max(0,increment))
    upper += increment
assert peak < Q(572,1000)
assert ln2hi < Q(6932,10000)
print('Certified upper bound', peak, 'approximately', float(peak))
