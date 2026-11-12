"""Gamma band probabilities via the positive lower-tail exponential series."""
from decimal import Decimal as D, getcontext
from math import factorial
getcontext().prec = 70
n = 64
def lower(t):
    t = D(t)
    term = t**n / D(factorial(n))
    total = term
    j = n
    while term:
        j += 1
        term *= t / j
        updated = total + term
        if updated == total:
            break
        total = updated
    return (-t).exp() * total
def prob(a,b):
    return lower(b)-lower(a)
def show(label,s,d,gamma,intervals):
    mass = sum(prob(a,b) for a,b in intervals)
    print(label,*(format(D(x),'.12g') for x in (s,d,gamma,mass,2**n*mass)))
w = D('0.00002')
show('corner',1,0,1,[(0,n*w)])
show('mode',prob(62,66),0,1,[(62,66)])
gamma = D(8)/9
show('mixed',gamma*prob(52,63),gamma*prob(65,76),gamma,[(116,127),(129,140)])
