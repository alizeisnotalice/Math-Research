"""Finite exact audit checks; no optimizer, search, or operator-bound evaluation."""
from fractions import Fraction as F
from math import comb, prod
import json
import platform


def s(x):
    return str(x)


# Coefficients of q^(6-j)(1-t^(-j)), obtained independently by
# expanding (q-1/(2*rho))^6 and integrating rho^(-j-1).
derived = [F(comb(6, j), j) * F(-1, 2) ** j for j in range(1, 7)]
printed = [F(-3), F(15, 8), F(-5, 6), F(15, 64), F(-3, 80), F(1, 384)]
assert derived == printed


def n1_direct(u, c, t):
    # Direct z integration: positive z never contributes. On [-1,0],
    # the surviving segment ends at u-c. On z=-rho<-1, solve
    # u+rho >= c*rho, then integrate (2*rho)^(-1).
    ball = max(F(0), min(F(0), u-c) + 1)
    cutoff = min(t, u / (c-1))
    return ball / 2, max(F(1), cutoff)


def n1_shell(u, c, t):
    # Independently evaluate L(1), and N(rho) on its only breakpoint.
    ball = 1 - min(F(2), 1+c-u) / 2
    cuts = sorted(set([F(1), t] + ([u/(c-1)] if 1 < u/(c-1) < t else [])))
    contributing = []
    for lo, hi in zip(cuts, cuts[1:]):
        rho = (lo+hi)/2
        count = 1 + int(u < (c-1)*rho)
        if 2-count:
            contributing.append(hi/lo)
    # Each surviving shell has coefficient 1/2 in front of log.
    return ball, prod(contributing, start=F(1))


cases = [
    (F(0), F(5,4), F(2)),
    (F(1,4), F(5,4), F(2)),  # threshold rho=1
    (F(1,2), F(5,4), F(2)),  # threshold rho=t
    (F(3,4), F(5,4), F(2)),  # entire shell contributes
    (F(3,5), F(7,5), F(2)),  # threshold rho=3/2 inside shell
    (F(9,10), F(2), F(2)),   # empty cap
]
n1_results = []
for u, c, t in cases:
    direct = n1_direct(u,c,t)
    shell = n1_shell(u,c,t)
    assert direct == shell
    n1_results.append({"u":s(u),"c":s(c),"t":s(t),
                       "H_form":"ball + (1/2)*log(log_argument)",
                       "ball":s(direct[0]),"log_argument":s(direct[1]),"match":True})


# Fixed face checks, including a face exactly at a threshold. The direct
# routine clips each OTHER coordinate to the two interval endpoints and
# separately checks each signed face against the open capture condition.
face_cases = [
    (F(1), F(5,4), [F(1,4)]),
    (F(3,2), F(4,3), [F(0),F(1,2)]),
    (F(1), F(5,4), [F(0),F(1,4),F(1,2),F(3,4),F(1,8),F(7,8),F(1,3)]),
]
face_results=[]
for rho,c,u in face_cases:
    formula=F(0)
    direct=F(0)
    for i in range(len(u)):
        formula += (1+int(u[i]<(c-1)*rho))*prod(
            (min(2*rho,(1+c)*rho-u[j]) for j in range(len(u)) if j!=i), start=F(1))
        for sign in [-1,1]:
            if abs(sign*rho-u[i]) < c*rho:
                direct += prod((max(F(0), min(rho,u[j]+c*rho)-max(-rho,u[j]-c*rho))
                                for j in range(len(u)) if j!=i), start=F(1))
    assert formula == direct
    face_results.append({"n":len(u),"rho":s(rho),"c":s(c),"u":list(map(s,u)),
                         "captured_face_sum":s(direct),"match":True})

print(json.dumps({"python":platform.python_version(),"arithmetic":"standard-library fractions.Fraction (exact)",
    "sixth_power_coefficients":list(map(s,derived)),"sixth_power_formula_match":True,
    "n1_direct_vs_shell":n1_results,"fixed_face_checks":face_results,
    "scope":"finite audit identities only; no numerical optimizer, parameter scan, scalar bound, or new research"}, indent=2))
