#!/usr/bin/env python3
"""New coefficient/trace guard only; no actual FIRST or weak simulation."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import hashlib
import json
import sys
sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
PREFIX = 'obstacle_positive_source_leak'
checks = []
rounds = []

def rec(label, ok):
    checks.append({'label': label, 'pass': bool(ok)})

def expminus(q):
    M = 4*((q.numerator+q.denominator-1)//q.denominator)+64
    term = F(1)
    low = term
    for k in range(1, M+1):
        term *= q/k
        low += term
    first = term*q/(M+1)
    high = low+first/(1-q/(M+2))
    return 1/high, 1/low

def affine_exp(a, b, elo, ehi):
    return (a+b*elo, a+b*ehi) if b >= 0 else (a+b*ehi, a+b*elo)

def G(u, i):
    return [(a+u[x^(1<<i)])/2 for x, a in enumerate(u)]

def Q(u, n):
    gs = [G(u, i) for i in range(n)]
    return [sum(g[x] for g in gs) for x in range(len(u))]

for round_id, n in enumerate([8, 32, 128]):
    start = len(checks)
    k = (n+1).bit_length()  # exactly ceil(log2(n+2))
    nodes = [F(0), F(1,n*n), F(1,2*n), F(1,n), F(2,n), min(F(4,n),F(1,2)), F(1,2), F(1)]
    samples = []
    for ri, r in enumerate(nodes):
        q = n*r
        elo, ehi = expminus(q)
        b = (1-r)**(n-1)
        d1lo, d1hi = affine_exp(r*b, -r, elo, ehi)
        blo, bhi = affine_exp(((n+1)*r-1)*b, 1-n*r, elo, ehi)
        tag = f'n{n}/r{ri}'
        rec(tag+'/degree1-base', -F(1,n) <= d1lo <= d1hi <= F(1,n))
        rec(tag+'/degree1-beta', -F(18,n) <= blo <= bhi <= F(18,n))
        for t in [F(1), F(k)]:
            tt = tag+f'/t{t.numerator}'
            a1lo, a1hi = affine_exp(b*((1+t*(n+1))*r-t), t-(1+t*n)*r, elo, ehi)
            directlo = d1lo+t*blo
            directhi = d1hi+t*bhi
            rec(tt+'/a1-interval-consistency', max(a1lo,directlo) <= min(a1hi,directhi))
            rec(tt+'/a1-trace', -F(18,n)-F(1,n)/t <= a1lo/t <= a1hi/t <= F(18,n)+F(1,n)/t)
            a2lo, a2hi = affine_exp(-t*r*b, t*r-(1+t*n)*r*r/2, elo, ehi)
            rec(tt+'/a2-nonpositive', a2hi <= 0)
            fact = 2
            for m in range(3,n+3):
                fact *= m
                coefficient = r**(m-1)*(t*m-(1+t*n)*r)/fact
                hi = coefficient*(ehi if coefficient >= 0 else elo)
                rec(tt+f'/repeat{m}', max(hi,F(0)) <= t/F(n**(m-1)))
            cutoff = n+2
            partial = n*sum(t/F(n**(m-1)) for m in range(3,cutoff+1))
            remaining = n*t/F(n**cutoff)/(1-F(1,n))
            rec(tt+'/one-face-infinite-envelope', partial+remaining == t/F(n-1))
            samples.append({'r': str(r), 't': str(t), 'a2_sign': 'nonpositive', 'gamma': str(F(18,n)+F(1,n)/t), 'one_face_L1_mass_factor': str(t/F(n-1))})
    rec(f'n{n}/constant18', F(1113,64) < 18)
    J = isqrt(n)//2
    rec(f'n{n}/short-J-range', 4*J*J <= n)
    rec(f'n{n}/short-exp-bound', F(2*J*J,n) <= F(1,2))
    # Independently verify the new Q^j trace, not any old clipping or
    # resolvent guard. Finite Markov cube is solely an algebra component.
    d = round_id+2
    u = [F(7) if x==0 else F(2) if x==1 else F(0) for x in range(1<<d)]
    qu = Q(u,d)
    su = [d*a-b for a,b in zip(u,qu)]
    cap = max(-a for a in su)
    qj = u[:]
    for j in range(1,12):
        qj = Q(qj,d)
        rec(f'n{n}/finite-d{d}/Q{j}-global', all(qj[x] <= d**j*u[x]+j*d**(j-1)*cap for x in range(len(u))))
        rec(f'n{n}/finite-d{d}/Q{j}-exterior', all(qj[x] <= j*d**(j-1)*cap for x,a in enumerate(u) if a==0))
    rounds.append({'n': n, 'checks': len(checks)-start, 'samples': samples, 'J': J, 'short_degree_trace_factor': 6*J*(J+1)})

registration = HERE/f'{PREFIX}_registration_20261007.json'
result = {'status': 'PASS' if all(c['pass'] for c in checks) else 'FAIL', 'seed': None,
          'scope': 'New original coefficient intervals and finite general trace components only; no actual FIRST/weak endpoint certification.',
          'total_checks': len(checks), 'failed_checks': [c for c in checks if not c['pass']],
          'rounds': rounds, 'checks': checks,
          'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'registration_sha256': hashlib.sha256(registration.read_bytes()).hexdigest()}
(HERE/f'{PREFIX}_results_20261007.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','total_checks','failed_checks','script_sha256']},ensure_ascii=False))
raise SystemExit(0 if result['status']=='PASS' else 1)
