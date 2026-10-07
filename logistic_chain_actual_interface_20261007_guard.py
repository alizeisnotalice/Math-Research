#!/usr/bin/env python3
"""Finite exact representation guard; not an original-space or actual-gate test."""
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PREFIX = 'logistic_chain_actual_interface_20261007'
REG = HERE / (PREFIX + '_registration.json')
OUT = HERE / (PREFIX + '_results.json')
if OUT.exists():
    raise SystemExit('Refusing to overwrite a completed guard.')
reg = json.loads(REG.read_text())
assert reg['rounds'] == [4, 5, 6]
assert reg['source_masses'] == [2, 3]
assert reg['Ch'] == '16/3'
assert reg['registered_before_execution'] is True

def pgf(probs):
    out = [F(1)]
    for p in probs:
        nxt = [F(0)] * (len(out) + 1)
        for k, value in enumerate(out):
            nxt[k] += value * (1-p)
            nxt[k+1] += value * p
        out = nxt
    return out

def fmt(q):
    return f'{q.numerator}/{q.denominator}'

round_results = []
for cells in reg['rounds']:
    # Keep the explicit count in one wrapper; no floats or tolerance.
    def verify(ok):
        nonlocal_box[0] += 1
        if not ok:
            raise AssertionError(f'round cells={cells}, assertion={nonlocal_box[0]}')
    nonlocal_box = [0]
    probs = [F(1, i+2) for i in range(cells)]
    sets = [{0}, {0, 1} | {i for i in range(4, cells) if i % 2 == 0},
            {0, 2} | {i for i in range(4, cells) if i % 2 == 1}, set(range(cells))]
    states = []
    for bits in product([0, 1], repeat=cells):
        probability = F(1)
        for i, bit in enumerate(bits):
            probability *= probs[i] if bit else 1-probs[i]
        states.append((bits, probability))
    verify(sum(p for _, p in states) == 1)
    for selected in sets:
        empirical = [F(0)] * (len(selected)+1)
        for bits, p in states:
            empirical[sum(bits[i] for i in selected)] += p
        expected = pgf([probs[i] for i in sorted(selected)])
        for left, right in zip(empirical, expected):
            verify(left == right)
    # Counts in an included set and in its disjoint new cells are independent.
    for a, b in [(0,1), (0,2), (0,3), (1,3), (2,3)]:
        verify(sets[a] <= sets[b])
        old, new = sets[a], sets[b]-sets[a]
        joint = defaultdict(F)
        for bits, p in states:
            joint[sum(bits[i] for i in old), sum(bits[i] for i in new)] += p
        po, pn = pgf([probs[i] for i in sorted(old)]), pgf([probs[i] for i in sorted(new)])
        for i in range(len(po)):
            for j in range(len(pn)):
                verify(joint[i,j] == po[i]*pn[j])
    intersection = sets[1] & sets[2]
    left_only, right_only = sets[1]-intersection, sets[2]-intersection
    verify(intersection == {0})
    verify(not sets[1] <= sets[2] and not sets[2] <= sets[1])
    joint = defaultdict(F)
    for bits, p in states:
        joint[sum(bits[i] for i in intersection), sum(bits[i] for i in left_only),
              sum(bits[i] for i in right_only)] += p
    polynomials = [pgf([probs[i] for i in sorted(s)])
                   for s in [intersection, left_only, right_only]]
    for i, j, k in product(*[range(len(v)) for v in polynomials]):
        verify(joint[i,j,k] == polynomials[0][i]*polynomials[1][j]*polynomials[2][k])
    mean_left = sum(probs[i] for i in sets[1])
    mean_right = sum(probs[i] for i in sets[2])
    mixed = sum(p*sum(bits[i] for i in sets[1])*sum(bits[i] for i in sets[2]) for bits,p in states)
    verify(mixed-mean_left*mean_right == sum(probs[i]*(1-probs[i]) for i in intersection))

    source_mass, W, Ch = [F(2), F(3)], F(5), F(16,3)
    def Q(y, bits, j, x):
        count = sum(bits[i] for i in sets[j])
        return F((x+y+count) % 4 + 1, 10)
    def acceptance(y, j):
        return F(2+y+j, 7+j)
    kernel = [[sum(p*Q(y,bits,j,j) for bits,p in states) for j in range(4)] for y in range(2)]
    for y in range(2):
        for bits, _ in states:
            for j in range(4):
                verify(sum(Q(y,bits,j,x) for x in range(4)) == 1)
                verify(0 < acceptance(y,j) < 1)
    T = [Ch*sum(source_mass[y]*kernel[y][j]*acceptance(y,j) for y in range(2)) for j in range(4)]
    tau = min(T)/2
    weak_g = [tau/t for t in T]
    verify(all(t > tau for t in T))
    verify(all(0 < g < 1 for g in weak_g))
    mode_results = {}
    for mode, weights in [('strong', [F(1)]*4), ('weak', weak_g)]:
        count_pgf = [F(0)]*5
        occurrences = [F(0)]*4
        verify(sum(source_mass[y]/W for y in range(2)) == 1)
        for y in range(2):
            for bits, p in states:
                mixing = source_mass[y]/W*p
                coin_probs = [Q(y,bits,j,j)*acceptance(y,j)*weights[j] for j in range(4)]
                verify(all(0 <= q <= 1 for q in coin_probs))
                conditional = pgf(coin_probs)
                verify(sum(conditional) == 1)
                for k in range(5):
                    count_pgf[k] += mixing*conditional[k]
                for j in range(4):
                    occurrences[j] += mixing*coin_probs[j]
        mean = sum(F(k)*count_pgf[k] for k in range(5))
        stop = 1-count_pgf[0]
        deficit = sum(F(k-1)*count_pgf[k] for k in range(2,5))
        verify(sum(count_pgf) == 1)
        verify(sum(occurrences) == mean)
        verify(mean == sum(weights[j]*T[j] for j in range(4))/(Ch*W))
        verify(mean == stop+deficit)
        verify(deficit == sum(sum(count_pgf[k] for k in range(m,5)) for m in range(2,5)))
        verify(0 <= stop <= 1)
        palm = [F(k)*count_pgf[k]/mean for k in range(5)]
        verify(sum(palm) == 1)
        verify(sum(palm[k]/k for k in range(1,5)) == stop/mean)
        verify(sum(palm[k]*F(k-1,k) for k in range(1,5)) == deficit/mean)
        if mode == 'weak':
            verify(mean == 4*tau/(Ch*W))
            for j in range(4):
                verify(occurrences[j]/mean == F(1,4))
        else:
            for j in range(4):
                verify(occurrences[j]/mean == T[j]/sum(T))
        mode_results[mode] = {'mean_count':fmt(mean), 'source_once_stop_probability':fmt(stop),
                              'unpaid_deficit':fmt(deficit), 'count_pgf':[fmt(p) for p in count_pgf]}
    round_results.append({'cells':cells, 'field_patterns':2**cells, 'checks':nonlocal_box[0],
                          'tau':fmt(tau), 'Ch':fmt(Ch), 'source_mass':fmt(W), 'modes':mode_results,
                          'status':'PASS'})

output = {'id':reg['id'], 'status':'PASS', 'rounds':round_results,
          'total_checks':sum(r['checks'] for r in round_results), 'scope':reg['scope'],
          'actual_gate_sample':False, 'central_order_sample':False, 'paid_fee':False,
          'registration_sha256':hashlib.sha256(REG.read_bytes()).hexdigest(),
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
OUT.write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({'status':output['status'], 'round_checks':[r['checks'] for r in round_results],
                  'total_checks':output['total_checks'], 'scope':reg['scope']}, ensure_ascii=False))
