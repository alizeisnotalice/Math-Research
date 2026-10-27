#!/usr/bin/env python3
"""Finite integer-root checks and 16 unrun spacing events; never computes a union."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import factorial, prod
from pathlib import Path
import json

D = 10**40

def iroot8(n):
    a, b = 0, 1
    while b**8 <= n:
        b *= 2
    while b-a > 1:
        c = (a+b)//2
        if c**8 <= n:
            a = c
        else:
            b = c
    assert a**8 <= n < (a+1)**8
    return a

T = tuple(combinations_with_replacement(range(1, 5), 7))
M = sorted({j*prod(t) for j in (1, 2) for t in T} | {96, 192})
R = {m: iroot8(m*D**8) for m in M}

def root_iv(m):
    q = R[m]
    return F(q, D), F(q if q**8 == m*D**8 else q+1, D)

def orbit_mult(t):
    return factorial(7)//prod(factorial(t.count(v)) for v in set(t))

def tail_bounds(j, t):
    # Coordinates in units 1/(4D); a is the fixed tail spacing.
    a = (2*(D+R[2]), 2*(D+R[2]+1))
    rad = (2*R[j*prod(t)], 2*(R[j*prod(t)]+1))
    inner, outer = [], []
    for ell in t:
        q = -(ell//2) if ell % 2 else 1-ell//2
        c = q+ell-1
        lo = (max(0, c*a[0]-rad[1]), max(0, c*a[1]-rad[0]))
        hi = (min(a[0]//2, q*a[1]+rad[0]), min(a[1]//2, q*a[0]+rad[1]))
        inner.append((lo[1], hi[0]))
        outer.append((lo[0], hi[1]))
    return inner, outer

def signature(terms):
    # terms (m,c) mean c*m^(1/8). Extract only exact eighth powers.
    # Equal signatures prove equality; no linear-independence claim is used.
    result = {}
    for m, c in terms:
        aa = bb = 0
        while m % 2 == 0:
            m //= 2
            aa += 1
        while m % 3 == 0:
            m //= 3
            bb += 1
        assert m == 1
        k = (aa % 8, bb % 8)
        result[k] = result.get(k, F(0))+c*2**(aa//8)*3**(bb//8)
    return tuple((a,b,c.numerator,c.denominator) for (a,b),c in sorted(result.items()))

def terms_iv(terms):
    return (sum(c*root_iv(m)[0] for m,c in terms),
            sum(c*root_iv(m)[1] for m,c in terms))

def encode_iv(pair):
    return [str(z) for z in pair]

def main():
    positives = {1: [], 2: []}
    for j in (1,2):
        for t in T:
            inn,out = tail_bounds(j,t)
            ip,op = all(l<h for l,h in inn),all(l<h for l,h in out)
            assert ip == op, (j,t,"unresolved tail-box sign")
            if ip:
                positives[j].append(t)
    P = {j: sorted({j*prod(t) for t in positives[j]}) for j in (1,2)}
    raw = []
    for q in P[2]:
        raw.append(("r",q,None,((q,F(1,2)),)))
    for m in sorted(set(P[1]+P[2])):
        raw.append(("2r",m,None,((m,F(1)),)))
    for p in P[1]:
        for q in P[2]:
            if p<q:
                raw.append(("sum",p,q,((p,F(1,2)),(q,F(1,2)))))
    events = {}
    for kind,p,q,terms in raw:
        key = signature(terms)
        iv = terms_iv(terms)
        alias = {"kind":kind,"p":p,"q":q}
        if key in events:
            ev = events[key]
            ev["interval"] = (max(iv[0],ev["interval"][0]),min(iv[1],ev["interval"][1]))
            assert ev["interval"][0] <= ev["interval"][1]
            ev["aliases"].append(alias)
        else:
            events[key] = {"signature":key,"interval":iv,"aliases":[alias]}
    ordered = sorted(events.values(),key=lambda e:e["interval"][0])
    # These strict separations certify that no unrecognised identities remain.
    for a,b in zip(ordered,ordered[1:]):
        assert a["interval"][1] < b["interval"][0], "unresolved event order"
    ref_terms = ((96,F(1,2)),(192,F(1,2)))
    ref = terms_iv(ref_terms)
    assert signature(ref_terms) not in events
    for ev in ordered:
        lo,hi = ev["interval"]
        if hi < ref[0]:
            ev["distance"] = (ref[0]-hi,ref[1]-lo)
            ev["side"] = "below"
        else:
            assert ref[1] < lo, "reference/event separation unresolved"
            ev["distance"] = (lo-ref[1],hi-ref[0])
            ev["side"] = "above"
    ranked = sorted(ordered,key=lambda e:e["distance"][0])
    for a,b in zip(ranked,ranked[1:]):
        assert a["distance"][1] < b["distance"][0], "unresolved distance order"
    chosen = []
    for rank,ev in enumerate(ranked[:16],1):
        chosen.append({"rank":rank,"aliases":ev["aliases"],"signature":ev["signature"],
            "spacing_interval":encode_iv(ev["interval"]),
            "distance_from_s96_interval":encode_iv(ev["distance"]),
            "side":ev["side"],"display_spacing":float(sum(ev["interval"])/2),
            "union_status":"NOT RUN","finite_L1_status":"NOT RUN"})
    predecessor = max((ev for ev in ordered if ev["interval"][1]<ref[0]),key=lambda e:e["interval"][0])
    successor = min((ev for ev in ordered if ref[1]<ev["interval"][0]),key=lambda e:e["interval"][0])
    output = {"scope":"fixed tail a=(1+2^(1/8))/2; labels {1,2} x {1,2,3,4}^7",
        "integer_root_scale_D":str(D),"root_floor_numerators":{str(m):str(R[m]) for m in M},
        "all_root_floor_inequalities_passed":True,
        "tail_orbits_checked":2*len(T),"raw_labels_represented":2*sum(orbit_mult(t) for t in T),
        "positive_tail_orbits":{str(j):len(positives[j]) for j in (1,2)},
        "positive_tail_raw_labels":{str(j):sum(orbit_mult(t) for t in positives[j]) for j in (1,2)},
        "positive_tail_mass_sets":{str(j):P[j] for j in (1,2)},
        "raw_safe_events":len(raw),"distinct_safe_events":len(ordered),
        "all_distinct_event_orders_strictly_certified":True,
        "all_distance_ranks_strictly_certified":True,
        "terminal_saturation_spacing":"32768^(1/8)",
        "terminal_saturation_interval":encode_iv(root_iv(32768)),
        "s96_interval":encode_iv(ref),"s96_is_safe_event":False,
        "s96_containing_open_event_cell":[predecessor["aliases"],successor["aliases"]],
        "candidate_generation":"16 distinct safe events of smallest absolute distance from s96; exact eighth-power identities merged before interval ranking",
        "candidates":chosen,"union_evaluations":0,"finite_L1_evaluations":0}
    target = Path(__file__).with_name("event_certificate.json")
    target.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({k:output[k] for k in ("positive_tail_orbits","positive_tail_raw_labels","raw_safe_events","distinct_safe_events","s96_containing_open_event_cell","union_evaluations")},indent=2))
    print("PASS: integer root, tail sign, all event order and all neighbour rank checks")

if __name__ == "__main__":
    main()
