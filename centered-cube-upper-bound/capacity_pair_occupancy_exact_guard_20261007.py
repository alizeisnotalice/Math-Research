#!/usr/bin/env python3
"""Exact finite mass-tree guards; not an actual FIRST/geometric sample.

Only the algebra connecting original source-fixed coins to conditional failure
acceptance is tested. The Lebesgue-column receipt estimate is an upstream proof,
not supplied by these finite surrogate pair weights. No randomness is used.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib
import json


CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def rat(value):
    value = F(value)
    return f"{value.numerator}/{value.denominator}"


def coin(numerator, denominator):
    return F(1) if denominator == 0 else min(F(1), numerator / denominator)


def ceil_log2(k):
    return (k - 1).bit_length()


def threshold_tests(records, thresholds):
    """Records are (p,r,positive base weight,conditional failing acceptance)."""
    union = sum((weight * (p + r - p*r) for p, r, weight, g in records), F(0))
    failure = sum((weight * (1-p)*(1-r)*g for p, r, weight, g in records), F(0))
    result = []
    for delta in thresholds:
        high = low = F(0)
        high_count = low_count = equality_count = 0
        for p, r, weight, g in records:
            ell = p + r - p*r
            w = (1-p)*(1-r)
            check(w == 1-ell, "joint union/failure")
            check(p*(1-r)+(1-p)*r+p*r == ell, "independent union outcomes")
            check(w+p*(1-r)+(1-p)*r+p*r == 1, "four outcome mass")
            check(0 <= g*w <= w, "original failing gate domination")
            if ell >= delta:
                high += weight*g*w
                high_count += 1
                equality_count += int(ell == delta)
                check(g*w <= (1-delta)*ell/delta, "pointwise high fee")
            else:
                low += weight*g*w
                low_count += 1
                check(p < delta and r < delta and w > 1-delta, "strict low contract")
        check(high+low == failure, "complementary current failure split")
        check(high <= (1-delta)*union/delta, "summed high amplification")
        check(union+high <= union/delta, "consolidated original union plus high failure")
        result.append({"delta": rat(delta), "union": rat(union),
                       "failure": rat(failure), "high_failure": rat(high),
                       "low_failure": rat(low), "high_count": high_count,
                       "low_count": low_count, "equality_count": equality_count})
    return result


def tree_model(n, arity, depth, skew, delta):
    rootn = isqrt(n)
    leaves = arity**depth
    raw = [F(2**(i % arity) if skew else 1) for i in range(leaves)]
    total = sum(raw)
    m = [v/total for v in raw]
    h = [mass if not skew or i % 3 else F(0) for i, mass in enumerate(m)]
    whole = sum(m)
    born = sum(h)
    check(whole == 1 and 0 < born <= whole, "whole source/born subsource")
    nodes = {}
    for level in range(depth+1):
        size = arity**level
        for node in range(leaves//size):
            lo, hi = node*size, (node+1)*size
            nodes[level, node] = (sum(m[lo:hi]), sum(h[lo:hi]))
    cap_rows = []
    for level in range(1, depth+1):
        forward = reverse = F(0)
        for parent in range(leaves//arity**level):
            parent_m, parent_h = nodes[level, parent]
            child_msum = child_hsum = F(0)
            for child in range(parent*arity, (parent+1)*arity):
                cm, ch = nodes[level-1, child]
                child_msum += cm
                child_hsum += ch
                p = coin(cm, rootn*(parent_h-ch))
                r = coin(ch, n*(parent_m-cm))
                forward += p*(parent_h-ch)
                reverse += r*(parent_m-cm)
                check(p*(parent_h-ch) <= cm/rootn, "child forward capacity")
                check(r*(parent_m-cm) <= ch/n, "child reverse capacity")
            check(child_msum == parent_m and child_hsum == parent_h, "full parent source partition")
        check(forward <= whole/rootn and reverse <= born/n, "level source-once capacities")
        cap_rows.append({"level": level, "forward": rat(forward), "reverse": rat(reverse)})
    records = []
    for i in range(leaves):
        for j in range(leaves):
            if i == j or h[j] == 0:
                continue  # outside legal distinct-child positive-source pair domain
            level = next(level for level in range(1, depth+1)
                         if i//arity**level == j//arity**level)
            parent = i//arity**level
            yi, zj = i//arity**(level-1), j//arity**(level-1)
            check(yi != zj, "unique LCA has distinct children")
            pm, ph = nodes[level, parent]
            ym, yh = nodes[level-1, yi]
            zm, zh = nodes[level-1, zj]
            p = coin(ym, rootn*(ph-yh))
            r = coin(zh, n*(pm-zm))
            ell = p+r-p*r
            if ell < delta:
                check(ym < delta*rootn*(ph-yh), "low full soft child mass")
                check(zh < delta*n*(pm-zm), "low full hard child mass")
            weight = m[i]*h[j]*F((i+3*j)%5+1, 5)
            g = F((2*i+j)%7, 6)  # pair-correlated gate, no independence assertion
            records.append((p, r, weight, g))
    endpoint = min(p+r-p*r for p, r, _, _ in records)
    results = threshold_tests(records, [delta, endpoint, F(1)])
    check(results[1]["equality_count"] > 0, "threshold equality assigned paid")
    return {"model": "skew partial-born tree" if skew else "uniform full-born tree",
            "arity": arity, "depth": depth, "whole_mass": rat(whole),
            "born_mass": rat(born), "ordered_positive_distinct_pairs": len(records),
            "capacities": cap_rows, "thresholds": results}


def compressed_star_and_chain(n, delta):
    b, depth = 512, 4
    rootn = isqrt(n)
    p, r = F(1, rootn*(b-1)), F(1, n*(b-1))
    ell, w = p+r-p*r, (1-p)*(1-r)
    check(ell < delta, "512-child model strictly low at main threshold")
    check(b*p*F(b-1,b) == F(1,rootn), "star exact forward capacity")
    check(b*r*F(b-1,b) == F(1,n), "star exact reverse capacity")
    same = b*b//2-b
    opposite = b*b//2
    check(same+opposite == b*(b-1), "exact ordered pair multiplicities")
    # Two exact orbit classes, no sampling or huge pair enumeration.
    records = [(p,r,F(same,b*b),F(1,2)),
               (p,r,F(opposite,b*b),F(1,3))]
    thresholds = threshold_tests(records, [delta, ell, F(1)])
    accepted = F(same,b*b)/2 + F(opposite,b*b)/3
    check(accepted == F(5,12)-F(1,2*b), "correlated gate orbit sum")
    rings = []
    for j in range(1, depth+1):
        child = F(1,b**(depth+1-j))
        ring = (b-1)*child
        parent = b*child
        check(child+ring == parent, "chain full parent source partition")
        check(coin(child,rootn*ring) == p and coin(child,n*ring) == r,
              "ancestor-independent coins with full masses")
        rings.append(ring)
    leaf = F(1,b**depth)
    check(leaf+sum(rings) == 1, "disjoint ring source normalization")
    check(sum(rings) == 1-leaf, "all sibling rings counted once")
    check(w*sum(rings) > F(99,100), "low coins do not supply small ring factor")
    return {"model": "exactly compressed 512-child star and depth-4 chain",
            "occupied_children": b, "chain_depth": depth,
            "p": rat(p), "r": rat(r), "ell": rat(ell), "w": rat(w),
            "strict_main_low": True, "star_thresholds": thresholds,
            "ring_masses": [rat(x) for x in rings], "inner_leaf_mass": rat(leaf),
            "ring_total": rat(sum(rings)), "failure_weighted_ring_total": rat(w*sum(rings))}


def main():
    rounds = []
    for n in [1024,4096,16384]:
        before = CHECKS
        rootn = isqrt(n)
        check(rootn*rootn == n, "exact square dimension")
        k = ceil_log2(n+2)
        check(2**(k-1) < n+2 <= 2**k, "exact predetermined binary log ceil")
        delta = F(1,k**4)
        models = [tree_model(n,2,5,False,delta), tree_model(n,4,2,True,delta),
                  compressed_star_and_chain(n,delta)]
        rounds.append({"n": n, "sqrt_n": rootn, "k_n": k, "delta_n": rat(delta),
                       "models": models, "exact_checks": CHECKS-before, "status": "PASS"})
    result = {"scope": "finite source-tree mass/union/failure algebra only; upstream spatial receipt proof is not simulated",
              "no_actual_FIRST_or_history_certification": True,
              "no_actual_geometric_counterexample": True,
              "seed": "none; deterministic exact Fraction arithmetic",
              "old_epsilon_unchanged": "ln(n+2)^(-4); current coin delta is a separate ceil(log2(n+2))^(-4)",
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "rounds": rounds, "total_exact_checks": CHECKS, "status": "PASS"}
    output = Path(__file__).with_name(Path(__file__).stem.replace("_exact_guard", "_exact_guard_results")+".json")
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status": "PASS", "rounds": 3, "checks": CHECKS,
                      "output": str(output), "scope": result["scope"]},ensure_ascii=False))


if __name__ == "__main__":
    main()
