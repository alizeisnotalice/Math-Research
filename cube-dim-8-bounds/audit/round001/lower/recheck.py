#!/usr/bin/env python3
"""Independent full-label, fixed-coordinate exact event sweep for n8 R001.

No imports or function calls from candidate production/recheck code.
Only Python standard library; stdout is the exact result; input files are read-only.
"""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from itertools import product
from math import isqrt, prod
from pathlib import Path
import hashlib
import json
import resource
import time

HERE = Path(__file__).resolve().parent
D = 10**35
S = 4 * D
LIMIT = 2**22

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def root8(m):
    scaled = m * D**8
    q = isqrt(isqrt(isqrt(scaled)))
    assert q**8 <= scaled < (q + 1)**8
    return q

def decimal(q):
    with localcontext() as ctx:
        ctx.prec = 65
        return str(Decimal(q.numerator) / Decimal(q.denominator))

class LimitReached(RuntimeError):
    pass

class FixedAxisSweep:
    """Each state is the set of physical boxes on the remaining fixed axes.

    Sweep first remaining axis. Event multiplicities maintain the union of
    suffix boxes. Fubini reduces to the next fixed axis; no orbit compression.
    """
    def __init__(self, budget):
        self.memo = {}
        self.budget = budget
        self.max_state_boxes = 0

    def volume(self, state):
        if not state:
            return 0
        if state in self.memo:
            return self.memo[state]
        if self.budget[0] >= LIMIT:
            raise LimitReached("2^22 total memo states exhausted; incomplete")
        self.budget[0] += 1
        self.max_state_boxes = max(self.max_state_boxes, len(state))
        if len(state[0]) == 2:
            total = 0
            left, right = state[0]
            for lo, hi in state[1:]:
                if lo > right:
                    total += right - left
                    left, right = lo, hi
                else:
                    right = max(right, hi)
            total += right - left
        else:
            events = defaultdict(list)
            for box in state:
                events[box[0]].append((box[2:], 1))
                events[box[1]].append((box[2:], -1))
            active = {}
            previous = None
            total = 0
            for x in sorted(events):
                if previous is not None and active:
                    child = tuple(sorted(active))
                    total += (x - previous) * self.volume(child)
                for suffix, delta in events[x]:
                    count = active.get(suffix, 0) + delta
                    assert count >= 0
                    if count:
                        active[suffix] = count
                    else:
                        del active[suffix]
                previous = x
            assert not active
        self.memo[state] = total
        return total

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    frozen = json.loads((HERE / "frozen_inputs.json").read_text())
    source = args.candidate_root
    def verify_inputs():
        for name, digest in frozen["candidate_files"].items():
            assert sha(source / name) == digest, "changed input " + name
        assert sha(args.reference) == frozen["reference_sha256"]
    verify_inputs()
    manifest = json.loads((source / "manifest.json").read_text())
    assert manifest["files"] == {k: v for k, v in frozen["candidate_files"].items() if k != "manifest.json"}
    saved = json.loads((source / "results.json").read_text())
    registration = json.loads((source / "registration.json").read_text())
    old = json.loads(args.reference.read_text())
    assert saved["reference_sha256"] == sha(args.reference)
    assert saved["registration_sha256"] == sha(source / "registration.json")
    assert saved["code_sha256"] == sha(source / "run.py")
    assert registration["dimension"] == 8 and registration["ordered_m"][0] == 64
    assert len(saved["trials"]) == 1 and saved["trials"][0]["m"] == 64
    expected = saved["trials"][0]
    reference = Q(old["full_construction_comparison_interval"][1])
    inherited = Q(old["finite_L1"]["exact_ratio_lower"])
    assert reference == Q(saved["comparison_upper"])
    assert inherited == Q(saved["inherited_finite_lower"])
    labels = list(product((1, 2), *[(1, 2, 3, 4)] * 7))
    assert len(labels) == 32768
    masses = sorted({prod(e) for e in labels})
    roots = {m: root8(m) for m in set(masses) | {1, 2, 64, 128}}
    assert {str(m): str(roots[m]) for m in masses} == saved["root_lower_numerators"]
    first_spacing = (2 * (roots[64] + roots[128]), 2 * (roots[64] + roots[128] + 2))
    tail_spacing = (2 * (D + roots[2]), 2 * (D + roots[2] + 1))
    spacings = [first_spacing] + [tail_spacing] * 7
    cell_upper = Q(prod(hi for lo, hi in spacings), S**8)
    assert cell_upper == Q(expected["cell_upper"]) and cell_upper > inherited
    sets = {"inner": [], "outer": []}
    positive = {"inner": 0, "outer": 0}
    widths = []
    checks = 0
    endpoint_checks = 0
    collar_checks = 0
    for lengths in labels:
        m = prod(lengths)
        radius_lo, radius_hi = 2 * roots[m], 2 * (roots[m] + 1)
        boxes = {"inner": [], "outer": []}
        source_ranges = []
        for ell, (alo, ahi) in zip(lengths, spacings):
            q = -(ell // 2) if ell % 2 else 1 - ell // 2
            p = q + ell - 1
            assert q <= 0 <= p
            source_ranges.append(range(q, p + 1))
            # Bound interval of max(0,p*a-r) and min(a/2,q*a+r).
            boxes["inner"].extend((max(0, p * ahi - radius_lo), min(alo // 2, q * ahi + radius_lo)))
            boxes["outer"].extend((max(0, p * alo - radius_hi), min(ahi // 2, q * alo + radius_hi)))
            for index in range(q, p + 1):
                for sign in (-1, 1):
                    assert -2 <= sign * index <= 2
                    collar_checks += 1
        for mode in ("inner", "outer"):
            box = tuple(boxes[mode])
            if all(box[2*i] < box[2*i+1] for i in range(8)):
                positive[mode] += 1
                sets[mode].append(box)
        box = boxes["inner"]
        if any(box[2*i] >= box[2*i+1] for i in range(8)):
            continue
        for axis, (alo, ahi) in enumerate(spacings):
            lo, hi = box[2*axis:2*axis+2]
            assert 0 <= lo < hi <= alo // 2
            widths.append(hi-lo)
            q, p = source_ranges[axis].start, source_ranges[axis].stop - 1
            assert lo >= p * ahi - radius_lo
            assert hi <= q * ahi + radius_lo
            endpoint_checks += 2
            # Every source coordinate lies within the SAME radius, not a
            # separately chosen radius for each coordinate or source point.
            for index in source_ranges[axis]:
                point_lo, point_hi = sorted((index * alo, index * ahi))
                assert hi - point_lo <= radius_lo
                assert point_hi - lo <= radius_lo
                checks += 2
    assert positive == {"inner": 882, "outer": 882}
    assert positive["inner"] == expected["positive_inward_raw_boxes"]
    assert positive["outer"] == expected["positive_outward_raw_boxes"]
    assert endpoint_checks == expected["containment_checks"]
    budget = [0]
    volumes = {}
    statistics = {}
    for mode in ("inner", "outer"):
        state = tuple(sorted(set(sets[mode])))
        engine = FixedAxisSweep(budget)
        t0 = time.time()
        volumes[mode] = Q(256 * engine.volume(state), S**8)
        statistics[mode] = {"positive_label_boxes": positive[mode], "distinct_physical_boxes": len(state), "memo_states": len(engine.memo), "seconds": time.time() - t0, "max_state_boxes": engine.max_state_boxes}
        print(json.dumps({"progress": mode, **statistics[mode]}), file=__import__('sys').stderr, flush=True)
    assert volumes["inner"] <= volumes["outer"] <= cell_upper
    assert [str(volumes[x]) for x in ("inner", "outer")] == expected["union_interval"]
    eps = Q(1, 10**25)
    eta = eps / 2
    theta = 1 - Q(1, 10**18)
    K = 10**19
    assert Q(min(widths), S) == Q(expected["minimum_width"])
    assert Q(min(widths), S) > 2 * eps
    assert all(Q(w, S) < 2 for w in widths)
    # A smooth normalized tensor bump supported in [-eta,eta]^8 exists;
    # shrinking each receiver face by eps provides >eta source margin.
    loss = 256 * positive["inner"] * 16 * eps * 2**7
    finite = theta * (volumes["inner"] - loss) * Q(2*K+1, 2*K+5)**8
    assert finite > reference
    finite_values = {"epsilon": str(eps), "eta": str(eta), "theta": str(theta), "K": K, "collar": 2, "loss_upper": str(loss), "lower": str(finite), "margin": str(finite-reference)}
    assert finite_values == expected["finite_L1"]
    assert expected["winner"] is True
    verify_inputs()
    output = {"status": "PASS", "method": "all 32768 labels, fixed-axis integer event sweep; no production imports or orbit compression", "candidate_m": 64, "source_manifest_sha256": frozen["candidate_files"]["manifest.json"], "reference_sha256": sha(args.reference), "checker_sha256": sha(__file__), "memo_limit_total": LIMIT, "memo_states_total": budget[0], "root_lower_numerators": {str(m): str(roots[m]) for m in sorted(roots)}, "statistics": statistics, "containment_endpoint_checks": endpoint_checks, "containment_all_source_coordinate_checks": checks, "reflection_collar_coordinate_checks": collar_checks, "minimum_width": str(Q(min(widths), S)), "cell_upper": str(cell_upper), "union_interval": [str(volumes[x]) for x in ("inner", "outer")], "finite_L1": finite_values, "comparison_upper": str(reference), "decimal_display": {"union_lower": decimal(volumes["inner"]), "union_upper": decimal(volumes["outer"]), "finite_lower": decimal(finite), "old_full_upper": decimal(reference), "margin": decimal(finite-reference)}, "wall_seconds": time.time()-started, "peak_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
