#!/usr/bin/env python3
"""Independent full-label, fixed-coordinate exact event sweep for n8 R003 left event.

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
    parser.add_argument("--theory-root", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    frozen = json.loads((HERE / "frozen_inputs.json").read_text())
    source = args.candidate_root
    def verify_inputs():
        for name, digest in frozen["candidate_files"].items():
            assert sha(source / name) == digest, "changed input " + name
        assert sha(args.reference) == frozen["reference_sha256"]
        for name, digest in frozen["theory_files"].items():
            assert sha(args.theory_root / name) == digest, "changed theory input " + name
    verify_inputs()
    manifest = json.loads((source / "manifest.json").read_text())
    assert manifest["files"] == {k: v for k, v in frozen["candidate_files"].items() if k != "manifest.json"}
    saved = json.loads((source / "candidate_left.json").read_text())
    registration = json.loads((source / "registration.json").read_text())
    old = json.loads(args.reference.read_text())
    lock = json.loads((source / "input_lock.json").read_text())
    assert saved["input_lock_sha256"] == sha(source / "input_lock.json")
    assert lock["r002_candidate_96_sha256"] == sha(args.reference)
    for name, digest in lock["local_sha256"].items():
        assert sha(source / name) == digest
    assert sha(source / "theory_manifest.json") == sha(args.theory_root / "manifest.json")
    assert sha(source / "theory_event_certificate.json") == sha(args.theory_root / "event_certificate.json")
    assert registration["dimension"] == 8 and registration["events"] == [{"name":"left","roots":[64,256]},{"name":"right","roots":[32,486]}]
    assert saved["event"] == {"name":"left","roots":[64,256]} and saved["status"] == "complete"
    expected = saved["result"]
    assert expected["event"] == "left"
    reference = Q(old["candidate"]["union_interval"][1])
    inherited = Q(old["candidate"]["finite_L1"]["lower"])
    assert reference == Q(saved["comparison_upper"])
    assert inherited == Q(saved["inherited_finite_lower"])
    theory = json.loads((args.theory_root / "event_certificate.json").read_text())
    DT = int(theory["integer_root_scale_D"])
    rt = {int(m): int(q) for m, q in theory["root_floor_numerators"].items()}
    for m, q in rt.items():
        assert q**8 <= m * DT**8 < (q+1)**8
        assert q == isqrt(isqrt(isqrt(m * DT**8)))
    left_iv = (Q(rt[64]+rt[256],2*DT),Q(rt[64]+rt[256]+2,2*DT))
    middle_iv = (Q(rt[96]+rt[192],2*DT),Q(rt[96]+rt[192]+2,2*DT))
    right_iv = (Q(rt[32]+rt[486],2*DT),Q(rt[32]+rt[486]+2,2*DT))
    assert left_iv[1] < middle_iv[0] < middle_iv[1] < right_iv[0]
    assert theory["candidates"][0]["aliases"] == [{"kind":"sum","p":32,"q":486}]
    assert theory["candidates"][2]["aliases"] == [{"kind":"sum","p":64,"q":256}]
    labels = list(product((1, 2), *[(1, 2, 3, 4)] * 7))
    assert len(labels) == 32768
    masses = sorted({prod(e) for e in labels})
    roots = {m: root8(m) for m in set(masses) | {1, 2, 64, 256}}
    first_spacing = (2 * (roots[64] + roots[256]), 2 * (roots[64] + roots[256] + 2))
    assert [str(roots[64]), str(roots[256])] == saved["event_root_lower_numerators"]
    assert [str(z) for z in first_spacing] == expected["spacing_endpoints"]
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
    tail_counts = {1:0,2:0}
    tail_masses = {1:set(),2:set()}
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
        tail_in = all(boxes["inner"][2*i] < boxes["inner"][2*i+1] for i in range(1,8))
        tail_out = all(boxes["outer"][2*i] < boxes["outer"][2*i+1] for i in range(1,8))
        assert tail_in == tail_out
        if tail_in:
            tail_counts[lengths[0]] += 1
            tail_masses[lengths[0]].add(m)
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
    assert {str(k):v for k,v in tail_counts.items()} == theory["positive_tail_raw_labels"]
    assert {str(k):sorted(v) for k,v in tail_masses.items()} == theory["positive_tail_mass_sets"]
    assert positive == {"inner": 875, "outer": 875}
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
    assert volumes["inner"] > reference
    right = json.loads((source / "candidate_right.json").read_text())
    assert right["event"] == {"name":"right","roots":[32,486]} and right["status"] == "complete"
    assert right["input_lock_sha256"] == sha(source / "input_lock.json")
    assert Q(right["comparison_upper"]) == reference
    right_lo, right_hi = map(Q, right["result"]["union_interval"])
    assert right_lo <= right_hi < reference
    saved_comparison = {"event":"right", "status":"saved rational comparison only; union NOT independently recomputed", "saved_lower":str(right_lo), "saved_upper":str(right_hi), "reference_minus_saved_upper":str(reference-right_hi)}
    for name in ("left","right"):
        original_receipt = json.loads((source / f"receipt_{name}.json").read_text())
        assert original_receipt["returncode"] == 0 and original_receipt["name"] == name
    theory_checks = {"source_manifest_sha256":sha(args.theory_root / "manifest.json"), "source_proof_sha256":sha(args.theory_root / "event_theory.md"), "root_isolations_verified":len(rt), "positive_tail_raw_labels":tail_counts, "positive_tail_mass_sets":{k:sorted(v) for k,v in tail_masses.items()}, "left_middle_right_strict_bracketing":True, "342_event_list_regenerated":False, "all_event_union_optimization_performed":False}
    verify_inputs()
    output = {"status": "PASS", "method": "all 32768 labels, fixed-axis integer event sweep; no production imports or orbit compression", "candidate_event": {"name":"left","roots":[64,256]}, "right_saved_comparison": saved_comparison, "theory_checks":theory_checks, "author_original_process_exit_codes": {"left":0,"right":0}, "source_manifest_sha256": frozen["candidate_files"]["manifest.json"], "reference_sha256": sha(args.reference), "checker_sha256": sha(__file__), "memo_limit_total": LIMIT, "memo_states_total": budget[0], "root_lower_numerators": {str(m): str(roots[m]) for m in sorted(roots)}, "statistics": statistics, "containment_endpoint_checks": endpoint_checks, "containment_all_source_coordinate_checks": checks, "reflection_collar_coordinate_checks": collar_checks, "minimum_width": str(Q(min(widths), S)), "cell_upper": str(cell_upper), "union_interval": [str(volumes[x]) for x in ("inner", "outer")], "finite_L1": finite_values, "comparison_upper": str(reference), "decimal_display": {"union_lower": decimal(volumes["inner"]), "union_upper": decimal(volumes["outer"]), "finite_lower": decimal(finite), "old_full_upper": decimal(reference), "margin": decimal(finite-reference)}, "wall_seconds": time.time()-started, "peak_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
