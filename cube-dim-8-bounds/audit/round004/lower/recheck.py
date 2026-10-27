#!/usr/bin/env python3
"""Independent full-label, fixed-coordinate exact event sweep for n8 R004 rank16 event.

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
        for item in frozen["external_readonly_sources"].values():
            assert sha(item["path"]) == item["sha256"], "changed external source"
    verify_inputs()
    manifest = json.loads((source / "manifest.json").read_text())
    assert manifest["files"] == {k: v for k, v in frozen["candidate_files"].items() if k != "manifest.json"}
    saved = json.loads((source / "candidate_16.json").read_text())
    registration = json.loads((source / "registration.json").read_text())
    old = json.loads(args.reference.read_text())
    lock = json.loads((source / "input_lock.json").read_text())
    assert saved["input_lock_sha256"] == sha(source / "input_lock.json")
    assert registration["baseline_r003_left_sha256"] == sha(args.reference)
    for name, digest in lock["files"].items():
        assert sha(source / name) == digest
    ranks = [e["rank"] for e in registration["selected"]]
    assert ranks == [2,4,5,6,7,8,9,10,11,12,13,14,15,16]
    event16 = next(e for e in registration["selected"] if e["rank"] == 16)
    assert event16["selected_alias"] == {"kind":"sum","p":16,"q":576}
    assert saved["event"] == event16 and saved["status"] == "complete"
    expected = saved["result"]
    assert expected["event"] == "16"
    reference = Q(old["result"]["union_interval"][1])
    inherited = Q(old["result"]["finite_L1"]["lower"])
    assert reference == Q(saved["comparison_upper"])
    assert inherited == Q(saved["inherited_finite_lower"])
    labels = list(product((1, 2), *[(1, 2, 3, 4)] * 7))
    assert len(labels) == 32768
    masses = sorted({prod(e) for e in labels})
    roots = {m: root8(m) for m in set(masses) | {1, 2, 16, 576}}
    first_spacing = (2 * (roots[16] + roots[576]), 2 * (roots[16] + roots[576] + 2))
    assert len(saved["event_roots"]) == 2
    for term, m in zip(saved["event_roots"], (16, 576)):
        assert term["radicand"] == m and Q(term["coefficient"]) == Q(1,2)
        q = roots[m]
        exact = q**8 == m * D**8
        assert int(term["lower_numerator"]) == q
        assert int(term["upper_numerator"]) == q + (not exact)
        assert term["perfect_eighth_power"] == exact
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
    assert volumes["inner"] > reference
    saved_rows = []
    saved_intervals = {}
    original_codes = {}
    for e in registration["selected"]:
        rank = e["rank"]
        other = json.loads((source / f"candidate_{rank:02}.json").read_text())
        receipt = json.loads((source / f"receipt_{rank:02}.json").read_text())
        assert other["event"] == e and other["status"] == "complete"
        assert other["input_lock_sha256"] == sha(source / "input_lock.json")
        assert Q(other["comparison_upper"]) == reference
        assert Q(other["inherited_finite_lower"]) == inherited
        assert receipt["rank"] == rank and receipt["returncode"] == 0 and receipt["stderr"] == ""
        original_codes[str(rank)] = receipt["returncode"]
        c = other["result"]
        lo, hi = map(Q, c["union_interval"])
        assert lo <= hi <= Q(c["cell_upper"])
        assert c["states"] == c["inner_states"] + c["outer_states"] <= registration["memo_limit_per_candidate"]
        fl = c["finite_L1"]
        eps_o, theta_o, K_o = Q(fl["epsilon"]), Q(fl["theta"]), fl["K"]
        loss_o = 256 * c["positive_inward_raw_boxes"] * 16 * eps_o * 2**7
        value = theta_o * (lo-loss_o) * Q(2*K_o+1,2*K_o+5)**8
        assert loss_o == Q(fl["loss_upper"]) and value == Q(fl["lower"])
        assert value-reference == Q(fl["margin"]) and c["winner"] == (value > reference)
        saved_rows.append({"rank":rank,"saved_finite_lower":str(value),"saved_finite_minus_R003_upper":str(value-reference),"saved_inner":str(lo),"saved_outer":str(hi),"status":"winner independently recomputed" if rank == 16 else "saved arithmetic/ranking only; union NOT independently recomputed"})
        saved_intervals[rank] = (lo,hi)
    saved_rows.sort(key=lambda row:Q(row["saved_finite_lower"]),reverse=True)
    assert saved_rows[0]["rank"] == 16
    assert all(Q(saved_rows[i]["saved_finite_lower"]) > Q(saved_rows[i+1]["saved_finite_lower"]) for i in range(13))
    author_order = json.loads((source / "verify_saved_receipt.json").read_text())
    assert [r["rank"] for r in saved_rows] == [r["rank"] for r in author_order["ordered"]]
    pairwise = []
    for i, a in enumerate(ranks):
        for b in ranks[i+1:]:
            al, ah = saved_intervals[a]
            bl, bh = saved_intervals[b]
            if ah < bl:
                item = {"ranks":[a,b],"relation":"first strictly below second"}
            elif bh < al:
                item = {"ranks":[a,b],"relation":"second strictly below first"}
            else:
                item = {"ranks":[a,b],"relation":"overlap; complete-union order NOT resolved"}
            pairwise.append(item)
    winner_interval_comparisons = []
    for rank in ranks:
        if rank == 16: continue
        other_upper = saved_intervals[rank][1]
        winner_interval_comparisons.append({"other_rank":rank,"relation":"winner independently recomputed lower > other saved upper" if volumes["inner"] > other_upper else "overlap; complete-union winner order NOT resolved", "winner_lower_minus_saved_upper":str(volumes["inner"]-other_upper)})
    saved_ranking = {"scope":"only rank16 full geometry/union independently recomputed; other13 saved exact arithmetic/interval comparisons", "ordered_saved_finite_lowers":saved_rows,"pairwise_saved_union_interval_relations":pairwise,"winner_vs_other_saved_intervals":winner_interval_comparisons,"saved_interval_overlap_pairs":[x["ranks"] for x in pairwise if x["relation"].startswith("overlap")]}
    verify_inputs()
    output = {"status": "PASS", "method": "all 32768 labels, fixed-axis integer event sweep; no production imports or orbit compression", "candidate_event": {"rank":16,"roots":[16,576]}, "saved_ranking": saved_ranking, "author_original_process_exit_codes": original_codes, "source_manifest_sha256": frozen["candidate_files"]["manifest.json"], "reference_sha256": sha(args.reference), "checker_sha256": sha(__file__), "memo_limit_total": LIMIT, "memo_states_total": budget[0], "root_lower_numerators": {str(m): str(roots[m]) for m in sorted(roots)}, "statistics": statistics, "containment_endpoint_checks": endpoint_checks, "containment_all_source_coordinate_checks": checks, "reflection_collar_coordinate_checks": collar_checks, "minimum_width": str(Q(min(widths), S)), "cell_upper": str(cell_upper), "union_interval": [str(volumes[x]) for x in ("inner", "outer")], "finite_L1": finite_values, "comparison_upper": str(reference), "decimal_display": {"union_lower": decimal(volumes["inner"]), "union_upper": decimal(volumes["outer"]), "finite_lower": decimal(finite), "old_full_upper": decimal(reference), "margin": decimal(finite-reference)}, "wall_seconds": time.time()-started, "peak_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
