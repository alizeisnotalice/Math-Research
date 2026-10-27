#!/usr/bin/env python3
"""Independent full-label, fixed-coordinate exact event sweep for n8 R005 E098 event.

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

def audit_catalog(args, registration, positive_masses, positive_counts, winner_lower, saved_intervals):
    """Independent event arithmetic/provenance; never evaluates an old union."""
    source = args.candidate_root
    theory = args.theory_root
    catalog = json.loads((source / "theory_catalog.json").read_text())
    snapshot = json.loads((theory / "input_snapshot.json").read_text())
    old_certificate = json.loads((theory / "frozen_r003_event_certificate.json").read_text())
    scale = int(old_certificate["integer_root_scale_D"])
    rootq = {}
    for sm, sq in old_certificate["root_floor_numerators"].items():
        m, q = int(sm), int(sq)
        assert q == isqrt(isqrt(isqrt(m * scale**8)))
        assert q**8 <= m*scale**8 < (q+1)**8
        rootq[m] = q
    assert {str(k):sorted(v) for k,v in positive_masses.items()} == old_certificate["positive_tail_mass_sets"]
    assert {str(k):v for k,v in positive_counts.items()} == old_certificate["positive_tail_raw_labels"]

    def terms(alias):
        if alias["kind"] == "sum":
            return [(alias["p"],Q(1,2)),(alias["q"],Q(1,2))]
        return [(alias["p"],Q(1,2) if alias["kind"] == "r" else Q(1))]

    def sig(term_list):
        accum = defaultdict(Q)
        for radicand, coefficient in term_list:
            powers = []
            remaining = radicand
            for prime in (2,3):
                power = 0
                while remaining % prime == 0:
                    power += 1
                    remaining //= prime
                powers.append(power)
            assert remaining == 1
            a,b = powers
            accum[(a%8,b%8)] += coefficient*2**(a//8)*3**(b//8)
        return tuple((a,b,c.numerator,c.denominator) for (a,b),c in sorted(accum.items()))

    def interval(term_list):
        lo = hi = Q(0)
        for m,c in term_list:
            q = rootq[m]
            top = q if q**8 == m*scale**8 else q+1
            lo += c*Q(q,scale)
            hi += c*Q(top,scale)
        return lo,hi

    # Rebuild every permitted endpoint crossing from independently checked
    # positive-tail masses; canonical equality plus strict intervals suffices.
    aliases = [{"kind":"r","p":m,"q":None} for m in sorted(positive_masses[2])]
    aliases += [{"kind":"2r","p":m,"q":None} for m in sorted(positive_masses[1] | positive_masses[2])]
    aliases += [{"kind":"sum","p":p,"q":q} for p in sorted(positive_masses[1]) for q in sorted(positive_masses[2]) if p<q]
    grouped = {}
    for alias in aliases:
        key = sig(terms(alias))
        lo,hi = interval(terms(alias))
        if key not in grouped:
            grouped[key] = {"lo":lo,"hi":hi,"aliases":[]}
        row = grouped[key]
        row["lo"] = max(row["lo"],lo)
        row["hi"] = min(row["hi"],hi)
        row["aliases"].append(alias)
        assert row["lo"] <= row["hi"]
    ordered = sorted(grouped.items(), key=lambda item:item[1]["lo"])
    assert len(aliases) == 345 and len(ordered) == 342
    for (_,a),(_,b) in zip(ordered,ordered[1:]):
        assert a["hi"] < b["lo"]
    assert len(catalog["events"]) == 342
    by_id = {}
    for i, ((key,row),e) in enumerate(zip(ordered,catalog["events"]),1):
        assert e["event_id"] == f"E{i:03d}"
        assert tuple(map(tuple,e["signature"])) == key
        assert tuple(map(Q,e["spacing_interval"])) == (row["lo"],row["hi"])
        assert {json.dumps(a,sort_keys=True) for a in e["aliases"]} == {json.dumps(a,sort_keys=True) for a in row["aliases"]}
        by_id[e["event_id"]] = e
    assert catalog["terminal_saturation_event_id"] == "E342"
    assert {"kind":"2r","p":32768,"q":None} in by_id["E342"]["aliases"]

    # Verify selected historical data against unchanged whole source files,
    # rather than treating an unbound JSON snapshot as new evidence.
    for key, record in snapshot["records"].items():
        path = source.parent / record["source_key"]
        assert sha(path) == record["source_sha256"]
        original = json.loads(path.read_text())
        if key == "m16":
            continue
        if key == "m64":
            assert record["data"]["result"] == original["trials"][0]
        elif key.startswith("m"):
            assert record["data"]["result"] == original["candidate"]
            assert record["data"]["m"] == original["m"]
            assert record["data"]["status"] == original["status"]
        else:
            assert set(record["data"]) == {"event","status","result"}
            assert all(value == original[field] for field,value in record["data"].items())

    old16_by_sig = {tuple(map(tuple,e["signature"])):e["rank"] for e in old_certificate["candidates"]}
    assert len(old16_by_sig) == 16
    actual_old16 = {}
    history_matches = {}
    for e in catalog["events"]:
        key = tuple(map(tuple,e["signature"]))
        rank = old16_by_sig.get(key)
        assert e["old16_distance_rank"] == rank
        assert e["outside_old16"] == (rank is None)
        if rank is not None:
            actual_old16[rank] = e["event_id"]
    assert set(actual_old16) == set(range(1,17))
    for sm,hist in catalog["historical_constructs"].items():
        m = int(sm)
        term_list = [(m,Q(1,2)),(2*m,Q(1,2))]
        key = sig(term_list)
        lo,hi = interval(term_list)
        assert tuple(map(tuple,hist["signature"])) == key
        assert tuple(map(Q,hist["spacing_interval"])) == (lo,hi)
        matches = [e["event_id"] for e in catalog["events"] if tuple(map(tuple,e["signature"])) == key]
        assert len(matches) <= 1 and hist["is_safe_event"] == bool(matches)
        assert hist["matching_event_id"] == (matches[0] if matches else None)
        if matches:
            history_matches[m] = matches[0]
        else:
            assert all(hi<Q(e["spacing_interval"][0]) or Q(e["spacing_interval"][1])<lo for e in catalog["events"])
    assert history_matches == {16:"E044",64:"E088",128:"E123",32:"E059"}
    old_ids = set(actual_old16.values())
    historic_ids = set(history_matches.values())
    assert not old_ids & historic_ids
    all_ids = set(by_id)
    new_ids = all_ids - old_ids - historic_ids
    assert len(old_ids) == 16 and len(historic_ids) == 4 and len(new_ids) == 322
    for e in catalog["events"]:
        eid = e["event_id"]
        matches = [m for m,x in history_matches.items() if x == eid]
        assert e["historical_m_matches"] == matches
        assert e["new_unique_event"] == (eid in new_ids)
        for binding in e["source_bindings"]:
            record = snapshot["records"][binding["snapshot_record"]]
            assert binding["source_key"] == record["source_key"] and binding["source_sha256"] == record["source_sha256"]
        if e["old16_distance_rank"] is not None:
            rank = e["old16_distance_rank"]
            record = snapshot["records"][f"rank{rank}"]["data"]
            source_key = sig([(m,Q(1,2)) for m in record["event"]["roots"]]) if rank in (1,3) else tuple(map(tuple,record["event"]["signature"]))
            assert source_key == tuple(map(tuple,e["signature"]))

    r2 = source.parent / "round002/candidate_96.json"
    screen_threshold = Q(json.loads(r2.read_text())["candidate"]["finite_L1"]["lower"])
    screen_info = catalog["cell_prescreen_reference"]
    assert screen_info["binding"]["source_sha256"] == sha(r2)
    assert Q(screen_info["finite_L1_lower"]) == screen_threshold
    audit2 = HERE.parents[1] / "round002/lower"
    audit2_binding = snapshot["audit_r002_m96_binding"]
    assert sha(audit2 / "manifest.json") == audit2_binding["manifest_sha256"]
    for name,digest in audit2_binding["files"].items():
        assert sha(audit2 / name) == digest
    assert Q(json.loads((audit2 / "exactresults.json").read_text())["finite_L1"]["lower"]) == screen_threshold
    tail_upper = Q(scale + rootq[2] + 1, 2*scale)
    screened = set()
    all_cell_upper = {}
    for e in catalog["events"]:
        cell = Q(e["spacing_interval"][1]) * tail_upper**7
        assert cell == Q(e["cell_volume_upper"])
        assert e["cell_upper_below_audited_r002_lower"] == (cell < screen_threshold)
        all_cell_upper[e["event_id"]] = cell
        if e["event_id"] in new_ids and cell < screen_threshold:
            screened.add(e["event_id"])
    assert len(screened) == 56 and all_cell_upper["E044"] < screen_threshold
    assert set(catalog["new_exact_cell_prescreen_event_ids"]) == screened == set(registration["prescreened_event_ids"])
    needed = new_ids - screened
    assert len(needed) == 266 and set(catalog["new_union_needed_after_cell_prescreen_event_ids"]) == needed
    assert set(catalog["new_unique_event_ids"]) == new_ids
    assert set(catalog["outside_old16_event_ids"]) == all_ids-old_ids
    assert {e["event_id"] for e in registration["tasks"]} == needed
    center = tuple(map(Q,by_id["E101"]["spacing_interval"]))
    previous_dist_hi = None
    for task in registration["tasks"]:
        e = by_id[task["event_id"]]
        for field in ("signature","aliases","spacing_interval"):
            assert task[field] == e[field]
        assert task["selected_alias"] in e["aliases"]
        lo,hi = map(Q,task["spacing_interval"])
        if hi < center[0]:
            dist = center[0]-hi,center[1]-lo
        else:
            assert center[1] < lo
            dist = lo-center[1],hi-center[0]
        assert tuple(map(Q,task["distance_to_R004_winner_interval"])) == dist
        if previous_dist_hi is not None:
            assert previous_dist_hi < dist[0]
        previous_dist_hi = dist[1]
    assert set(json.loads((source / "task_state.json").read_text())["excluded"]) == screened

    # All 342 are mapped, but an accurate comparison of saved upper fractions
    # is not an independent recomputation of those unions.
    coverage = []
    for e in catalog["events"]:
        eid = e["event_id"]
        if eid in screened or eid == "E044":
            upper = all_cell_upper[eid]
            evidence = "independently verified strict cell screen"
        elif eid in saved_intervals:
            upper = saved_intervals[eid][1]
            evidence = "independent full union" if eid == "E098" else "saved full union only; not independently recomputed"
        else:
            binding = e["source_bindings"][0]
            rec = snapshot["records"][binding["snapshot_record"]]["data"]
            result = rec["result"]
            upper = Q(result["union_interval"][1])
            evidence = "earlier independently audited winner" if eid in ("E088","E101","E110") else "earlier saved full union only; not independently recomputed"
        coverage.append({"event_id":eid,"evidence":evidence,"winner_lower_above_upper":winner_lower>upper if eid != "E098" else None})
    for m in (32,128):
        eid = history_matches[m]
        assert "wrapper failed" in by_id[eid]["evidence_note"]
        original_text = (source.parent / f"round002/receipt_{m}.txt").read_text()
        assert "TypeError" in original_text and "flush" in original_text
    return {"raw_expressions_independently_regenerated":345,"distinct_events_independently_regenerated":342,"all_catalog_intervals_strictly_ordered":True,"old16_count":16,"additional_history_count":4,"new_unique_count":322,"new_strict_cell_screens":56,"historical_m16_fresh_cell_screen":True,"new_union_tasks":266,"unknown":0,"history_event_ids":history_matches,"screen_threshold":str(screen_threshold),"largest_new_screen_cell_upper":str(max(all_cell_upper[e] for e in screened)),"source_snapshot_records_verified":len(snapshot["records"]),"scope":"catalog geometry, identities, mapping, cell screens and saved comparisons; only E098 new union independently recomputed", "evidence_inventory":coverage,"unseparated_or_higher_saved_event_upper_ids":[e["event_id"] for e in coverage if e["winner_lower_above_upper"] is False],"independent_full_family_optimum_certified":False,"original_m32_m128_status":"saved complete JSON; original wrapper TypeError after save; actual original exit code unknown; no new independent union", "global_C8_optimum_certified":False}

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
        for item in frozen["external_readonly_sources"].values():
            assert sha(item["path"]) == item["sha256"], "changed external source"
    verify_inputs()
    manifest = json.loads((source / "manifest.json").read_text())
    assert manifest["files"] == {k: v for k, v in frozen["candidate_files"].items() if k != "manifest.json"}
    saved = json.loads((source / "candidate_E098.json").read_text())
    registration = json.loads((source / "registration.json").read_text())
    old = json.loads(args.reference.read_text())
    lock = json.loads((source / "input_lock.json").read_text())
    assert saved["input_lock_sha256"] == sha(source / "input_lock.json")
    assert registration["reference_candidate_sha256"] == sha(args.reference)
    for name, digest in lock["files"].items():
        assert sha(source / name) == digest
    assert sha(source / "theory_catalog.json") == sha(args.theory_root / "catalog.json")
    assert sha(source / "theory_manifest.json") == sha(args.theory_root / "manifest.json")
    assert len(registration["tasks"]) == 266
    selected = next(e for e in registration["tasks"] if e["event_id"] == "E098")
    assert selected["selected_alias"] == {"kind":"sum","p":1,"q":1944}
    assert saved["event"] == selected and saved["status"] == "complete"
    expected = saved["result"]
    assert expected["event"] == "E098"
    reference = Q(old["result"]["union_interval"][1])
    inherited = Q(old["result"]["finite_L1"]["lower"])
    assert reference == Q(saved["comparison_upper"])
    assert inherited == Q(saved["inherited_finite_lower"])
    labels = list(product((1, 2), *[(1, 2, 3, 4)] * 7))
    assert len(labels) == 32768
    masses = sorted({prod(e) for e in labels})
    roots = {m: root8(m) for m in set(masses) | {1, 2, 1944}}
    first_spacing = (2 * (roots[1] + roots[1944]), 2 * (roots[1] + roots[1944] + 1))
    assert len(saved["event_roots"]) == 2
    for term, m in zip(saved["event_roots"], (1, 1944)):
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
        ti = all(boxes["inner"][2*i] < boxes["inner"][2*i+1] for i in range(1,8))
        to = all(boxes["outer"][2*i] < boxes["outer"][2*i+1] for i in range(1,8))
        assert ti == to
        if ti:
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
    state_file = json.loads((source / "task_state.json").read_text())
    assert state_file["status"] == "finished" and state_file["completed"] == 266 and state_file["unknown"] == 0
    assert state_file["registration_sha256"] == sha(source / "registration.json")
    assert set(state_file["tasks"]) == {e["event_id"] for e in registration["tasks"]}
    for e in registration["tasks"]:
        eid = e["event_id"]
        other = json.loads((source / f"candidate_{eid}.json").read_text())
        receipt = json.loads((source / f"receipt_{eid}.json").read_text())
        st = state_file["tasks"][eid]
        assert st["status"] == "complete" and st["returncode"] == 0
        assert st["candidate"] == f"candidate_{eid}.json" and st["receipt"] == f"receipt_{eid}.json"
        assert other["event"] == e and other["status"] == "complete"
        assert other["input_lock_sha256"] == sha(source / "input_lock.json")
        assert Q(other["comparison_upper"]) == reference and Q(other["inherited_finite_lower"]) == inherited
        assert receipt["event_id"] == eid and receipt["returncode"] == 0 and receipt["stderr"] == ""
        original_codes[eid] = receipt["returncode"]
        c = other["result"]
        lo, hi = map(Q, c["union_interval"])
        assert lo <= hi <= Q(c["cell_upper"])
        assert c["states"] == c["inner_states"] + c["outer_states"] <= registration["cache_states_per_candidate"]
        fl = c["finite_L1"]
        eps_o, theta_o, K_o = Q(fl["epsilon"]), Q(fl["theta"]), fl["K"]
        loss_o = 256 * c["positive_inward_raw_boxes"] * 16 * eps_o * 2**7
        value = theta_o * (lo-loss_o) * Q(2*K_o+1,2*K_o+5)**8
        assert loss_o == Q(fl["loss_upper"]) and value == Q(fl["lower"])
        assert value-reference == Q(fl["margin"]) and c["winner"] == (value > reference)
        saved_rows.append({"event_id":eid,"value":value,"lower":lo,"upper":hi,"source_sha256":sha(source / f"candidate_{eid}.json")})
        saved_intervals[eid] = (lo,hi)
    saved_rows.sort(key=lambda row:(-row["value"],row["event_id"]))
    assert saved_rows[0]["event_id"] == "E098"
    assert saved_rows[0]["value"] > saved_rows[1]["value"]
    assert state_file["best_event_id"] == "E098"
    author_order = json.loads((source / "verify_saved_receipt.json").read_text())
    assert [r["event_id"] for r in saved_rows] == [r["event_id"] for r in author_order["ordered"]]
    for a,b in zip(saved_rows,author_order["ordered"]):
        assert a["value"] == Q(b["finite_lower"])
    overlap_pairs = []
    separated_pairs = 0
    ids = sorted(saved_intervals)
    for i, a in enumerate(ids):
        for b in ids[i+1:]:
            al,ah = saved_intervals[a]; bl,bh = saved_intervals[b]
            if ah < bl or bh < al: separated_pairs += 1
            else: overlap_pairs.append([a,b])
    saved_ranking = {"scope":"only E098 independent full union;other265 saved arithmetic/ranking only", "ordered_saved_finite_ids":[r["event_id"] for r in saved_rows], "exact_finite_tie_groups":[[a["event_id"],b["event_id"]] for a,b in zip(saved_rows,saved_rows[1:]) if a["value"] == b["value"]], "top10_saved_values":[{"event_id":r["event_id"],"finite_lower":str(r["value"]),"union_interval":[str(r["lower"]),str(r["upper"])]} for r in saved_rows[:10]], "saved_pairwise_strictly_separated":separated_pairs,"saved_overlap_pairs":overlap_pairs,"independent_winner_inner_not_above_saved_upper":[eid for eid,(lo,hi) in saved_intervals.items() if eid != "E098" and not volumes["inner"] > hi],"saved_finite_above_R004_upper":[r["event_id"] for r in saved_rows if r["value"] > reference],"saved_outer_not_below_R004_upper":[r["event_id"] for r in saved_rows if not r["upper"] < reference]}
    catalog_checks = audit_catalog(args,registration,tail_masses,tail_counts,volumes["inner"],saved_intervals)
    verify_inputs()
    output = {"status": "PASS", "method": "all 32768 labels, fixed-axis integer event sweep; no production imports or orbit compression", "candidate_event": {"event_id":"E098","roots":[1,1944]}, "catalog_checks":catalog_checks,"saved_ranking": saved_ranking, "author_original_process_exit_codes": original_codes, "source_manifest_sha256": frozen["candidate_files"]["manifest.json"], "reference_sha256": sha(args.reference), "checker_sha256": sha(__file__), "memo_limit_total": LIMIT, "memo_states_total": budget[0], "root_lower_numerators": {str(m): str(roots[m]) for m in sorted(roots)}, "statistics": statistics, "containment_endpoint_checks": endpoint_checks, "containment_all_source_coordinate_checks": checks, "reflection_collar_coordinate_checks": collar_checks, "minimum_width": str(Q(min(widths), S)), "cell_upper": str(cell_upper), "union_interval": [str(volumes[x]) for x in ("inner", "outer")], "finite_L1": finite_values, "comparison_upper": str(reference), "decimal_display": {"union_lower": decimal(volumes["inner"]), "union_upper": decimal(volumes["outer"]), "finite_lower": decimal(finite), "old_full_upper": decimal(reference), "margin": decimal(finite-reference)}, "wall_seconds": time.time()-started, "peak_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
