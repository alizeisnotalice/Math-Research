"""Exact finite local-source graph clock guards, final registered batch.

No receiver, maximal, Phi, Monte Carlo or old script imports.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import datetime
import gzip
import hashlib
import json
import time

BASE = Path(__file__).resolve().parent
PREFIX = "local_clock_graph_guard_20261007"
REG = BASE / (PREFIX + "_registration.json")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


def save(path, obj):
    with Path(path).open("x") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)
        handle.write("\n")


def labels(mask, N):
    return [j for j in range(N) if mask & (1 << j)]


def harmonic(m):
    return sum((F(1, k) for k in range(1, m + 1)), F(0))


def make_case(rd, n, raw, graph):
    N, d = len(raw), F(1, 8 * n)
    if graph == "clique":
        first = [j * d for j in range(N)]
    elif graph == "path":
        first = [F(3, 2) * j for j in range(N)]
    elif graph == "cluster":
        first = [F(4 * (j // 3)) + (j % 3) * d for j in range(N)]
    else:
        raise ValueError(graph)
    source = {"n": n, "N": N, "raw_weights": raw,
              "weights": list(map(str, [F(w, sum(raw)) for w in raw])), "W": "1",
              "atoms": [list(map(str, [value] + [F(0)] * (n - 1))) for value in first]}
    return {"case_id": f"r{rd}_{graph}_N{N}_n{n}", "round": rd, "graph_recipe": graph,
            "N": N, "n": n, "a": "1", "b": "2", "d": str(d),
            "source": source, "source_sha256": hashlib.sha256(canonical(source)).hexdigest()}


def register():
    cases = [make_case(rd, n, raw, graph)
             for rd, n, raw in [(1, 1, [1, 2, 4]),
                                (2, 16, [1, 2, 3, 5, 11]),
                                (3, 256, [1, 1, 2, 3, 5, 8, 13, 55])]
             for graph in ["clique", "path", "cluster"]]
    spine = BASE / "stopped_source_spine_guard_20261007_registration.json"
    save(REG, {
        "status": "registered_before_local_clock_graph_tables",
        "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": digest(__file__), "cases": cases,
        "same_weight_reference": {"path": str(spine), "sha256": digest(spine),
                                  "scope": "Same complete raw masses only; source positions newly registered, no spine/main rerun."},
        "edge": "Distinct complete half-open d cells; undirected conservative edge iff full coordinate l-infinity distance <= b+2d, no self loop. Graph generated from Fraction coordinates, not stipulated adjacency.",
        "active": "For alive A, a_j=1 iff j has another alive graph neighbor; otherwise0. Under each fixed immortal tag i in A, non-tag active j dies with rate1, inactive with rate0, tag rate0.",
        "recurrence": "r=sum_{active j != i}1. r0 terminal, requires a_i=0 and induced graph independent. H_i(A)=a_i/r+sum_{active j!=i} H_i(A-j)/r; terminal0.",
        "target": "H_i(A)=harmonic(deg_A(i)) for every A/tag, and all terminal A independent set. Work is integrated active indicator of tag, not elapsed time of every disconnected component.",
        "source_average": "Full initial source once: sum_i w_i H_i(full)/W; static initial source weights, not sampling source as receiver and not summing all subset budgets.",
        "budget": {"cases": 9, "nonempty_A_graph_rows": 879,
                   "tag_state_rows": 3348, "all_pair_distance_checks": 123},
        "seeds": "No random seeds; explicit positions and weights, all finite subsets/tags exhaustively integrated by recursion.",
        "stopping": "One complete pass through nine fixed graphs; retain isolated tags/non-tag components, no scope expansion or old oracle execution.",
        "evidence": "All coordinates, rates, recursion and averages Fraction; no log or time approximation.",
        "scope": "Finite legal local-clock complexity implementation guard only; not realmax/Phi, nearmax/geom qualification, receiver integration, or dimension-order/sqrt evidence.",
    })
    print(json.dumps({"registration_sha256": digest(REG)}))


def process(case):
    n, N, d, b = case["n"], case["N"], F(case["d"]), F(case["b"])
    atoms = [[F(value) for value in y] for y in case["source"]["atoms"]]
    weights = list(map(F, case["source"]["weights"]))
    assert all(w > 0 for w in weights) and sum(weights) == 1
    cells = [tuple(int(coord // d) for coord in y) for y in atoms]
    assert len(set(cells)) == N
    assert all(all(coord == d * index for coord, index in zip(y, cell)) for y, cell in zip(atoms, cells))
    adjacency = [0] * N
    pairs = []
    for i in range(N):
        for j in range(i + 1, N):
            distance = max(abs(yi - yj) for yi, yj in zip(atoms[i], atoms[j]))
            edge = distance <= b + 2 * d
            if edge:
                adjacency[i] |= 1 << j
                adjacency[j] |= 1 << i
            pairs.append({"i": i, "j": j, "l_infinity_distance": str(distance),
                          "edge_threshold": str(b + 2 * d), "edge": edge})
    for i in range(N):
        for j in range(N):
            expected = i != j if case["graph_recipe"] == "clique" else (abs(i - j) == 1 if case["graph_recipe"] == "path" else i != j and i // 3 == j // 3)
            assert bool(adjacency[i] & (1 << j)) == expected
    masks = sorted(range(1, 2 ** N), key=lambda A: (A.bit_count(), A))
    H, rows = {}, []
    terminal_tag_rows, isolated_nonterminal_tags, contribution_count = 0, 0, 0
    for A in masks:
        js = labels(A, N)
        degrees = {j: (adjacency[j] & A).bit_count() for j in js}
        active = {j: int(degrees[j] > 0) for j in js}
        independent = all(degrees[j] == 0 for j in js)
        wA = sum(weights[j] for j in js)
        tags = []
        for i in js:
            active_nontag = [j for j in js if j != i and active[j]]
            r = len(active_nontag)
            contributions = []
            if r == 0:
                assert active[i] == 0 and independent
                value = F(0)
                terminal_tag_rows += 1
            else:
                assert not independent
                value = F(active[i], r)
                isolated_nonterminal_tags += active[i] == 0
                for j in active_nontag:
                    child = A ^ (1 << j)
                    assert child & (1 << i)
                    contribution = H[(i, child)] / r
                    value += contribution
                    contributions.append({"deleted_j": j, "child_mask": child,
                                          "conditional_jump_probability": str(F(1, r)),
                                          "child_H": str(H[(i, child)]),
                                          "recurrence_contribution": str(contribution)})
                    contribution_count += 1
            expected = harmonic(degrees[i])
            assert value == expected
            H[(i, A)] = value
            tags.append({"tag": i, "tag_immortal_rate": "0", "alive_degree": degrees[i],
                         "active_tag": active[i], "active_nontag_labels": active_nontag,
                         "non_tag_rates": [{"j": j, "rate": active[j]} for j in js if j != i],
                         "total_rate": r, "terminal": r == 0,
                         "tag_first_hold_work": "0" if r == 0 else str(F(active[i], r)),
                         "H_tag_active_work": str(value), "harmonic_degree": str(expected),
                         "H_minus_harmonic": "0", "contributions": contributions})
        average = sum(weights[i] * H[(i, A)] for i in js) / wA
        rows.append({"A_mask": A, "A_labels": js, "wA": str(wA),
                     "active": [{"j": j, "a_j": active[j], "degree": degrees[j]} for j in js],
                     "independent_set": independent, "tags": tags,
                     "conditional_initial_source_weight_average": str(average)})
    full = 2 ** N - 1
    fullwork = sum(weights[i] * H[(i, full)] for i in range(N))
    summary = {"case_id": case["case_id"], "round": case["round"], "N": N, "n": n,
               "graph_recipe": case["graph_recipe"], "source_cells": N,
               "edges": sum(row["edge"] for row in pairs), "pair_distance_checks": len(pairs),
               "nonempty_A_rows": len(rows), "tag_state_rows": sum(len(row["tags"]) for row in rows),
               "harmonic_identity_passes": sum(len(row["tags"]) for row in rows),
               "independent_nonempty_states": sum(row["independent_set"] for row in rows),
               "terminal_tag_rows": terminal_tag_rows, "isolated_nonterminal_tag_rows": isolated_nonterminal_tags,
               "recurrence_child_contributions": contribution_count,
               "full_source_weighted_tag_work": str(fullwork),
               "full_tag_work": [str(H[(i, full)]) for i in range(N)], "failures": 0}
    return {"case": case, "source_cell_indices": [list(cell) for cell in cells],
            "pair_distance_graph": pairs, "neighbor_masks": adjacency,
            "states": rows, "summary": summary}


def run():
    start = time.monotonic()
    reg = json.loads(REG.read_text())
    assert reg["script_sha256"] == digest(__file__)
    ref = reg["same_weight_reference"]
    assert digest(ref["path"]) == ref["sha256"]
    spine_cases = json.loads(Path(ref["path"]).read_text())["cases"]
    manifests, summaries = [], []
    for case in reg["cases"]:
        assert hashlib.sha256(canonical(case["source"])).hexdigest() == case["source_sha256"]
        assert case["source"]["weights"] == next(c for c in spine_cases if c["N"] == case["N"])["source"]["weights"]
        data = process(case)
        path = BASE / (PREFIX + "_" + case["case_id"] + ".json.gz")
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                zipped.write(canonical(data))
        manifests.append({"path": str(path), "sha256": digest(path)})
        summaries.append(data["summary"])
        print(json.dumps({"done": case["case_id"], "tag_states": data["summary"]["tag_state_rows"]}), flush=True)
    resultpath = BASE / (PREFIX + "_results.json")
    save(resultpath, {"status": "complete_registered_local_graph_clock_guard", "summaries": summaries,
                      "profiles": manifests, "elapsed_seconds": time.monotonic() - start,
                      "registration_sha256": digest(REG), "script_sha256": digest(__file__),
                      "scope": reg["scope"], "old_oracle_calls": 0, "Phi_evaluations": 0})
    receipt = BASE / (PREFIX + "_receipt.json")
    save(receipt, {"status": "complete", "exit_status": 0, "main_runs": 1,
                   "registration": {"path": str(REG), "sha256": digest(REG)}, "script": {"path": __file__, "sha256": digest(__file__)},
                   "results": {"path": str(resultpath), "sha256": digest(resultpath)},
                   "profiles": manifests, "old_oracle_calls": 0, "Phi_evaluations": 0,
                   "scope": reg["scope"]})
    print(json.dumps({"results_sha256": digest(resultpath), "receipt_sha256": digest(receipt)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["register", "run"])
    mode = parser.parse_args().mode
    register() if mode == "register" else run()
