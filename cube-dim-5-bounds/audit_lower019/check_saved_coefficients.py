#!/usr/bin/env python3
"""One bounded independent arithmetic audit, plus one polynomial toy family.

Does not import or execute the author's field, moving-box union, or C union.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def integer_fifth_root(n):
    a, b = 0, 1
    while b**5 <= n:
        b *= 2
    while a+1 < b:
        c = (a+b)//2
        if c**5 <= n:
            a = c
        else:
            b = c
    assert a**5 <= n < (a+1)**5
    return a


def coefficients(spec):
    result = [Q(0)]*25
    for key, value in spec.items():
        i, j = map(int, key.split(","))
        assert 0 <= i < 5 and 0 <= j < 5
        result[5*i+j] = Q(value)
    return tuple(result)


def enclosure(cs, roots):
    total_low = total_high = Q(0)
    (ul, uh), (vl, vh) = roots
    for k, c in enumerate(cs):
        i, j = divmod(k, 5)
        low, high = ul**i*vl**j, uh**i*vh**j
        total_low += c*(low if c >= 0 else high)
        total_high += c*(high if c >= 0 else low)
    return total_low, total_high


def trim(p):
    p = list(map(Q, p))
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def plus(a, b):
    return trim([(a[k] if k < len(a) else 0)+(b[k] if k < len(b) else 0)
                 for k in range(max(len(a), len(b)))])


def minus(a, b):
    return plus(a, tuple(-x for x in b))


def times(a, b):
    if not a or not b:
        return ()
    p = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p[i+j] += x*y
    return trim(p)


def compare_right(a, b):
    """Polynomial order in a sufficiently small right neighborhood of zero."""
    for coefficient in minus(a, b):
        if coefficient:
            return 1 if coefficient > 0 else -1
    return 0


def pmin(polys):
    result = polys[0]
    for p in polys[1:]:
        if compare_right(p, result) < 0:
            result = p
    return result


def pmax(polys):
    result = polys[0]
    for p in polys[1:]:
        if compare_right(p, result) > 0:
            result = p
    return result


def union_polynomial(boxes):
    """Exact inclusion-exclusion of four polynomial rectangles, not jet slicing."""
    total = ()
    for size in range(1, len(boxes)+1):
        for selection in combinations(boxes, size):
            volume = (Q(1),)
            for axis in range(2):
                lower = pmax([box[axis][0] for box in selection])
                upper = pmin([box[axis][1] for box in selection])
                length = minus(upper, lower)
                if compare_right(length, ()) <= 0:
                    volume = ()
                    break
                volume = times(volume, length)
            total = plus(total, volume if size % 2 else tuple(-x for x in volume))
    return total


def evaluate(p, theta):
    return sum((c*theta**i for i, c in enumerate(p)), Q(0))


def numeric_union(boxes, theta):
    total = Q(0)
    for size in range(1, len(boxes)+1):
        for selection in combinations(boxes, size):
            volume = Q(1)
            for axis in range(2):
                lower = max(evaluate(box[axis][0], theta) for box in selection)
                upper = min(evaluate(box[axis][1], theta) for box in selection)
                volume *= max(Q(0), upper-lower)
            total += volume if size % 2 else -volume
    return total


def toy(spec):
    boxes = [tuple((trim(a), trim(b)) for a, b in box) for box in spec["boxes"]]
    full = union_polynomial(boxes)
    linearized = union_polynomial([tuple((a[:2], b[:2]) for a, b in box) for box in boxes])
    dropped = union_polynomial([boxes[i] for i in (0, 2, 3)])
    assert full == tuple(map(Q, (2, -1, 3, -1)))
    assert linearized == tuple(map(Q, (2, -1, 2)))
    assert full[:2] == linearized[:2]
    assert dropped[1] == -2
    theta = Q(spec["theta"])
    area = numeric_union(boxes, theta)
    assert area == evaluate(full, theta)
    return {"family": "one four-rectangle degeneracy family",
            "exact_area_polynomial": list(map(str, full)),
            "linearized_endpoint_union_polynomial": list(map(str, linearized)),
            "right_derivative": str(full[1]),
            "incorrect_derivative_if_first_order_zero_width_box_dropped": str(dropped[1]),
            "theta": str(theta), "exact_area_at_theta": str(area),
            "higher_order_discard_changes_no_first_order_term": True}


def main():
    started = time.monotonic()
    regpath = HERE/"registration.json"
    reg = json.loads(regpath.read_text())
    assert sha(__file__) == reg["script_sha256"]
    for source in reg["locked_sources"]:
        assert sha(source["path"]) == source["sha256"], source["path"]
    author = Path(reg["author_directory"])
    source_hash = sha(author/"first_variation.py")
    Cpath = Path(reg["standing_C_path"])
    C = json.loads(Cpath.read_text())
    denominator = int(C["root_isolation_denominator"])
    ns = [integer_fifth_root(k*denominator**5) for k in (2, 3)]
    assert ns == C["root_lower_numerators"]
    assert all(n**5 < k*denominator**5 < (n+1)**5 for n, k in zip(ns, (2, 3)))
    roots = [(Q(n, denominator), Q(n+1, denominator)) for n in ns]
    cinterval = enclosure(coefficients(C["coefficients"]), roots)
    assert cinterval == tuple(map(Q, C["ratio_interval"]))
    receipt = json.loads((author/"completion_receipt.json").read_text())
    assert all(sha(author/name) == expected for name, expected in receipt["files"].items())
    rows = []
    for name in reg["modes"]:
        d = json.loads((author/(name+"_derivative.json")).read_text())
        assert d["source_sha256"] == source_hash
        assert d["C_source_sha256"] == sha(Cpath)
        assert d["period_mass"] == "4" and d["period_mass_derivative"] == "0"
        cs = coefficients(d["right_derivative_coefficients"])
        interval = enclosure(cs, roots)
        assert interval == tuple(map(Q, d["derivative_interval"]))
        assert interval[1] < 0
        faces = [coefficients(f["coefficients"]) for f in d["coordinate_face_contributions"]]
        assert tuple(sum(v) for v in zip(*faces)) == cs
        assert all(enclosure(cs2, roots) == tuple(map(Q, face["interval"]))
                   for cs2, face in zip(faces, d["coordinate_face_contributions"]))
        if name == "first_tail":
            assert faces[2] == faces[3] == faces[4]
            assert all(enclosure(cs2, roots)[1] < 0 for cs2 in faces)
        else:
            assert not any(faces[0]) and faces[1] == faces[2] and faces[3] == faces[4]
            assert all(enclosure(cs2, roots)[1] < 0 for cs2 in faces[1:])
        inn = json.loads((author/(name+"_inner_theta_results.json")).read_text())
        out = json.loads((author/(name+"_outer_theta_results.json")).read_text())
        inner, outer = Q(inn["ratio"]), Q(out["ratio"])
        assert inner <= outer < cinterval[0]
        rows.append({"name": name, "saved_coefficients_independently_negative": True,
                     "right_derivative_interval": list(map(str, interval)),
                     "display_only": list(map(float, interval)),
                     "exact_coordinate_coefficient_sum_passes": True,
                     "saved_finite_theta_gap_to_C_lower": str(cinterval[0]-outer),
                     "saved_finite_theta_enclosure_width": str(outer-inner),
                     "geometry_to_derivative_reconstructed": False})
    out = {"status": "PASS_METHOD_AND_LOCKED_COEFFICIENT_ARITHMETIC_WITH_EXPLICIT_SCOPE_LIMIT",
           "python": sys.executable, "python_version": platform.python_version(),
           "arithmetic": "standard-library Fraction and integer comparisons only",
           "registration_sha256": sha(regpath), "script_sha256": sha(__file__),
           "locked_sources": reg["locked_sources"], "rows": rows,
           "root_lower_numerators": ns, "root_denominator": str(denominator),
           "C_coefficient_interval_exactly_reproduced": True,
           "one_toy_family": toy(reg["toy"]),
           "full_author_union_replayed": False, "old_C_union_replayed": False,
           "elapsed_seconds": time.monotonic()-started}
    path = HERE/"results.json"
    path.write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps({"status": out["status"], "toy_right_derivative": out["one_toy_family"]["right_derivative"],
                      "elapsed_seconds": out["elapsed_seconds"], "results_sha256": sha(path)}))


if __name__ == "__main__":
    main()
