#!/usr/bin/env python3
"""Compare Python r724 outputs with Wolfram exports at matching precision."""
import json
import re
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
mp.mp.dps = 110


def wolfram_mpf(value):
    """Parse Mathematica `precision` suffixes and *^ exponents at high dps."""
    s = str(value).strip().replace("−", "-")
    # In Mathematica the exponent may follow the precision annotation, e.g.
    # 2.75`75.5*^-7; retain it before stripping the annotation.
    exponent = 0
    match = re.search(r"\*\^\s*([+-]?\d+)", s)
    if match:
        exponent = int(match.group(1))
        s = s[:match.start()]
    s = s.split("`", 1)[0].strip()
    if not s or s in ("+", "-"):
        raise ValueError(f"invalid Wolfram numeric token: {value!r}")
    return mp.mpf(s) * mp.power(10, exponent)


def nstr(x):
    return mp.nstr(x, 70)


def main():
    py = json.loads((HERE / "micro.json").read_text())
    wl = json.loads((HERE / "wolfram.json").read_text())
    py_by_key = {(x["n"], x["precision_dps"]): x for x in py["results"]}
    wl_by_key = {(x["n"], x["precision"]): x for x in wl["results"]}
    fields = {
        "c": ("c", lambda x: x["baseline"]["c"]),
        "r": ("r", lambda x: x["baseline"]["r"]),
        "a_n": ("face_energy_coefficient", lambda x: x["quantities"]["a_n"]),
        "kappa_n": ("Schur_energy_coefficient", lambda x: x["quantities"]["kappa_n"]),
        "micro_energy_slope": ("curvature_slope", lambda x: x["quantities"]["micro_energy_slope"]),
    }
    rows = []
    missing = []
    for key in sorted(wl_by_key):
        n, precision = key
        w = wl_by_key[key]
        p = py_by_key.get(key)
        if p is None:
            missing.append([n, precision])
            continue
        entry = {"n": n, "precision_dps": precision, "status_python": p["status"],
                 "status_wolfram": w.get("status"), "fields": {}}
        for label, (wfield, pgetter) in fields.items():
            wv = wolfram_mpf(w[wfield])
            pv = mp.mpf(pgetter(p))
            entry["fields"][label] = {
                "wolfram_value": nstr(wv), "python_value": nstr(pv),
                "wolfram_minus_python": nstr(wv - pv),
                "abs_difference": nstr(abs(wv - pv)),
            }
        rows.append(entry)

    within_precision = []
    for n in py["dimensions"]:
        p50, p80 = py_by_key.get((n, 50)), py_by_key.get((n, 80))
        w50, w80 = wl_by_key.get((n, 50)), wl_by_key.get((n, 80))
        if not all((p50, p80, w50, w80)):
            missing.append([n, "within_software_precision_pair"])
            continue
        record = {"n": n, "python_80_minus_50": {}, "wolfram_80_minus_50": {}}
        for label, (wfield, pgetter) in fields.items():
            pdelta = mp.mpf(pgetter(p80)) - mp.mpf(pgetter(p50))
            wdelta = wolfram_mpf(w80[wfield]) - wolfram_mpf(w50[wfield])
            record["python_80_minus_50"][label] = nstr(pdelta)
            record["wolfram_80_minus_50"][label] = nstr(wdelta)
        within_precision.append(record)

    maxima = {}
    for field in fields:
        candidates = [(abs(mp.mpf(row["fields"][field]["wolfram_minus_python"])),
                       row["n"], row["precision_dps"])
                      for row in rows]
        worst = max(candidates, default=(mp.mpf(0), None, None))
        maxima[field] = {"max_abs_wolfram_minus_python": nstr(worst[0]),
                         "n": worst[1], "precision_dps": worst[2]}

    result = {
        "status": "complete" if len(rows) == len(wl["results"]) and not missing else "incomplete",
        "purpose": "High-precision numerical comparison of the independently computed micro-slope data; not a rigorous interval certificate.",
        "matched_rows": len(rows), "missing": missing,
        "parser": "Remove Mathematica precision annotation after backtick; convert *^ scientific exponent to e; parse with mpmath at 110 dps.",
        "fields": list(fields), "per_case_comparisons": rows,
        "within_software_precision_differences": within_precision,
        "maximum_same_precision_wolfram_python_differences": maxima,
    }
    (HERE / "crosscheck.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "matched_rows": len(rows),
                      "missing": missing, "maxima": maxima}, indent=2))


if __name__ == "__main__":
    main()
