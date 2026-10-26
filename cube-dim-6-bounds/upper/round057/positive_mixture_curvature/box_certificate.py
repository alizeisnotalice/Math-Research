#!/usr/bin/env python3
"""A057 conditional whole-face-box certificate; no kernel evaluator is present."""
from fractions import Fraction as Rat
from itertools import product
import json
import sys
from flint import arb, ctx

MODEL = "body32_unit_cauchy_v1"
MODEL_DEFINITION = "V={(1,s2,...,s6):si=+-1}; f_t(z)=t/(pi*(t^2+z^2)); 32 independent unit-weight factors"


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise ValueError("Exact rational strings or integers are required")
    return Rat(value)


def interval(raw):
    if not isinstance(raw, list) or len(raw) != 2:
        raise ValueError("Each interval needs exactly two rational endpoints")
    lo, hi = map(rational, raw)
    if lo > hi:
        raise ValueError("Reversed interval")
    return lo, hi


def ball(q):
    return arb(q.numerator) / arb(q.denominator)


def quadratic(d, mode):
    ss = sum(x*x for x in d)
    if mode == "frame32":
        return ss / 32
    if mode == "basis6":
        return ss / 4 - sum(d)**2 / 32
    raise ValueError("Use one fixed proved quadratic: frame32 or basis6")


def certify(data, precision_bits=384):
    """Consumes a jet interval box; PASS remains conditional on the jet proof."""
    if data.get("model") != MODEL:
        raise ValueError("Mismatched model")
    if data.get("model_definition") != MODEL_DEFINITION:
        raise ValueError("Mismatched model definition/normalization")
    window = data["window"]
    a = rational(window["a"])
    if a <= 0:
        raise ValueError("Strictly positive time lower endpoint required")
    if window["kind"] == "H":
        if "T" in window:
            raise ValueError("H window must not carry a finite T")
    elif window["kind"] == "B":
        if rational(window["T"]) <= a:
            raise ValueError("Band requires T>a")
    else:
        raise ValueError("Only H_a or B_(a,T) is supported")
    c = list(map(rational, data["center"]))
    if len(c) != 6 or c[0] != 1:
        raise ValueError("Center must be six-dimensional with x1=1")
    box = list(map(interval, data["face_box"]))
    if len(box) != 5:
        raise ValueError("Exactly five face intervals required")
    if any(not lo <= ci <= hi for (lo, hi), ci in zip(box, c[1:])):
        raise ValueError("Registered center must belong to the face box")
    jet = data["jet"]
    if jet["model"] != data["model"] or jet["model_definition"] != data["model_definition"]:
        raise ValueError("Jet model mismatch")
    if jet["window"] != window or list(map(rational, jet["center"])) != c:
        raise ValueError("Jet/window/center mismatch")
    if not isinstance(jet.get("source_proof"), str) or not jet["source_proof"]:
        raise ValueError("The supplied jet must identify its external proof source")
    lo_f, hi_f = interval(jet["F_interval"])
    if lo_f <= 0:
        raise ValueError("Strictly positive center-value lower endpoint required")
    grad = list(map(interval, jet["face_gradient_intervals"]))
    if len(grad) != 5:
        raise ValueError("Exactly five matching face derivative intervals required")
    target = rational(data["target"])
    if target <= 0:
        raise ValueError("Positive target required")
    mode = data.get("quadratic", "frame32")
    quadratic([Rat(0)]*5, mode)
    if not 64 <= precision_bits <= 4096:
        raise ValueError("Precision must be between 64 and 4096 bits")
    ctx.prec = precision_bits
    rows = []
    lower_endpoints = []
    upper_endpoints = []
    for bits in product((0, 1), repeat=5):
        d = [box[i][bit]-c[i+1] for i, bit in enumerate(bits)]
        b = sum(min(gl*di, gu*di) for (gl, gu), di in zip(grad, d))
        f = min(hi_f, max(lo_f, b))
        q = quadratic(d, mode)
        exponent = b/f-q/a**2
        enclosure = ball(f)*ball(exponent).exp()
        if not enclosure.is_finite():
            raise ValueError("Nonfinite exponential enclosure")
        lower = Rat(str(enclosure.lower().fmpq()))
        upper = Rat(str(enclosure.upper().fmpq()))
        lower_endpoints.append(lower)
        upper_endpoints.append(upper)
        rows.append({"corner_bits": list(bits), "displacement": list(map(str, d)),
                     "gradient_dot_lower": str(b), "minimizing_center_value": str(f),
                     "quadratic": str(q), "exponent": str(exponent),
                     "certificate_lower": str(lower), "certificate_upper": str(upper),
                     "passes_target": lower >= target})
    lower = min(lower_endpoints)
    upper = min(upper_endpoints)
    return {"status": "PASS_CONDITIONAL_ON_MATCHING_CENTER_JET" if lower >= target else "UNRESOLVED",
            "conditional_on_jet_source": jet["source_proof"],
            "scope": "one whole box on x1=1; no global coverage assertion",
            "model": MODEL, "window": window, "quadratic": mode,
            "target": str(target), "precision_bits": precision_bits,
            "tangent_bound_lower": str(lower), "tangent_bound_upper": str(upper),
            "corners": rows, "kernel_calls": 0, "global_g_certificate": False}


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: box_certificate.py INPUT.json")
    with open(sys.argv[1]) as stream:
        data = json.load(stream)
    print(json.dumps(certify(data), indent=2))
