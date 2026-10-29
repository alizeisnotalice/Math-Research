"""No new queries: rational layered variance caps from saved kernel moments."""
from pathlib import Path
from fractions import Fraction as F
import math
import json
import hashlib
import numpy as np

HERE = Path(__file__).resolve().parent


def sqrt_interval(value, scale=10**40):
    assert value >= 0
    t = math.isqrt(value.numerator*scale**2//value.denominator)
    lower = F(t, scale)
    upper = lower if lower**2 == value else F(t+1, scale)
    assert lower**2 <= value <= upper**2
    return lower, upper


def main():
    result_path = HERE/"LAYERED_RESULTS.json"
    if result_path.exists():
        raise RuntimeError("Refusing replacement of deterministic postprocessing.")
    source = HERE/"RESULTS.json"
    original = json.loads(source.read_text())
    output = []
    for case in original["cases"]:
        n, r = case["n"], case["R"]
        moments = [{(x["b"], x["class"]): x for x in rows}
                   for rows in case["kernel_moments"]]
        table_means, table_caps = [], []
        for b, segment in enumerate(case["segments"]):
            coeff = [F(1)]
            factors = []
            for i in range(n):
                minus = moments[i][b,"minus"]
                plus = moments[i][b,"plus"]
                mu_minus, mu_plus = F(minus["mean"]), F(plus["mean"])
                factor = max(F(minus["relative_variance"]), F(plus["relative_variance"]))
                factors.append(factor)
                nxt = [F(0)]*(len(coeff)+1)
                for s, c in enumerate(coeff):
                    nxt[s] += c*mu_minus
                    nxt[s+1] += c*mu_plus
                coeff = nxt
            lo, hi = segment["lower"], segment["upper"]
            mean = sum(coeff[lo:hi+1], F(0)) if lo<=hi else F(0)
            table_means.append(mean)
            table_caps.append(math.prod(1+x for x in factors)-1)
        total = sum(table_means)
        assert total == F(case["conditional_mean"])
        intervals = [sqrt_interval(v) for v in table_caps]
        layer_lower = (sum((w*a for w,(a,z) in zip(table_means,intervals)),F(0))/total)**2
        layer_upper = (sum((w*z for w,(a,z) in zip(table_means,intervals)),F(0))/total)**2
        actual = F(case["relative_variance"])
        assert actual <= layer_upper <= F(case["relative_variance_cap"])
        output.append({"n": n, "R": r, "L": case["L"],
            "table_gated_conditional_means": list(map(str,table_means)),
            "table_relative_variance_caps": list(map(str,table_caps)),
            "layered_cap_lower_rational": str(layer_lower),
            "layered_cap_upper_rational": str(layer_upper),
            "actual_le_certified_layered_upper": True,
            "layered_upper_le_global_cap": True,
            "layered_cap_float": float(layer_upper),
            "global_cap_float": case["cap_float"],
            "actual_relative_variance_float": float(actual)})
    out = {"scope": "saved-kernel postprocessing; zero new queries or samples",
           "input_SHA256": hashlib.sha256(source.read_bytes()).hexdigest(),
           "script_SHA256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "sqrt_rational_scale": "10^40", "numpy": np.__version__, "cases": output}
    result_path.write_text(json.dumps(out,separators=(",",":")))
    print(json.dumps([{k:c[k] for k in ["n","R","L","layered_cap_float",
        "global_cap_float","actual_relative_variance_float"]} for c in output],indent=2))


if __name__ == "__main__":
    main()
