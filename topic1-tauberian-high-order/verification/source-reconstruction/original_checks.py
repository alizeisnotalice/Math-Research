"""Small exact checks for the original-PDF comparison, not geometric proofs."""
from fractions import Fraction as F
from pathlib import Path
import json

checks = {}
# A counterexample to the *invalid inference*, not a step used in the PDF.
centres = [F(0), F(3, 2)]
radii = [F(2), F(1)]
points = [F(0), F(3, 2), F(12, 5)]
checks["finite_centres_covered_by_first"] = all(abs(x) <= radii[0] for x in centres)
checks["original_two_boxes_cover_K"] = all(any(abs(x-c) <= r for c, r in zip(centres, radii)) for x in points)
checks["first_box_does_not_cover_K"] = any(abs(x) > radii[0] for x in points)
checks["descending_radii"] = radii[0] >= radii[1]

# Strict finite-partition budget, a rational sample of the displayed identity.
m, s, d = 7, F(2, 5), F(1, 5)
bound = d/m + (s-d)/(2*m)
checks["partition_margin_identity"] = s/m-bound == (s-d)/(2*m)
checks["partition_margin_strict"] = bound < s/m

# The periodic finite geometry, exact for these parameters.
N, rho = 13, F(2, 7)
delta, r = rho/(N*(1+2*rho)), 1/(N*(1+2*rho))
checks["own_hole_contained"] = 2*delta < r
checks["next_hole_only_touch"] = r + 2*delta == F(1, N)
checks["boundary_guard"] = r + delta < F(1, N)
checks["density_ratio"] = delta/r == rho

out = {"tool": "Python 3 fractions.Fraction", "scope": "finite rational identities only; general-parameter arguments are in original_comparison.md", "checks": checks, "all_true": all(checks.values())}
path = Path(__file__).with_suffix(".json")
path.write_text(json.dumps(out, ensure_ascii=False, indent=2)+"\n")
print(json.dumps(out, ensure_ascii=False, indent=2))
