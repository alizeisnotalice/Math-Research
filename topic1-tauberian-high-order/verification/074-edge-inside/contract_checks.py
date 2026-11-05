"""Portable R074 byte/finite arithmetic verification.

Only file discovery and in-memory inverse reconstruction differ from the
reviewer's original local script; monomial/derivative assertions are unchanged.
No continuous geometry is encoded.
"""
from pathlib import Path
import hashlib,json
from fractions import Fraction as F
from math import comb
base=Path(__file__).resolve().parent
sha=lambda data:hashlib.sha256(data).hexdigest()
v3=(base.parents[1]/"proofs/074-edge-inside.md").read_bytes()
assert sha(v3)=="57514c9109718eb904bb47088fb995b8733ebf302200f269a4287c90b42c6424"
states={"ninth_mechanisms_and_inside.v3.md":v3}
out={"reverse_patches":[]}
for patch_name in ["v3_source_contract_patch.json","v2_boundary_patch.json"]:
    patch=json.loads((base/patch_name).read_text())
    new=states[Path(patch["new_path"]).name]
    assert sha(new)==patch["new_sha256"]
    recovered=new.decode()
    for block in reversed(patch["blocks"]):
        assert recovered.count(block["new"])==1
        recovered=recovered.replace(block["new"],block["old"],1)
    old=recovered.encode();assert sha(old)==patch["old_sha256"]
    states[Path(patch["old_path"]).name]=old
    out["reverse_patches"].append({"patch":patch_name,"inverse_byte_identical":True})
out["sources"]={name:{"sha256":sha(data),"bytes":len(data),"LF":data.count(b"\n")} for name,data in states.items()}
# Exact monomial integration in the two shallow-axis orderings:
# integral_0^Q P(Q-y)dy = P Q^2(1-1/2), and the other ordering swaps P,Q.
# The dictionary uses powers (P,Q); f is left as a common denominator.
one_endpoint_numerator = {(1, 2): F(1)-F(1, 2),
                          (2, 1): F(1)-F(1, 2)}
assert one_endpoint_numerator == {(1, 2): F(1, 2), (2, 1): F(1, 2)}
# Two endpoints, free source width vR and f=3v^2R: after division by v^8R^3.
per_edge = {powers: F(2,3)*coefficient
            for powers, coefficient in one_endpoint_numerator.items()}
assert per_edge == {(1, 2): F(1,3), (2, 1): F(1,3)}
assert sum(per_edge.values()) == F(2, 3)
assert 3*sum(per_edge.values()) == 2
out["nominal_n3_whole_pair_physical_coefficient"] = {
    f"p^{i}q^{j}": str(c) for (i,j),c in per_edge.items()}
out["canonical_coefficient_times_three_slots"] = "2"

# n=3 exact deletion derivative: write m=R^3u^3 and expand f=R^3(1+u)^3.
# Divide df/du by dm/du=3R^3u^2, then multiply by v^2. With w=v/u,
# the coefficients of v^k w^(2-k) must be those of (v+w)^2.
relative_coefficients = {}
for power in range(1,4):
    coefficient = F(power*comb(3,power),3)
    relative_coefficients[(power-1,3-power)] = coefficient
assert relative_coefficients == {(0,2): F(1), (1,1): F(2), (2,0): F(1)}
out["n3_derivative_relative_to_v_over_a"] = "(v+z^(-1/3))^2"
out["derivative_exact_polynomial_coefficients_v_w"] = {
    f"v^{i}w^{j}": str(c) for (i,j),c in relative_coefficients.items()}
out["uniform_interval"] = "z in [1-A*v,1] gives relative derivative 1+O_A(v)"

out["script_sha256"] = sha(Path(__file__).read_bytes())
out["all_assertions_passed"] = True
target = Path(__file__).with_suffix(".json")
target.write_text(json.dumps(out, indent=2)+"\n")
print(json.dumps(out, indent=2))
