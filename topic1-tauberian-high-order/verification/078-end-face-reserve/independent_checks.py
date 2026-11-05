"""Portable exact scalar checks; no continuous geometry is certified.

Adapted from reviewer script SHA256 bcf1ee759357155acc7afedd2cd80604b125a19153eb063b0143558deadad9a4.
The source/diff assertions are supplied by this package's read-only verifier.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import runpy
from contextlib import redirect_stdout
from io import StringIO

here=Path(__file__).resolve().parent
# Execute the package-local verifier without its stdout transcript; it reverses
# the staged unified diff in memory and returns old/final bytes in its namespace.
with redirect_stdout(StringIO()):
    inverse=runpy.run_path(str(here/"verify_inverse.py"))
old,new=inverse["old"],inverse["proof"]
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(old)=="e5ee441e2506b06d198b0f43572cbc8a5e51bb5a2b62aa454042c77cf1ba4b7a"
assert sha(new)=="792c6a7edd5ff3f4dc49db99bcdae91e4a500895194eec95936641d3637a0b1f"
out={"sources":{
    "v1":{"sha256":sha(old),"bytes":len(old),"LF":old.count(b"\n")},
    "v2":{"sha256":sha(new),"bytes":len(new),"LF":new.count(b"\n")}},
    "three_literal_patches_inverse_byte_identical":True,
    "face_separation_sufficient_constants":[],
    "portable_rewrite":True,
    "original_script_sha256":"bcf1ee759357155acc7afedd2cd80604b125a19153eb063b0143558deadad9a4",
    "inverse_checker_sha256":inverse["sha"](Path(inverse["__file__"]).read_bytes()),
    "geometry_certified":False}
for tau in (F(1,128),F(1,200),F(1,1000)):
    eps=tau/100
    rlo=(F(1,2)-eps)/(F(1,2)+eps)
    rhi=(F(1,2)+eps)/(F(1,2)-eps)
    # Center separation against i and membership in X_k force rlo<r_k<rhi.
    # Tangential centers cannot separate if eps<tau/4; original B_k is within
    # eps*r_k of z, still inside p's tangential margin.
    assert eps<tau/4 and eps*rhi<tau/4
    # If B_k is left of B_p, its near-unit length cannot fit into the gap.
    assert rlo>eps
    # If B_k is right of B_p, its X_k extension cannot reach the z strip.
    right_separation=(F(3,2)-eps)-eps*rhi-(F(1,2)+eps)
    assert right_separation>0
    out["face_separation_sufficient_constants"].append({
        "tau":str(tau),"epsilon":str(eps),
        "possible_root_lower":str(rlo),"possible_root_upper":str(rhi),
        "right_case_certified_positive_gap":str(right_separation)})

# Exact scalar telescope and reverse-ratio checks, independent of author WL.
b,c,v=F(3,7),F(1,4),F(1,256)
theta=b*v/c
w={h:theta**(9-h) for h in range(1,10)}
omega=sum(w.values())
f={h:F(h*h+3,h+5) for h in range(10)}
assert sum(w[h]*(c*f[h]-b*v*f[h-1]) for h in w)==c*f[9]-b*v*theta**8*f[0]
assert omega<2 and 1/omega>F(1,2)
Av,a,Q,V,M=F(3,2),F(1,64),F(1,512),F(7,3),F(2,7)
D=Av-a*(1+Q)
assert (1+Q)*(V+M)-Av*M/a == (D/a)*V*(Av/D-(1+M/V))
out["exact_telescope_and_same_M_reverse"] = True
out["terminal_reserve_at_least_one_half"] = True
out["initial_power"] = 10
out["endpoint_normalization"] = []
for n in range(2,9):
    radius=F(5,7)
    assert (v/v**n)*(v*radius)**(n-1)==radius**(n-1)
    out["endpoint_normalization"].append(n)
out["bad_owner_orders"]=[{
    "q":q,"long_young_with_bad_weight_v3":str(F(q*(2*q-1),q-1)+3),
    "good_intrinsic":3*q,"new_gap_if_E9":q-1+9}
    for q in range(3,8)]
assert out["bad_owner_orders"][0]["long_young_with_bad_weight_v3"]=="21/2"
assert out["bad_owner_orders"][0]["good_intrinsic"]==9
out["bad_q3_short_order"] = 10
out["all_assertions_passed"] = True
out["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
Path(__file__).with_suffix(".json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
