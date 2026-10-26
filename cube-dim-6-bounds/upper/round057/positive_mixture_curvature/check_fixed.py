#!/usr/bin/env python3
"""Run the seven registered algebraic configurations, without any kernel calls."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
import copy
import hashlib
import json
import sys
import time
import sympy as s
import flint
from flint import arb, ctx
from box_certificate import MODEL, MODEL_DEFINITION, certify, quadratic

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    started = time.time()
    inputs = json.loads((ROOT / "input_manifest.json").read_text())
    for path, digest in inputs["source_sha256"].items():
        assert sha(path) == digest, path
    assert not (ROOT / "results.json").exists(), "Do not overwrite a completed batch"
    checks = []

    # 1: one symbolic scalar configuration; the gaps prove global inequalities.
    z, t = s.symbols("z t", real=True, positive=True)
    score = -2*z/(t*t+z*z)
    tau = s.diff(score, z)
    gaps = [1/t**2-score**2, tau+2/t**2, 1/(4*t**2)-tau]
    expected = [(z*z-t*t)**2/(t*t*(t*t+z*z)**2),
                2*z*z*(z*z+3*t*t)/(t*t*(t*t+z*z)**2),
                (z*z-3*t*t)**2/(4*t*t*(t*t+z*z)**2)]
    assert all(s.factor(a-b) == 0 for a,b in zip(gaps, expected))
    assert tau.subs(z, 0) == -2/t**2
    assert s.simplify(tau.subs(z, s.sqrt(3)*t)) == 1/(4*t**2)
    checks.append({"configuration": 1, "gap_identities": list(map(str, expected)), "pass": True})

    # 2: exactly one six-body basis.
    N = (s.ones(6)-2*s.eye(6))*s.diag(-1,1,1,1,1,1)
    bodies = [(1,)+v for v in product((1,-1), repeat=5)]
    assert all(tuple(N.col(i)) in bodies for i in range(6))
    assert len(set(tuple(N.col(i)) for i in range(6))) == 6
    assert N.det() == 128 and N*N.T == 4*s.eye(6)+2*s.ones(6)
    Q = (N*N.T).inv()[1:6,1:6]
    assert Q == s.eye(5)/4-s.ones(5)/32
    assert Q.eigenvals() == {s.Rational(3,32):1, s.Rational(1,4):4}
    checks.append({"configuration": 2, "N_rows": [list(map(str,N.row(i))) for i in range(6)],
                   "determinant": "128", "Q_rows": [list(map(str,Q.row(i))) for i in range(5)],
                   "eigenvalues": {str(k):v for k,v in Q.eigenvals().items()}, "pass":True})

    # 3: one general symbolic face direction plus all sign corners of the same box.
    d = s.Matrix(s.symbols("d1:6", real=True))
    expr = (d.T*Q*d)[0]
    assert s.expand(expr-(sum(di**2 for di in d)/4-sum(d)**2/32)) == 0
    corner_q = [quadratic(list(map(F,bits)), "basis6") for bits in product((-1,1),repeat=5)]
    assert max(corner_q) == F(39,32) and min(corner_q) == F(15,32)
    checks.append({"configuration":3, "symbolic_quadratic":str(s.factor(expr)),
                   "unit_box_max_q6":"39/32", "unit_box_min_q6":"15/32", "pass":True})

    # 4: a fixed anisotropic box and rational log-gradient, exact corner objective.
    halfwidth = list(map(F, ["1/8","1/10","1/12","1/14","1/16"]))
    v = list(map(F, ["1/3","-1/5","2/7","-3/11","1/13"]))
    a = F(1,2)
    rows = []
    for signs in product((-1,1), repeat=5):
        delta = [sign*h for sign,h in zip(signs,halfwidth)]
        linear = sum(vi*di for vi,di in zip(v,delta))
        q = quadratic(delta,"basis6")
        rows.append((linear-q/a**2, signs, q))
    joint, winner, _ = min(rows)
    separated = -sum(abs(vi)*h for vi,h in zip(v,halfwidth))-max(row[2] for row in rows)/a**2
    assert joint >= separated
    checks.append({"configuration":4,"halfwidth":list(map(str,halfwidth)),"log_gradient":list(map(str,v)),
                   "a":str(a),"joint_corner_min":str(joint),"winning_signs":winner,
                   "separate_payment":str(separated),"improvement":str(joint-separated),"pass":True})

    # 5: algebraic interval minimization and synthetic passing interface only.
    L,U = F(2),F(5)
    branches = [(F(-1),L),(F(3),F(3)),(F(7),U)]
    assert all(min(U,max(L,b)) == expected_f for b,expected_f in branches)
    x,b = s.symbols("x b", positive=True)
    assert s.factor(s.diff(s.log(x)+b/x,x)-(x-b)/x**2) == 0
    center = ["1","1/2","1/2","1/2","1/2","1/2"]
    window = {"kind":"H","a":"1"}
    synthetic = {"model":MODEL,"model_definition":MODEL_DEFINITION,"window":window,
                 "center":center,"face_box":[["2/5","3/5"] for _ in range(5)],
                 "target":"1","quadratic":"frame32",
                 "jet":{"model":MODEL,"model_definition":MODEL_DEFINITION,
                        "window":window,"center":center,"F_interval":["2","3"],
                        "face_gradient_intervals":[["-1/10","1/10"] for _ in range(5)],
                        "source_proof":"SYNTHETIC_INTERFACE_TEST_NOT_A_BODY32_JET"}}
    passing = certify(synthetic)
    assert passing["status"] == "PASS_CONDITIONAL_ON_MATCHING_CENTER_JET"
    assert len(passing["corners"]) == 32
    assert F(passing["tangent_bound_lower"]) <= F(passing["tangent_bound_upper"])
    checks.append({"configuration":5,"shared_denominator_cases":[{"dot":str(b),"minimizer":str(f)} for b,f in branches],
                   "synthetic_only":True,"interface_status":passing["status"],
                   "bound_lower":passing["tangent_bound_lower"],"bound_upper":passing["tangent_bound_upper"],"pass":True})

    # 6: same synthetic configuration with a failing target plus metadata guards.
    failing = copy.deepcopy(synthetic)
    failing["target"] = "4"
    assert certify(failing)["status"] == "UNRESOLVED"
    rejected = []
    for label, alter in [("zero_cutoff",lambda data:data["window"].update(a="0")),
                         ("model_mismatch",lambda data:data.update(model="augmented_other_model"))]:
        invalid = copy.deepcopy(synthetic)
        alter(invalid)
        try:
            certify(invalid)
        except ValueError:
            rejected.append(label)
        else:
            raise AssertionError(label)
    checks.append({"configuration":6,"synthetic_only":True,"status":"UNRESOLVED", "rejected":rejected,"pass":True})

    # 7: fixed 32-body frame and its one orthogonal right inverse; no basis search.
    A = s.Matrix.hstack(*map(s.Matrix,bodies))
    R = A.T/32
    assert A*A.T == 32*s.eye(6)
    assert A*R == s.eye(6) and R.T*R == s.eye(6)/32
    difference = Q-s.eye(5)/32
    assert difference.eigenvals() == {s.Rational(1,16):1,s.Rational(7,32):4}
    assert F(39,32)/F(5,32) == F(39,5)
    checks.append({"configuration":7,"frame_gram":"32 I6","R":"A.T/32","face_Q":"I5/32",
                   "Q6_minus_Q32_eigenvalues":{str(k):v for k,v in difference.eigenvals().items()},
                   "isotropic_curvature_cost_ratio_Q6_over_Q32":"39/5","pass":True})

    assert len(checks) == 7
    out = {"status":"PASS_SEVEN_FIXED_CONFIGURATIONS_NO_GLOBAL_CERTIFICATE", "fixed_configurations":7,
           "checks":checks,"kernel_G_H_p_gradient_Hessian_calls":0,"fullfan_calls":0,"matlab_calls":0,
           "global_g_certificate":False,"new_C6_bound":False,"synthetic_jets_are_kernel_evidence":False,
           "python":sys.version,"sympy":s.__version__,"python_flint":flint.__version__,
           "precision_bits":ctx.prec,"elapsed_seconds":time.time()-started,
           "source_sha256":{name:sha(ROOT/name) for name in ["registration.md","theorem.md","box_certificate.py","check_fixed.py","input_manifest.json"]}}
    (ROOT/"results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["status","fixed_configurations","kernel_G_H_p_gradient_Hessian_calls","elapsed_seconds"]},indent=2))


if __name__ == "__main__":
    main()
