#!/usr/bin/env python3
"""r060: split-branch, fixed-c dual reoptimization for 33 receivers.

Floating-point diagnostic only.  The t2 integral is split at S=1, n and
Phi'(cap); the outer t1 integral is composite Gauss--Legendre.  No tensor
grid is used.  The dual is scaled by exp(tau), so dD/dw_g is
coeff_g * (q - average_g).
"""
from __future__ import annotations

import json
import math
import sys
import time
import argparse
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammainc, gammaln


# Set from --output-dir in main.  A temporary default prevents an ordinary
# invocation from overwriting checked-in or shared-worktree artifacts.
OUT = Path(tempfile.gettempdir()) / "r060-compute"
N = 64
B = 32
CAP = 2.0
Q = 1.00001
C = 1.2360667450799148
TAU = 66.43897513452868
GAMMAS = (0.2, 0.4, 0.6, 0.8)
V = np.arange(17, dtype=float)
COEFF = np.r_[1.0, np.full(16, 2.0)]
BINOM = np.array([math.comb(B, k) for k in range(B + 1)], dtype=float)
KVAL = np.arange(B + 1, dtype=float)
POLYNORM = math.exp(-gammaln(B)) / (1.0 + math.exp(-C)) ** B
ALPHA_CAP = (1.0 + math.log(CAP)) ** (N - 2) * (N + math.log(CAP))


def jsonable(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.integer):
        return int(x)
    raise TypeError(type(x).__name__)


def gl_composite(a: float, b: float, panels: int, order: int = 16):
    xx, ww = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(a, b, panels + 1)
    mid = (edges[:-1] + edges[1:]) / 2.0
    half = (edges[1:] - edges[:-1]) / 2.0
    return (mid[:, None] + half[:, None] * xx).ravel(), (half[:, None] * ww).ravel()


@dataclass
class BranchEvaluator:
    gamma: float
    outer_nodes: int
    soft_order: int = 16

    def __post_init__(self):
        d = self.gamma * V * V / N
        self.u = np.vstack((np.array([[0.0, 0.0]]),
                            np.column_stack((d[1:] + V[1:], d[1:] - V[1:])),
                            np.column_stack((d[1:] - V[1:], d[1:] + V[1:]))))
        # This is the same finite range design used by r032, with a generous
        # 16-sigma upper edge.  The reported values are not tail-certified.
        p = 1.0 / (1.0 + math.exp(C))
        mz = 1.0 + C * p
        vz = 1.0 + C * C * p * (1.0 - p)
        self.lo = float(np.min(self.u))
        self.hi = float(np.max(self.u) + B * mz + 16.0 * math.sqrt(B * vz))
        panels = self.outer_nodes // 16
        assert panels * 16 == self.outer_nodes
        self.t1, self.t1w = gl_composite(self.lo, self.hi, panels)
        self.gx, self.gw = np.polynomial.legendre.leggauss(self.soft_order)
        self.eval_count = 0

    def poly(self, x):
        """Polynomial in e^t f(t), evaluated for scalar/vector x."""
        x = np.asarray(x, dtype=float)
        z = np.maximum(x[..., None] - C * KVAL, 0.0)
        return POLYNORM * np.sum(BINOM * z ** (B - 1), axis=-1)

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        z = np.maximum(x[..., None] - C * KVAL, 0.0)
        weights = BINOM * np.exp(-C * KVAL) / (1.0 + math.exp(-C)) ** B
        return np.sum(weights * gammainc(B, z), axis=-1)

    @staticmethod
    def profile(s):
        s = np.asarray(s, dtype=float)
        ans = np.zeros_like(s)
        middle = (s >= 1.0) & (s <= N)
        ans[middle] = 1.0
        soft = (s > N) & (s < ALPHA_CAP)
        if np.any(soft):
            ls = np.log(s[soft])
            y = np.maximum(np.exp((ls - math.log(N)) / (N - 1)) - 1.0, 0.0)
            for _ in range(14):
                one = 1.0 + y
                f = (N - 2) * np.log(one) + np.log(N + y) - ls
                y = np.maximum(y - f / ((N - 2) / one + 1.0 / (N + y)), 0.0)
            ans[soft] = np.exp(y)
        ans[s >= ALPHA_CAP] = CAP
        return ans

    def like(self, t2, amp):
        return self.poly(np.asarray(t2)[..., None] - self.u[:, 1]) @ amp

    def boundary(self, level, amp):
        left, right = self.lo, self.hi
        if self.like(left, amp) >= level:
            return left
        if self.like(right, amp) <= level:
            return right
        for _ in range(52):
            mid = (left + right) / 2.0
            if self.like(mid, amp) < level:
                left = mid
            else:
                right = mid
        return (left + right) / 2.0

    def constant_piece(self, t1, left, right, h, f1):
        if right <= left:
            return 0.0, np.zeros(33)
        cost = math.exp(-t1) * (math.exp(-left) - math.exp(-right)) * float(self.phi(h))
        av = h * f1 * (self.cdf(right - self.u[:, 1]) - self.cdf(left - self.u[:, 1]))
        return cost, av

    @staticmethod
    def phi(h):
        h = np.asarray(h, dtype=float)
        return h * (1.0 + np.maximum(np.log(np.maximum(h, 1.0)), 0.0)) ** (N - 1)

    def evaluate(self, w, details=False):
        w = np.maximum(np.asarray(w, dtype=float), 0.0)
        fullw = np.r_[w[0], w[1:], w[1:]]
        total_cost = 0.0
        full_av = np.zeros(33)
        for t1, wt1 in zip(self.t1, self.t1w):
            p1 = self.poly(t1 - self.u[:, 0])
            f1 = np.exp(self.u[:, 0] - t1) * p1
            amp = fullw * p1 * np.exp(np.sum(self.u, axis=1) - TAU)
            roots = [self.boundary(level, amp) for level in (1.0, float(N), ALPHA_CAP)]
            c0, a0 = self.constant_piece(t1, roots[0], roots[1], 1.0, f1)
            c2, a2 = self.constant_piece(t1, roots[2], self.hi, CAP, f1)
            total_cost += wt1 * (c0 + c2)
            full_av += wt1 * (a0 + a2)
            if roots[2] > roots[1]:
                t2 = (roots[2] - roots[1]) * self.gx / 2.0 + (roots[2] + roots[1]) / 2.0
                wt2 = (roots[2] - roots[1]) * self.gw / 2.0
                s = self.like(t2, amp)
                h = self.profile(s)
                total_cost += wt1 * np.sum(wt2 * np.exp(-t1 - t2) * self.phi(h))
                f2 = np.exp(self.u[:, 1][None, :] - t2[:, None]) * self.poly(t2[:, None] - self.u[:, 1][None, :])
                full_av += wt1 * np.sum((wt2 * h)[:, None] * f2 * f1[None, :], axis=0)
        groups = np.r_[full_av[0], full_av[1:17] + full_av[17:]] / COEFF
        dual = Q * float(COEFF @ w) + math.exp(TAU) * total_cost - float(fullw @ full_av)
        self.eval_count += 1
        ans = {"direct_cost": float(total_cost), "logJ": float(math.log(total_cost)),
               "averages_by_group": groups, "dual_scaled": float(dual),
               "dual_gradient": COEFF * (Q - groups), "full_averages": full_av,
               # The two mirrored families are accumulated separately; they
               # are reported as a direct numerical mirror check, not copied.
               "plus_averages_direct": full_av[1:17],
               "minus_averages_direct": full_av[17:],
               "mirror_average_max_abs_difference": float(np.max(np.abs(full_av[1:17] - full_av[17:])))}
        return ans if details else (dual, COEFF * (Q - groups))

    def volume(self):
        """Exact finite union of upper orthants under exp(-t1-t2) dt."""
        pts = self.u[np.argsort(self.u[:, 0], kind="mergesort")]
        vol = 0.0
        floor = math.inf
        i = 0
        while i < len(pts):
            x = pts[i, 0]
            j = i
            newfloor = floor
            while j < len(pts) and pts[j, 0] == x:
                newfloor = min(newfloor, float(pts[j, 1]))
                j += 1
            nextx = float(pts[j, 0]) if j < len(pts) else math.inf
            floor = newfloor
            vol += (math.exp(-x) - (math.exp(-nextx) if math.isfinite(nextx) else 0.0)) * math.exp(-floor)
            i = j
        return vol


def run_gamma(gamma, start, outer=128):
    ev = BranchEvaluator(gamma, outer)
    calls = 0
    def fun(w):
        nonlocal calls
        d, grad = ev.evaluate(w)
        calls += 1
        return -d, -grad
    res = minimize(fun, np.maximum(start, 0.0), jac=True, method="L-BFGS-B",
                   bounds=[(0.0, None)] * 17,
                   options={"maxiter": 110, "ftol": 1e-12, "gtol": 3e-7,
                            "maxls": 30, "maxfun": 360})
    det = ev.evaluate(res.x, details=True)
    vol = ev.volume()
    record = {
        "gamma": gamma, "outer_nodes_optimized": outer, "soft_order_optimized": ev.soft_order,
        "finite_window": [ev.lo, ev.hi], "weights": res.x,
        "directcost": det["direct_cost"], "logJ": det["logJ"],
        "allaverages": det["averages_by_group"], "min_average": float(np.min(det["averages_by_group"])),
        "plus_averages_direct": det["plus_averages_direct"], "minus_averages_direct": det["minus_averages_direct"],
        "mirror_average_max_abs_difference": det["mirror_average_max_abs_difference"],
        "dual_scaled": det["dual_scaled"], "dualgradient": det["dual_gradient"],
        "gradient_inf_norm": float(np.max(np.abs(det["dual_gradient"]))),
        "optimizerstatus": {"success": bool(res.success), "status": int(res.status),
                            "message": str(res.message), "nit": int(res.nit),
                            "nfev": int(res.nfev), "objective_calls": calls},
        "coneunionvolume_exact_finitepoints": vol, "logV": math.log(vol),
        "score_logV_minus_logJ": math.log(vol) - det["logJ"],
        "gap_to_71_79": 71.79 - (math.log(vol) - det["logJ"]),
    }
    return record


def refine(gamma, w):
    """Refine at 256, common-rescale to q if necessary, then replay at 512."""
    ev256 = BranchEvaluator(gamma, 256)
    pre = ev256.evaluate(w, details=True)
    scale = 1.0
    if float(np.min(pre["averages_by_group"])) < Q:
        low, high = 1.0, 1.02
        highdet = ev256.evaluate(high * w, details=True)
        while float(np.min(highdet["averages_by_group"])) < Q and high < 8.0:
            low, high = high, 1.5 * high
            highdet = ev256.evaluate(high * w, details=True)
        if float(np.min(highdet["averages_by_group"])) < Q:
            raise RuntimeError("could not bracket q in common-weight rescale")
        for _ in range(42):
            mid = (low + high) / 2.0
            if float(np.min(ev256.evaluate(mid * w, details=True)["averages_by_group"])) < Q:
                low = mid
            else:
                high = mid
        scale = high
    corrected = scale * w
    rows = []
    for outer in (256, 512):
        ev = ev256 if outer == 256 else BranchEvaluator(gamma, outer)
        det = ev.evaluate(corrected, details=True)
        vol = ev.volume()
        rows.append({"outer_nodes": outer, "soft_order": ev.soft_order,
                     "same_sourceweights": True, "directcost": det["direct_cost"],
                     "logJ": det["logJ"], "allaverages": det["averages_by_group"],
                     "min_average": float(np.min(det["averages_by_group"])),
                     "dual_scaled": det["dual_scaled"], "dualgradient": det["dual_gradient"],
                     "gradient_inf_norm": float(np.max(np.abs(det["dual_gradient"]))),
                     "plus_averages_direct": det["plus_averages_direct"],
                     "minus_averages_direct": det["minus_averages_direct"],
                     "mirror_average_max_abs_difference": det["mirror_average_max_abs_difference"],
                     "coneunionvolume_exact_finitepoints": vol,
                     "score_logV_minus_logJ": math.log(vol) - det["logJ"],
                     "gap_to_71_79": 71.79 - (math.log(vol) - det["logJ"])})
    return {"initial_256_min_average": float(np.min(pre["averages_by_group"])),
            "common_weight_scale_at_256": scale, "corrected_weights": corrected,
            "rows": rows}


def write_checkpoint(payload):
    (OUT / "compute60.json").write_text(json.dumps(payload, indent=2, default=jsonable) + "\n")


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path,
                    help="directory for compute60.json/.md (default: a new temporary directory)")
    ap.add_argument("--replay-savedweights", type=Path,
                    help="load a prior compute60 JSON and replay corrected weights at 512 outer nodes with soft GL16/64")
    args = ap.parse_args()
    OUT = args.output_dir if args.output_dir else Path(tempfile.mkdtemp(prefix="r060-compute-"))
    OUT.mkdir(parents=True, exist_ok=True)
    if args.replay_savedweights:
        replay(args.replay_savedweights)
        return
    started = time.time()
    payload = {"status": "running floating-point split-branch fixed-c reoptimization; not a certificate",
               "parameters": {"n": N, "b": B, "cap": CAP, "q": Q, "c": C, "tau": TAU,
                              "v_nonnegative": V, "d_formula": "gamma*v^2/64",
                              "outer_rule": "composite Gauss-Legendre order 16", "soft_order": 16},
               "results": [], "started_unix": started}
    # Required deterministic initial point: a center-only multiplier.  Later
    # gamma values warm-start from the preceding solution only.
    start = np.zeros(17)
    start[0] = 1.0
    for gamma in GAMMAS:
        rec = run_gamma(gamma, start)
        payload["results"].append(rec)
        payload["elapsed_seconds"] = time.time() - started
        write_checkpoint(payload)  # durable per-gamma checkpoint
        print(json.dumps({"gamma": gamma, "score": rec["score_logV_minus_logJ"],
                          "minavg": rec["min_average"], "nit": rec["optimizerstatus"]["nit"]}), flush=True)
        start = np.asarray(rec["weights"])
    # Refine the best selected input; these runs deliberately do NOT optimize.
    best = max(payload["results"], key=lambda r: r["score_logV_minus_logJ"])
    best["fixed_weight_refinements"] = refine(best["gamma"], np.asarray(best["weights"]))
    payload["best_by_128node_score"] = best["gamma"]
    payload["status"] = "complete floating-point diagnostic; finite-window quadrature and optimizer results are not a certificate"
    payload["elapsed_seconds"] = time.time() - started
    write_checkpoint(payload)
    (OUT / "compute60.md").write_text(report(payload))
    print(json.dumps({"best_gamma": best["gamma"], "best_128_score": best["score_logV_minus_logJ"],
                      "refined": best["fixed_weight_refinements"]}), flush=True)


def report(p):
    lines = ["# r060 split-branch fixed-c reoptimization", "",
             "Floating-point evidence only: no interval enclosure, no certified tail bound, and no claim of global optimization.", "",
             f"Fixed `N=64`, `b=32`, `cap=2`, `q={Q}`, `c={C:.16g}`, `tau={TAU:.16g}`. "
             "For each gamma, 33 receivers use `v=0,...,16`, `d=gamma v^2/64`, and symmetric `(d+v,d-v)` points.", "",
             "The evaluator splits every t2 fibre at `S=1`, `S=N`, and `S=Phi'(cap)`; constant branches use Gamma--binomial CDF differences and only the soft branch uses GL16. The `256` and `512` rows below preserve, rather than reoptimize, the selected 128-node weights.", "",
             "| gamma | log J (128) | min average | score logV-logJ | gap to 71.79 | optimizer |", "|---:|---:|---:|---:|---:|---|"]
    for r in p["results"]:
        st = r["optimizerstatus"]
        lines.append(f"| {r['gamma']:.1f} | {r['logJ']:.10f} | {r['min_average']:.10f} | {r['score_logV_minus_logJ']:.10f} | {r['gap_to_71_79']:.10f} | {st['status']}: {st['message']} |")
    best = max(p["results"], key=lambda r: r["score_logV_minus_logJ"])
    lines += ["", f"Selected by the 128-node floating score: gamma `{best['gamma']}`. Its exact finite-point cone-union formula (no spatial grid) gave volume `{best['coneunionvolume_exact_finitepoints']:.16g}`.", "",
              ]
    corr = best["fixed_weight_refinements"]
    lines += ["", f"At 256 nodes the selected weights had minimum average `{corr['initial_256_min_average']:.12f}`. A common multiplier `{corr['common_weight_scale_at_256']:.15g}` was bisectioned against actual branch-split averages to reach `q`; the displayed 256 and 512 values use those corrected weights."]
    lines += ["", "| preserved weights outer nodes | log J | min average | score | gap to 71.79 |", "|---:|---:|---:|---:|---:|"]
    for r in corr["rows"]:
        lines.append(f"| {r['outer_nodes']} | {r['logJ']:.12f} | {r['min_average']:.12f} | {r['score_logV_minus_logJ']:.12f} | {r['gap_to_71_79']:.12f} |")
    lines += ["", "A positive reported gap means the floating computation is strictly below 71.79 at the displayed quadrature; it is not a numerical certificate.", ""]
    return "\n".join(lines)


def replay(saved_path: Path):
    """Independent, no-optimizer replay of saved corrected weights."""
    src = json.loads(saved_path.read_text())
    gamma = src["best_by_128node_score"]
    rec = next(x for x in src["results"] if float(x["gamma"]) == float(gamma))
    ref = rec.get("fixed_weight_refinements", {})
    w = np.asarray(ref.get("corrected_weights", rec["weights"]), dtype=float)
    rows = []
    for inner in (16, 64):
        ev = BranchEvaluator(float(gamma), 512, inner)
        det = ev.evaluate(w, details=True)
        vol = ev.volume()
        rows.append({"outer_nodes": 512, "soft_order": inner, "directcost": det["direct_cost"],
                     "logJ": det["logJ"], "allaverages": det["averages_by_group"],
                     "min_average": float(np.min(det["averages_by_group"])), "dual_scaled": det["dual_scaled"],
                     "dualgradient": det["dual_gradient"], "gradient_inf_norm": float(np.max(np.abs(det["dual_gradient"]))),
                     "plus_averages_direct": det["plus_averages_direct"], "minus_averages_direct": det["minus_averages_direct"],
                     "mirror_average_max_abs_difference": det["mirror_average_max_abs_difference"],
                     "score_logV_minus_logJ": math.log(vol)-det["logJ"], "gap_to_71_79": 71.79-(math.log(vol)-det["logJ"])})
    payload = {"status": "no-optimizer fixed-weight replay; floating diagnostic only", "source": str(saved_path),
               "gamma": gamma, "weights": w, "rows": rows}
    (OUT / "compute60_replay.json").write_text(json.dumps(payload, indent=2, default=jsonable) + "\n")
    print(json.dumps({"output_dir": str(OUT), "replay": rows}, default=jsonable), flush=True)


if __name__ == "__main__":
    main()
