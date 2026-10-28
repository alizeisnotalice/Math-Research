"""Exact multinomial bin labels + analytic conditional means.

These are not direct high-dimensional point samples. True cube-tail lower
estimates are obtained by conditional variance bounds and Cantelli.
"""
import os
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
from pathlib import Path
from fractions import Fraction
import sys
import json
import math
import time
import numpy as np

OUT = Path(__file__).resolve().parent
SEED = 20261029032
P = .5
BN = 4
LNUM = np.arange(80, 89)
LENGTHS = LNUM / 72
LOG_THRESHOLDS = np.arange(-2, 12.00001, .25)
THRESHOLDS = np.exp(LOG_THRESHOLDS)


def lse(x, axis=-1):
    m = np.max(x, axis=axis, keepdims=True)
    return np.squeeze(m + np.log(np.sum(np.exp(x - m), axis=axis,
                                       keepdims=True)), axis)


def field_labels(d, polynomial):
    q = 1 << d
    k = q - 1
    ex = []
    z = 1
    for _ in range(k):
        ex.append(z)
        z <<= 1
        if z & q:
            z ^= polynomial
    assert z == 1 and len(set(ex)) == k
    ex = np.array(ex, dtype=np.int64)
    log = np.full(q, -1, dtype=np.int64)
    log[ex] = np.arange(k)
    w = [1, 6, 7]
    products = [int(ex[(log[a] + log[b]) % k])
                for i, a in enumerate(w) for b in w[i:]]
    assert len(set(products)) == 6
    defects = {(int(log[a]) - int(log[b])) % k
               for a in w for b in w if a != b}
    assert len(defects) == 6
    assert not any(j % k in defects for j in [1, -1, 2, -2])
    expanded = defects | {0}
    assert not any((j + 1) % k in expanded for j in expanded)
    def trace(a):
        value, cur = a.copy(), a.copy()
        for _ in range(d - 1):
            cur = ex[2 * log[cur] % k]
            value ^= cur
        assert np.all((value == 0) | (value == 1))
        return value
    labels = trace(ex) + 2 * trace(ex[(log[ex] + log[6]) % k])
    histogram = np.bincount(labels, minlength=4)
    assert np.array_equal(histogram, [q // 4 - 1, q // 4, q // 4, q // 4])
    return k, labels, {"degree": d, "polynomial_hex": hex(polynomial),
                       "basis": [1, 6], "defects": sorted(defects),
                       "Sidon_verified": True, "short_lags_verified": True,
                       "nonadjacency_verified": True,
                       "q_gt_H4_129": q > 129,
                       "histogram": histogram.tolist()}


class Source:
    def __init__(self, degree, polynomial, v, rmax):
        self.k, labels, self.field = field_labels(degree, polynomial)
        self.v, self.rmax = v, rmax
        self.n = 116 * v * (self.k // 3) ** 2
        self.delta = 1 / (4 * self.k * v)
        self.eps = P / self.k
        self.m = math.ceil(math.sqrt(v))
        self.q = math.lcm(8 * self.k * v, 144, rmax)
        self.g = self.q // (4 * self.k * v)  # delta*Q, integer
        assert self.q % rmax == 0 and self.q % 144 == 0
        assert self.g >= 2
        atom_num = 8 * np.arange(self.k) + 2 * labels + 1
        self.centers2 = atom_num * (self.q // (4 * self.k))
        self.left2 = self.centers2 - self.g
        self.right2 = self.centers2 + self.g
        assert self.left2.min() > 0 and self.right2.max() < 2 * self.q
        assert np.all(self.right2[:-1] < self.left2[1:])
        self.mid2 = 2 * np.arange(self.q, dtype=np.int64) + 1
        mf = np.arange(1, self.m + 1, dtype=float)
        b = 4 * (729 / 29) * (
            (np.sin(np.pi * mf * 10 / 27) / (np.pi * mf)) ** 4
            + 3 * (np.sin(np.pi * mf / 9) / (np.pi * mf)) ** 4)
        self.c = np.sqrt(2 * b / 3)
        self.cosines = np.cos(np.pi * self.mid2[:, None] * mf / self.q)
        # pi < 22/7 and c_m < 20/(9m^2), hence exact rational upper bound.
        db_fraction = Fraction(880, 63) * sum(
            (Fraction(1, m) for m in range(1, self.m + 1)), Fraction(0))
        self.db = float(db_fraction)
        self.t0 = (self.eps * self.db / (2 * self.q)) ** 2
        self.t1 = (self.eps * (self.db + 1 / self.delta)
                   / (2 * self.q)) ** 2
        harmonic = sum(1 / m for m in range(1, self.m + 1))
        self.bg_lower = 1 - P - (20 * P / (3 * self.k)) * harmonic
        assert self.bg_lower > 0
        self.r_values = [r for r in [1, 4, 8, 16, 32] if r <= rmax]

    def primitive_numerator(self, f2):
        """P_delta=fraction numerator/(2Q), integer operations before division."""
        done = np.searchsorted(self.right2, f2, side="right")
        candidate = np.minimum(done, self.k - 1)
        partial = np.clip(f2 - self.left2[candidate], 0, 2 * self.g)
        partial[done == self.k] = 0
        scaled_count = done * (2 * self.g) + partial
        peak = ((done < self.k) & (f2 > self.left2[candidate])
                & (f2 < self.right2[candidate]))
        pnum = scaled_count * (self.q // self.g) - self.k * f2
        return pnum, peak

    def base_query(self, li):
        lq = int(LNUM[li]) * (self.q // 72)
        upper = (self.mid2 + lq) % (2 * self.q)
        lower = (self.mid2 - lq) % (2 * self.q)
        pu, iu = self.primitive_numerator(upper)
        pl, il = self.primitive_numerator(lower)
        h = (pu - pl) / (2 * self.q)
        freq = np.arange(1, self.m + 1, dtype=float)
        # Mean cosine in a bin: center cosine * sinc(pi*m/Q).
        coefficient = (self.c * np.sin(np.pi * freq * LENGTHS[li])
                       * np.sinc(freq / self.q))
        bmean = self.cosines @ coefficient
        means = 1 + self.eps / LENGTHS[li] * (h + bmean)
        assert means.min() >= 1 / LENGTHS[li] - 1e-12
        return means, iu != il

    def query(self, li):
        base, transition = self.base_query(li)
        shifts = np.arange(self.rmax) * (self.q // self.rmax)
        return (np.column_stack([np.roll(base, int(s)) for s in shifts]),
                np.column_stack([np.roll(transition, int(s)) for s in shifts]))

    def metadata(self):
        return {"K": self.k, "v": self.v, "n": self.n, "M": self.m, "B": 4,
                "p": P, "w": 4, "delta": self.delta, "R_max": self.rmax,
                "R_values": self.r_values, "bins_Q": self.q,
                "background_lower_bound": self.bg_lower, "field": self.field,
                "v_ge_16logn": self.v >= 16 * math.log(self.n),
                "B_ge_v": 4 >= self.v,
                "original_matched_lower_bound_claimed": False,
                "conditional_variance_DB_upper": self.db}


def draw_counts(src, rng, samples):
    counts = np.empty((samples, src.q), dtype=float)
    uniform = np.full(src.q, 1 / src.q)
    proposal_types = []
    probability_sum_residuals = []
    normalization_residuals = []
    for row in range(samples):
        if rng.random() < .5:
            probs = uniform
            proposal_types.append([-1, -1])
        else:
            li = int(rng.integers(len(LENGTHS)))
            phase = int(rng.integers(src.rmax))
            base, _ = src.base_query(li)
            means = np.roll(base, phase * (src.q // src.rmax))
            normalization_residuals.append(float(means.mean() - 1))
            probs = means / means.sum()
            proposal_types.append([li, phase])
        probability_sum_residuals.append(float(probs.sum() - 1))
        counts[row] = rng.multinomial(src.n, probs)
    assert np.all(counts.sum(axis=1) == src.n)
    labels = np.zeros((len(LENGTHS), src.rmax), dtype=int)
    for li, phase in proposal_types:
        if li >= 0:
            labels[li, phase] += 1
    uniform_rows = sum(li < 0 for li, _ in proposal_types)
    tail = sum(math.comb(samples, k) for k in
               range(min(uniform_rows, samples-uniform_rows) + 1)) / (2 ** samples)
    audit = {"uniform_rows": uniform_rows, "tilt_rows": samples-uniform_rows,
             "label_counts_L_by_phase": labels.tolist(),
             "uniform_branch_exact_two_sided_binomial_p":
                 min(1., 2 * tail),
             "max_abs_probability_sum_residual":
                 max(abs(x) for x in probability_sum_residuals),
             "max_abs_sampled_one_coordinate_tilt_normalizer_minus_one":
                 max([abs(x) for x in normalization_residuals] + [0.]),
             "weight_mean_pointwise_95pct_Hoeffding_radius":
                 2 * math.sqrt(math.log(40) / (2 * samples))}
    return counts, audit


def variance_mix(log_means, relative_variances):
    weights = np.exp(log_means - lse(log_means)[:, None])
    return np.sum(weights * np.sqrt(relative_variances), axis=1) ** 2


def tail_lower(log_mean, relative_variance, weights):
    out = []
    for threshold in THRESHOLDS:
        gap = np.maximum(0, -np.expm1(math.log(threshold) - log_mean))
        factor = np.divide(gap * gap, relative_variance + gap * gap,
                           out=np.zeros_like(gap), where=gap > 0)
        z = weights * factor
        total, sq = z.sum(), np.dot(z, z)
        out.append({"threshold": float(threshold),
                    "probability_lower_MC": float(z.mean()),
                    "coefficient_lower_MC": float(threshold * z.mean()),
                    "coefficient_se_descriptive":
                        float(threshold * z.std(ddof=1) / math.sqrt(len(z))),
                    "maximum_observed_coefficient_contribution":
                        float(threshold * z.max()),
                    "contribution_ess": float(total * total / sq) if sq else 0,
                    "counts_active": int(np.count_nonzero(z)),
                    "mean_conditional_factor": float(factor.mean())})
    return out


def proxy_scan(log_values, weights):
    order = np.argsort(log_values)[::-1]
    cumul = np.cumsum(weights[order])
    score = np.exp(log_values[order]) * cumul / len(weights)
    sq = np.cumsum(weights[order] ** 2)
    ess = cumul ** 2 / sq
    selected = int(np.argmax(np.where(ess >= 16, score, -np.inf)))
    return {"coefficient_proxy_only": float(score[selected]),
            "threshold": float(np.exp(log_values[order[selected]])),
            "tail_ess": float(ess[selected])}


def bounded_chernoff_lower(observed_mean, upper, samples, alpha):
    """One-sided lower confidence via Bernoulli KL mgf domination, no normal fit."""
    if observed_mean <= 0:
        return 0.
    x = min(1., observed_mean / upper)
    target = math.log(1 / alpha) / samples
    def divergence(mu):
        first = x * math.log(x / mu) if x else 0.
        second = ((1-x) * math.log((1-x)/(1-mu))) if x < 1 else 0.
        return first + second
    low, high = 0., x
    for _ in range(70):
        mid = (low + high) / 2
        if divergence(mid) > target:
            low = mid
        else:
            high = mid
    return upper * low


def bin_interface_check():
    """Implementation check only: physical smooth Y vs analytic bin means."""
    rng = np.random.default_rng(SEED + 200)
    nodes, quadrature_weights = np.polynomial.legendre.leggauss(8)
    output = []
    for d, poly, v, rmax in [(8, 0x11d, 4, 8), (10, 0x409, 16, 32)]:
        src = Source(d, poly, v, rmax)
        atoms = src.centers2 / (2 * src.q)
        def primitive(t):
            f = t % 1
            done = np.searchsorted(atoms + src.delta/2, f, side="right")
            candidate = np.minimum(done, src.k-1)
            part = np.clip((f-atoms[candidate]+src.delta/2)/src.delta, 0, 1)
            part[done == src.k] = 0
            return done+part-src.k*f
        errors, ratios = [], []
        for li in [0, 4, 8]:
            base, trans = src.base_query(li)
            transition_bins = np.flatnonzero(trans)
            selected = np.concatenate([rng.integers(src.q, size=20),
                transition_bins[rng.integers(len(transition_bins), size=20)]])
            x = ((selected[:, None] + .5) / src.q
                 + nodes[None, :] / (2*src.q))
            length = LENGTHS[li]
            h = primitive(x+length/2)-primitive(x-length/2)
            freq = np.arange(1,src.m+1)
            b = np.sum(src.c*np.sin(np.pi*freq*length)
                       *np.cos(2*np.pi*x[:, :, None]*freq), axis=2)
            y = 1+src.eps/length*(h+b)
            mean = y @ (quadrature_weights/2)
            variance = ((y-mean[:,None])**2) @ (quadrature_weights/2)
            upper = base[selected]**2*np.where(trans[selected],src.t1,src.t0)
            errors.extend(np.abs(mean-base[selected]).tolist())
            ratios.extend((variance/upper).tolist())
        output.append({"v":v,"K":src.k,"queries_checked":120,
                       "max_abs_bin_mean_difference":max(errors),
                       "max_numerical_variance_over_analytic_upper":max(ratios)})
    return {"seed":SEED+200,"note":"8-point quadrature checks implementation; "
            "the analytic mean and range proof do not rely on quadrature.",
            "results":output}


def merge_results():
    filenames = ["case_v4.json", "case_v8.json", "case_v16.json",
                 "case_v16_diagnostic_seed.json"]
    merged = {"scope": "conditional multinomial interface / synchronization; "
              "not matched sqrt(n) lower-bound sources",
              "seeds": {"main_base":SEED,"diagnostic":SEED+100,
                        "bin_check":SEED+200},
              "reproduction": {"script":Path(__file__).name,"python":sys.executable,
                               "python_version":sys.version,
                               "numpy_version":np.__version__,
                               "no_Gaussian_surrogate":True,
                               "counts_or_large_tables_archived":False},
              "fixed_log_thresholds": [0,1,2,4],
              "fixed_columns":["lambda","coefficient_MC","SE_descriptive",
                               "contribution_ESS","active_counts","maximum_observed_coefficient_contribution"],
              "holdout_columns":["lambda","coefficient_MC","SE_descriptive",
                                 "contribution_ESS","Chernoff_lower",
                                 "integrand_bound"],
              "holdout_alpha_each":.05/51, "cases":[]}
    initial_audit = json.loads((OUT/"v16_original_proposal_audit.json").read_text())
    for filename in filenames:
        full = json.loads((OUT/filename).read_text())
        audit = full.get("proposal_audit", {})
        if filename == "case_v16.json" and not audit:
            audit = initial_audit
        if audit:
            sparse = [[li,ri,count] for li,row in enumerate(
                audit["label_counts_L_by_phase"]) for ri,count in enumerate(row)
                if count]
            audit = {key:value for key,value in audit.items()
                     if key!="label_counts_L_by_phase"}
            audit["nonzero_label_counts_L_phase_count"] = sparse
        else:
            audit = {"initial_branch_counts":"not logged; no extra rerun"}
        radius = 2*math.sqrt(math.log(40)/(2*full["samples"]))
        row = {"local_full_file_not_for_archive":filename,
               "metadata":full["metadata"],"seed":full["seed"],
               "samples":full["samples"],"runtime_seconds":full["runtime_seconds"],
               "mass":{"mean":full["weight_mean"],"descriptive_se":full["weight_mean_se"],
                       "weight_max":full["weight_max"],
                       "pointwise_95pct_Hoeffding_radius":radius,
                       "interval":[full["weight_mean"]-radius,
                                   full["weight_mean"]+radius]},
               "proposal_audit":audit,
               "one_coordinate_normalizer_max_mass_residual":
                   full["max_one_coordinate_mean_mass_residual"],
               "max_n_log_normalizer_sum_order_discrepancy":
                   full.get("max_n_log_normalizer_sum_order_discrepancy"),
               "known_normalizer_contribution_inflation":full["known_normalizer_contribution_inflation"],
               "R_results":{}}
        for r,z in full["R_results"].items():
            fixed = {}
            for target,key in [("V","V_true_tail_lower_by_Cantelli"),
                               ("U","U_true_tail_lower_by_Cantelli"),
                               ("independent",
                                "independent_true_tail_lower_by_Cantelli")]:
                fixed[target] = []
                for log_t in [0,1,2,4]:
                    entry = z[key][int((log_t+2)*4)]
                    fixed[target].append([entry[k] for k in [
                        "threshold","coefficient_lower_MC",
                        "coefficient_se_descriptive","contribution_ess","counts_active","maximum_observed_coefficient_contribution"]])
            holdout = {}
            for target,entry in z["discovery_holdout_Cantelli_lower"].items():
                entry = entry["holdout"]
                holdout[target] = [entry[k] for k in [
                    "threshold","coefficient_lower_MC",
                    "coefficient_se_descriptive","contribution_ess",
                    "coefficient_lower_confidence_Chernoff",
                    "MC_integrand_upper_bound"]]
            row["R_results"][r] = {
                "fixed_threshold_Cantelli_integral_MC":fixed,
                "holdout":holdout,
                "V_proxy":z["V_conditional_mean_proxy"],
                "U_proxy":z["U_conditional_mean_proxy"],
                "independent_proxy":z["independent_conditional_mean_proxy"],
                "proposal_sample_log_Ubar_over_Vbar_quantiles":
                    z["conditional_mean_log_U_over_V_quantiles"],
                "V_relative_variance_upper_quantiles":
                    z["V_relative_variance_upper_quantiles"],
                "U_relative_variance_upper_quantiles":
                    z["U_relative_variance_upper_quantiles"],
                "V_U_coefficient_integrand_upper":
                    z["V_U_IS_integrand_upper"],
                "per_threshold_max_contribution_initially_saved":
                    "maximum_observed_coefficient_contribution" in
                    z["V_true_tail_lower_by_Cantelli"][0]}
        merged["cases"].append(row)
    merged["bin_interface_checks"] = json.loads(
        (OUT/"bin_interface_checks.json").read_text())
    merged["numerical_caution"] = (
        "This execution uses harmonized 1D normalizers in sampling and RN density. "
        "The known raw-mean normalizer exponent is retained in contribution bounds. "
        "The initial-seed v16 branch anomaly and one fixed independent diagnostic "
        "are both retained. No interval floating-point certification is provided.")
    (OUT/"large_dimension_results.json").write_text(json.dumps(merged,separators=(",",":")))
    return merged


def run_case(degree, polynomial, v, rmax, samples, seed, output_suffix=""):
    start = time.time()
    src = Source(degree, polynomial, v, rmax)
    rng = np.random.default_rng(seed)
    counts, proposal_audit = draw_counts(src, rng, samples)
    log_a = np.empty((samples, len(LENGTHS), rmax))
    rel_a = np.empty_like(log_a)
    log_ind = {r: np.empty((samples, len(LENGTHS))) for r in src.r_values}
    rel_ind = {r: np.empty((samples, len(LENGTHS))) for r in src.r_values}
    norms = []
    normalizer_order_discrepancies = []
    for li in range(len(LENGTHS)):
        means, trans = src.query(li)
        old_mean_norm = means.mean(axis=0)
        mean_norm = np.array([np.ascontiguousarray(means[:, r]).mean()
                              for r in range(rmax)])
        normalizer_order_discrepancies.append(np.max(np.abs(
            np.log(old_mean_norm) - np.log(mean_norm))) * src.n)
        norms.append(mean_norm)
        log_a[:, li] = counts @ np.log(means)
        # Binary transition counts are integer sums < 2^53.
        nt = counts @ trans.astype(float)
        exponent = (src.n * math.log1p(src.t0)
                    + nt * (math.log1p(src.t1) - math.log1p(src.t0)))
        rel_a[:, li] = np.expm1(exponent)
        for r in src.r_values:
            indices = np.arange(r) * (rmax // r)
            means_ind = means[:, indices].mean(axis=1)
            log_ind[r][:, li] = counts @ np.log(means_ind)
            transition_fraction = trans[:, indices].mean(axis=1)
            per_bin_bound = (src.eps / (2 * src.q)) ** 2 * (
                src.db + transition_fraction / src.delta) ** 2
            rel_ind[r][:, li] = np.expm1(counts @ np.log1p(per_bin_bound))
        del means, trans
    norms = np.array(norms)
    normalizer_inflation = math.exp(max(0., float(src.n * np.log(norms).max())))
    normalized_log = log_a - src.n * np.log(norms)[None, :, :]
    log_tilt = lse(normalized_log.reshape(samples, -1)) - math.log(
        len(LENGTHS) * rmax)
    log_q = np.logaddexp(math.log(.5), math.log(.5) + log_tilt)
    weights = np.exp(-log_q)
    assert weights.max() <= 2 + 1e-12
    result = {"metadata": src.metadata(), "seed": seed, "samples": samples,
              "numpy_version": np.__version__, "python": sys.executable,
              "proposal": "half uniform multinomial; half all L,phase bin tilts",
              "weight_mean": float(weights.mean()),
              "weight_mean_se": float(weights.std(ddof=1) / math.sqrt(samples)),
              "weight_max": float(weights.max()),
              "proposal_audit": proposal_audit,
              "known_normalizer_contribution_inflation": normalizer_inflation,
              "max_n_log_normalizer_sum_order_discrepancy":
                  float(max(normalizer_order_discrepancies)),
              "max_one_coordinate_mean_mass_residual":
                  float(np.max(np.abs(norms - 1))),
              "R_results": {}}
    rows = np.arange(samples)
    for r in src.r_values:
        indices = np.arange(r) * (rmax // r)
        z = log_a[:, :, indices]
        rv = rel_a[:, :, indices]
        log_v_all = lse(z, axis=2) - math.log(r)
        best_l = np.argmax(log_v_all, axis=1)
        log_v = log_v_all[rows, best_l]
        rv_v = variance_mix(z[rows, best_l], rv[rows, best_l])
        best_per_phase = np.argmax(z, axis=1)
        selected_z = np.take_along_axis(z, best_per_phase[:, None, :], axis=1)[:, 0]
        selected_rv = np.take_along_axis(rv, best_per_phase[:, None, :], axis=1)[:, 0]
        log_u = lse(selected_z) - math.log(r)
        rv_u = variance_mix(selected_z, selected_rv)
        assert np.min(log_u - log_v) >= -1e-10
        best_ind = np.argmax(log_ind[r], axis=1)
        zi = log_ind[r][rows, best_ind]
        ri = rel_ind[r][rows, best_ind]
        lower_v = tail_lower(log_v, rv_v, weights)
        lower_u = tail_lower(log_u, rv_u, weights)
        lower_i = tail_lower(zi, ri, weights)
        cut = samples // 2
        holdout = {}
        for name, zz, rr in [("V", log_v, rv_v), ("U", log_u, rv_u),
                              ("independent", zi, ri)]:
            discovery = tail_lower(zz[:cut], rr[:cut], weights[:cut])
            chosen = int(np.argmax([x["coefficient_lower_MC"] for x in discovery]))
            hold = tail_lower(zz[cut:], rr[cut:], weights[cut:])[chosen]
            bound = 2 * hold["threshold"]
            if name != "independent":
                bound = min(bound, 2 * len(LENGTHS) * rmax / r * normalizer_inflation)
            hold["MC_integrand_upper_bound"] = bound
            hold["one_sided_MC_alpha"] = .05 / 51
            hold["coefficient_lower_confidence_Chernoff"] = bounded_chernoff_lower(
                hold["coefficient_lower_MC"], bound, samples-cut, .05 / 51)
            holdout[name] = {"discovery_threshold": discovery[chosen]["threshold"],
                             "holdout": hold}
        result["R_results"][str(r)] = {
            "V_conditional_mean_proxy": proxy_scan(log_v, weights),
            "U_conditional_mean_proxy": proxy_scan(log_u, weights),
            "independent_conditional_mean_proxy": proxy_scan(zi, weights),
            "conditional_mean_log_U_over_V_quantiles":
                np.quantile(log_u-log_v, [0, .5, .9, .99, 1]).tolist(),
            "V_relative_variance_upper_quantiles":
                np.quantile(rv_v, [0, .5, .9, 1]).tolist(),
            "U_relative_variance_upper_quantiles":
                np.quantile(rv_u, [0, .5, .9, 1]).tolist(),
            "V_true_tail_lower_by_Cantelli": lower_v,
            "U_true_tail_lower_by_Cantelli": lower_u,
            "independent_true_tail_lower_by_Cantelli": lower_i,
            "discovery_holdout_Cantelli_lower": holdout,
            "V_U_IS_integrand_upper": 2 * len(LENGTHS) * rmax / r * normalizer_inflation,
        }
    result["runtime_seconds"] = time.time() - start
    (OUT / f"case_v{v}{output_suffix}.json").write_text(json.dumps(result, indent=2))
    # Small query summaries suffice for archival; no counts / large tables saved.
    print(json.dumps({"v": v, "K": src.k, "n": src.n, "Q": src.q,
                      "seconds": result["runtime_seconds"],
                      "mass": result["weight_mean"],
                      "Rmax_holdout": result["R_results"][str(rmax)]
                          ["discovery_holdout_Cantelli_lower"]}), flush=True)
    return result


def main():
    configs = [(8, 0x11d, 4, 8, 256),
               (10, 0x409, 8, 16, 192),
               (10, 0x409, 16, 32, 128)]
    summaries = []
    for j, (d, poly, v, r, samples) in enumerate(configs):
        res = run_case(d, poly, v, r, samples, SEED + j)
        if v == 16:
            (OUT / "v16_original_proposal_audit.json").write_text(
                json.dumps(res["proposal_audit"], indent=2))
        summaries.append({"v": v, "metadata": res["metadata"],
                          "seconds": res["runtime_seconds"],
                          "mass": res["weight_mean"], "mass_se": res["weight_mean_se"]})
    (OUT / "summary.json").write_text(json.dumps(summaries, indent=2))
    # Fixed-size independent diagnostic, retained alongside the anomalous seed.
    run_case(10, 0x409, 16, 32, 128, SEED + 100,
             output_suffix="_diagnostic_seed")
    (OUT / "bin_interface_checks.json").write_text(
        json.dumps(bin_interface_check(), indent=2))
    merge_results()


if __name__ == "__main__":
    main()
