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
SEED = 20261029033
BASELINE_SHA256 = 'd9efad5873ebd3c5e667351f740bd716d5b2d91106e82ce436399bc3c819b1e1'
P = .5
BN = 4
LOG_THRESHOLDS = np.arange(-2, 10.00001, .25)
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
        assert self.q == 24*self.k*v
        self.lengths = 10/9 + np.arange(self.k*v//9+1)/(self.k*v)
        self.length_q = 80*(self.k*v//3)+24*np.arange(len(self.lengths))
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
        lq = int(self.length_q[li])
        upper = (self.mid2 + lq) % (2 * self.q)
        lower = (self.mid2 - lq) % (2 * self.q)
        pu, iu = self.primitive_numerator(upper)
        pl, il = self.primitive_numerator(lower)
        h = (pu - pl) / (2 * self.q)
        freq = np.arange(1, self.m + 1, dtype=float)
        # Mean cosine in a bin: center cosine * sinc(pi*m/Q).
        coefficient = (self.c * np.sin(np.pi * freq * self.lengths[li])
                       * np.sinc(freq / self.q))
        bmean = self.cosines @ coefficient
        means = 1 + self.eps / self.lengths[li] * (h + bmean)
        assert means.min() >= 1 / self.lengths[li] - 1e-12
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
                "conditional_variance_DB_upper": self.db,
                "natural_scale_count":len(self.lengths),
                "natural_scale_first":str(Fraction(10,9)),
                "natural_scale_last":str(Fraction(10,9)+Fraction(self.k*self.v//9,self.k*self.v)),
                "gap_to_11_over_9":str(Fraction(11,9)-Fraction(10,9)-Fraction(self.k*self.v//9,self.k*self.v))}


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
            li = int(rng.integers(len(src.lengths)))
            phase = int(rng.integers(src.rmax))
            base, _ = src.base_query(li)
            means = np.roll(base, phase * (src.q // src.rmax))
            normalization_residuals.append(float(means.mean() - 1))
            probs = means / means.sum()
            proposal_types.append([li, phase])
        probability_sum_residuals.append(float(probs.sum() - 1))
        counts[row] = rng.multinomial(src.n, probs)
    assert np.all(counts.sum(axis=1) == src.n)
    labels = np.zeros((len(src.lengths), src.rmax), dtype=int)
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


def cantelli(log_mean, relative_variance, threshold):
    gap = np.maximum(0, -np.expm1(math.log(threshold)-log_mean))
    return np.divide(gap*gap, relative_variance+gap*gap,
                     out=np.zeros_like(gap), where=gap>0)


def product_upper(log_means, relative_variances, threshold):
    # Upper-tail contributions may be positive when all means <= threshold.
    # Therefore only the 2*lambda weight bound applies to this object.
    gap = np.maximum(0, np.expm1(math.log(threshold)-log_means))
    each = np.divide(relative_variances, relative_variances+gap*gap,
                     out=np.ones_like(gap), where=gap>0)
    markov = np.exp(np.minimum(0,log_means-math.log(threshold)))
    return np.minimum(1, np.minimum(each,markov).sum(axis=1))


def statistics(z, upper=None, alpha=.05/100):
    total, squares = z.sum(), np.dot(z,z)
    result={"mean":float(z.mean()),
            "SE_descriptive":float(z.std(ddof=1)/math.sqrt(len(z))),
            "ESS":float(total*total/squares) if squares else 0.,
            "active":int(np.count_nonzero(z)),
            "max_observed":float(z.max())}
    if np.any(z<0):
        result["signed_cancellation_ratio_not_ESS"]=result.pop("ESS")
    if upper is not None:
        assert z.min()>=-1e-10 and z.max()<=upper*(1+1e-9)+1e-10
        result["contribution_upper"]=float(upper)
        result["Chernoff_lower"]=bounded_chernoff_lower(
            result["mean"],upper,len(z),alpha)
    return result


def run_natural(v,rmax,samples,seed):
    started=time.time()
    src=Source(8,0x11d,v,rmax)
    jcount=len(src.lengths)
    assert jcount==src.k*v//9+1
    rng=np.random.default_rng(seed)
    counts,audit=draw_counts(src,rng,samples)
    log_a=np.empty((samples,jcount,rmax))
    rel_a=np.empty_like(log_a)
    log_i={r:np.empty((samples,jcount)) for r in src.r_values}
    rel_i={r:np.empty((samples,jcount)) for r in src.r_values}
    norms=[];sum_order=[]
    for li in range(jcount):
        means,transition=src.query(li)
        norm=np.array([np.ascontiguousarray(means[:,r]).mean()
                       for r in range(rmax)])
        norms.append(norm)
        sum_order.append(src.n*np.max(np.abs(np.log(means.mean(axis=0))
                                             -np.log(norm))))
        log_a[:,li]=counts@np.log(means)
        nt=counts@transition.astype(float)
        rel_a[:,li]=np.expm1(src.n*math.log1p(src.t0)
            +nt*(math.log1p(src.t1)-math.log1p(src.t0)))
        for r in src.r_values:
            ix=np.arange(r)*(rmax//r)
            log_i[r][:,li]=counts@np.log(means[:,ix].mean(axis=1))
            fraction=transition[:,ix].mean(axis=1)
            bound=(src.eps/(2*src.q))**2*(src.db+fraction/src.delta)**2
            rel_i[r][:,li]=np.expm1(counts@np.log1p(bound))
        if li in [0,jcount//2,jcount-1]:
            print(json.dumps({"v":v,"completed_L":li+1,"J":jcount,
                              "seconds":time.time()-started}),flush=True)
    norms=np.array(norms)
    c_z=math.exp(max(0.,float(src.n*np.log(norms).max())))
    normalized=log_a-src.n*np.log(norms)[None,:,:]
    log_d=np.logaddexp(math.log(.5),math.log(.5)
              +lse(normalized.reshape(samples,-1))-math.log(jcount*rmax))
    weights=np.exp(-log_d)
    assert weights.max()<=2+1e-12
    cut=samples//2
    rows=np.arange(samples)
    result={"metadata":src.metadata(),"seed":seed,"N":samples,
            "python":sys.executable,"numpy":np.__version__,
            "baseline_SHA256":BASELINE_SHA256,
            "proposal":"half uniform; half ALL natural L x global phase bin tilts",
            "mass_mean":float(weights.mean()),
            "mass_SE_descriptive":float(weights.std(ddof=1)/math.sqrt(samples)),
            "weight_max":float(weights.max()),"C_Z":c_z,
            "max_n_log_normalizer_sum_order_discrepancy":float(max(sum_order)),
            "max_raw_bin_mean_mass_residual":float(np.abs(norms-1).max()),
            "proposal_audit":audit,
            "conditional_mean_V_rule":"max L of conditional mean",
            "real_V_tail_lower_rule":"max L of Cantelli at EACH threshold",
            "fixed_log_thresholds":[0,.5,1,1.5,2,3,4,5,6],
            "R_results":{}}
    prod_z=log_a[:,:,0];prod_rv=rel_a[:,:,0]
    fixed_log=[0,.5,1,1.5,2,3,4,5,6]
    fixed_t=np.exp(fixed_log)
    lower_product={float(t):cantelli(prod_z,prod_rv,float(t)).max(axis=1)
                   for t in fixed_t}
    upper_product={float(t):product_upper(prod_z,prod_rv,float(t))
                   for t in fixed_t}
    for r in src.r_values:
        ix=np.arange(r)*(rmax//r)
        z=log_a[:,:,ix];rv=rel_a[:,:,ix]
        lv=lse(z,axis=2)-math.log(r)
        phase_weights=np.exp(z-lse(z,axis=2)[:,:,None])
        rvv=np.sum(phase_weights*np.sqrt(rv),axis=2)**2
        best=np.argmax(lv,axis=1);log_v=lv[rows,best]
        best_phase=np.argmax(z,axis=1)
        selected=np.take_along_axis(z,best_phase[:,None,:],axis=1)[:,0]
        selected_rv=np.take_along_axis(rv,best_phase[:,None,:],axis=1)[:,0]
        lu=lse(selected)-math.log(r)
        rvu=variance_mix(selected,selected_rv)
        assert (lu-log_v).min()>-1e-10
        cap=2*jcount*rmax/r*c_z
        def factors(name,t):
            if name=="V":
                return cantelli(lv,rvv,t).max(axis=1)
            if name=="U":
                return cantelli(lu,rvu,t)
            return cantelli(log_i[r],rel_i[r],t).max(axis=1)
        fixed={}
        for t in fixed_t:
            t=float(t);bp=lower_product[t];up=upper_product[t]
            table={}
            for name in ["V","U","independent"]:
                b=factors(name,t);bound=2*t if name=="independent" else min(2*t,cap)
                table[name]=statistics(t*weights*b,bound)
            lower_difference=t*weights*(factors("V",t)-bp)
            strict_difference=t*weights*(factors("V",t)-up)
            table["V_lower_minus_product_lower_descriptive"]=statistics(
                lower_difference)
            table["V_lower_minus_product_upper"]=statistics(strict_difference)
            table["product_upper"]=statistics(t*weights*up,2*t)
            table["product_upper"].pop("Chernoff_lower",None)
            table["product_upper"]["role"]=(
                "conditional upper integral MC; no statistical upper CI")
            table["strict_difference_one_sided_Hoeffding_lower"]=float(
                strict_difference.mean()
                -4*t*math.sqrt(math.log(100/.05)/(2*samples)))
            fixed[str(t)]=table
        holdout={}
        for name in ["V","U","independent"]:
            discovery=[]
            for t in THRESHOLDS:
                b=factors(name,float(t))
                discovery.append(float(t*np.mean(weights[:cut]*b[:cut])))
            chosen=int(np.argmax(discovery));t=float(THRESHOLDS[chosen])
            b=factors(name,t)
            bound=2*t if name=="independent" else min(2*t,cap)
            stat=statistics(t*weights[cut:]*b[cut:],bound)
            stat["threshold"]=t
            stat["discovery_coefficient"]=discovery[chosen]
            stat["discovery_at_scan_boundary"]=chosen in [0,len(THRESHOLDS)-1]
            if name=="V":
                pu=product_upper(prod_z[cut:],prod_rv[cut:],t)
                delta=t*weights[cut:]*(b[cut:]-pu)
                stat["V_lower_minus_product_upper_holdout"]=statistics(delta)
            holdout[name]=stat
        # Store only a few common curve points, not all threshold tables.
        result["R_results"][str(r)]={
            "V_mean_proxy":proxy_scan(log_v,weights),
            "U_mean_proxy":proxy_scan(lu,weights),
            "proposal_log_Ubar_over_Vbar_quantiles":
                np.quantile(lu-log_v,[0,.5,.9,.99,1]).tolist(),
            "V_relative_variance_cap_quantiles":
                np.quantile(rvv[rows,best],[0,.5,.9,1]).tolist(),
            "fixed_thresholds":fixed,"holdout":holdout,
            "C_Z_corrected_lower_IS_cap":cap}
    sparse=[[li,phase,count]for li,row in enumerate(
        result["proposal_audit"].pop("label_counts_L_by_phase"))
        for phase,count in enumerate(row)if count]
    result["proposal_audit"]["nonzero_label_counts_L_phase_count"]=sparse
    result["runtime_seconds"]=time.time()-started
    result["theoretical_counts_plus_single_query_memory_bytes"]=(
        counts.nbytes+src.q*rmax*8*3+log_a.nbytes+rel_a.nbytes)
    (OUT/f"natural_v{v}.json").write_text(json.dumps(result,separators=(",",":")))
    print(json.dumps({"DONE_v":v,"seconds":result["runtime_seconds"],
                      "mass":result["mass_mean"],"Rmax_holdout":
                      result["R_results"][str(rmax)]["holdout"]}),flush=True)
    return result


def main():
    first=run_natural(4,8,256,SEED)
    results=[first]
    # Deterministic cost gate, unrelated to whether the first result is positive.
    if first["runtime_seconds"]<=60 and (
            first["theoretical_counts_plus_single_query_memory_bytes"]<500_000_000):
        results.append(run_natural(8,16,256,SEED+1))
    compact={"scope":"natural finite grid true-tail conditional lower estimates; "
                      "legal positive sources, not original matched lower-bound family",
             "reproduction":{"script":"natural_probe.py","seed_base":SEED,
                             "python":sys.executable,"numpy":np.__version__,
                             "baseline_SHA256":BASELINE_SHA256},
             "cases":results}
    (OUT/"RESULTS.json").write_text(json.dumps(compact,separators=(",",":")))


if __name__=="__main__":
    main()

