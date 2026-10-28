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
SEED = 202610290340
ARCHIVE33_SHA256 = '3a0ba8bbe15112ab68e6841800ac1e2a4b6355d56198ef328bb46e865b292aed'
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


GROUPS=[1,2,4,8]
R=4
ALPHA=.05/1000
FIXED_LOG=[0,1,2,3,4,5,6]


def configure_source(v):
    src=Source(8,0x11d,v,R)
    natural={Fraction(10,9)+Fraction(j,src.k*v)
             for j in range(src.k*v//9+1)}
    old={Fraction(80+j,72)for j in range(9)}
    grid=sorted(natural|old)
    src.lengths=np.array([float(x)for x in grid])
    shifts=[x*src.q for x in grid]
    assert all(x.denominator==1 for x in shifts)
    src.length_q=np.array([int(x)for x in shifts],dtype=np.int64)
    assert src.n%8==0 and src.q==24*src.k*v
    assert len(grid)=={8:233,16:460}[v]
    assert grid[0]==Fraction(10,9)and grid[-1]==Fraction(11,9)
    return src,grid


def physical_y(src,x,length):
    atoms=src.centers2/(2*src.q)
    def prim(t):
        f=t%1
        done=np.searchsorted(atoms+src.delta/2,f,side="right")
        cand=np.minimum(done,src.k-1)
        part=np.clip((f-atoms[cand]+src.delta/2)/src.delta,0,1)
        part[done==src.k]=0
        return done+part-src.k*f
    h=prim(x+length/2)-prim(x-length/2)
    freq=np.arange(1,src.m+1)
    b=np.sum(src.c*np.sin(np.pi*freq*length)
             *np.cos(2*np.pi*x[...,None]*freq),axis=-1)
    return 1+src.eps/length*(h+b)


def direct_models(y):
    # y.shape=(sample, coordinate, phase), one COMMON length already fixed.
    n=y.shape[1];result={}
    for g in GROUPS:
        phase_product=np.prod(y.reshape(len(y),g,n//g,R),axis=2)
        result[g]=np.prod(phase_product.mean(axis=2),axis=1)
    return result


def sample_rho(src,rng,shape):
    size=math.prod(shape)
    values=np.empty(size)
    atom=rng.random(size)<P
    num=int(atom.sum())
    centers=src.centers2/(2*src.q)
    values[atom]=(centers[rng.integers(src.k,size=num)]
                  +src.delta*(rng.random(num)-.5))
    remaining=np.flatnonzero(~atom)
    freq=np.arange(1,src.m+1)
    amp=src.eps*np.pi*freq*src.c/(1-P)
    envelope=1+np.sum(np.abs(amp))
    while len(remaining):
        candidate=rng.random(len(remaining))
        density=1+np.cos(2*np.pi*candidate[:,None]*freq)@amp
        accepted=rng.random(len(remaining))<density/envelope
        values[remaining[accepted]]=candidate[accepted]
        remaining=remaining[~accepted]
    return values.reshape(shape)


def small_control():
    import itertools
    src,grid=configure_source(8)
    rng=np.random.default_rng(SEED)
    samples=512;n=8
    # Direct small-n point proposal using exact rho convolution sampling.
    x=rng.random((samples,n))
    latent=[]
    for row in range(samples):
        if rng.random()<.5:
            latent.append([-1,-1]);continue
        g=int(rng.choice(GROUPS));li=int(rng.integers(len(grid)))
        phases=np.repeat(rng.integers(R,size=g),n//g)/R
        x[row]=(sample_rho(src,rng,(n,))
                -(rng.random(n)-.5)*src.lengths[li]+phases)%1
        latent.append([g,li])
    point_log_d=[]
    for li,length in enumerate(src.lengths):
        y=physical_y(src,(x[:,:,None]-np.arange(R)[None,None,:]/R)%1,length)
        models=direct_models(y)
        point_log_d.append(np.column_stack([models[g]for g in GROUPS]))
    target=np.stack(point_log_d,axis=1)
    density=.5+.5*target.mean(axis=(1,2))
    weights=1/density
    assert weights.max()<=2
    mass={"mean":float(weights.mean()),
          "SE":float(weights.std(ddof=1)/math.sqrt(samples)),
          "weight_max":float(weights.max())}
    # Independent explicit latent enumeration at one real centered cube query.
    sample_x=np.array([[.03,.14,.27,.41,.53,.66,.79,.93]])
    y=physical_y(src,(sample_x[:,:,None]-np.arange(R)[None,None,:]/R)%1,7/6)
    direct=direct_models(y)
    enumeration={}
    for g in GROUPS:
        block_phase=np.prod(y.reshape(1,g,n//g,R),axis=2)[0]
        codes=np.array(list(itertools.product(range(R),repeat=g)))
        expanded=np.prod(block_phase[np.arange(g)[None,:],codes],axis=1).mean()
        enumeration[str(g)]={"formula":float(direct[g][0]),
                             "explicit_R_to_G_mixture":float(expanded),
                             "absolute_difference":abs(float(expanded-direct[g][0]))}
    # Conditional means/variance versus actual within-bin coordinates.
    midpoint=grid.index(Fraction(7,6))
    _,transition=src.base_query(midpoint)
    transition_bins=np.flatnonzero(transition)
    bins=np.concatenate([transition_bins[rng.integers(len(transition_bins),size=4)],
                          rng.integers(src.q,size=4)])
    conditional=[]
    for li in [0,midpoint,len(grid)-1]:
        means,trans=src.query(li)
        mm=means[bins];tt=np.where(trans[bins],src.t1,src.t0)
        xx=(bins[None,:]+rng.random((4096,n)))/src.q
        yy=physical_y(src,(xx[:,:,None]-np.arange(R)[None,None,:]/R)%1,
                      src.lengths[li])
        observed=direct_models(yy)
        for g in GROUPS:
            block_m=[];block_r=[]
            for block in range(g):
                part=slice(block*n//g,(block+1)*n//g)
                a=np.prod(mm[part],axis=0)
                rr=np.expm1(np.log1p(tt[part]).sum(axis=0))
                block_m.append(a.mean())
                block_r.append((np.sum(a*np.sqrt(rr))/a.sum())**2)
            m=float(np.prod(block_m))
            rv=float(np.expm1(np.log1p(block_r).sum()))
            obs=observed[g]
            se=float(obs.std(ddof=1)/math.sqrt(len(obs)))
            error=float(obs.mean()-m)
            assert abs(error)<=1e-10+6*se
            ratio=float(obs.var(ddof=1)/(m*m*rv))if rv else 0.
            assert ratio<=1.1
            conditional.append({"L":str(grid[li]),"G":g,"analytic_mean":m,
                "direct_within_bin_mean":float(obs.mean()),"mean_error":error,
                "MC_SE":se,"numerical_variance_over_proved_cap":ratio})
    assert max(x["absolute_difference"]for x in enumeration.values())<1e-12
    assert abs(mass["mean"]-1)<1e-6+6*mass["SE"]
    return {"seed":SEED,"n":8,"point_samples":samples,
            "point_proposal":"half uniform + half all G/all target L block mixtures",
            "mass":mass,"latent_enumeration":enumeration,
            "within_bin_controls":conditional,"passed":True}


def draw_fine_counts(src,rng,samples):
    fine_n=src.n//8
    counts=np.empty((samples,8,src.q),dtype=float)
    uniform=np.full(src.q,1/src.q)
    latent=[];residuals=[]
    for row in range(samples):
        if rng.random()<.5:
            for h in range(8):counts[row,h]=rng.multinomial(fine_n,uniform)
            latent.append({"G":0});continue
        g=int(rng.choice(GROUPS));li=int(rng.integers(len(src.lengths)))
        phases=rng.integers(R,size=g)
        means,_=src.base_query(li)
        probabilities={}
        for phase in set(phases.tolist()):
            m=np.roll(means,phase*(src.q//R))
            probs=m/m.sum();probabilities[phase]=probs
            residuals.append(float(probs.sum()-1))
        for coarse,phase in enumerate(phases):
            for h in range(coarse*(8//g),(coarse+1)*(8//g)):
                counts[row,h]=rng.multinomial(fine_n,probabilities[int(phase)])
        latent.append({"G":g,"L_index":li,"coarse_phases":phases.tolist()})
    assert np.all(counts.sum(axis=2)==fine_n)
    uniform_rows=sum(x["G"]==0 for x in latent)
    tail=sum(math.comb(samples,j)for j in
        range(min(uniform_rows,samples-uniform_rows)+1))/(2**samples)
    audit={"uniform_rows":uniform_rows,
           "latent_G_counts":{str(g):sum(x["G"]==g for x in latent)for g in GROUPS},
           "uniform_branch_exact_two_sided_binomial_p":min(1.,2*tail),
           "max_probability_sum_residual":max([abs(x)for x in residuals]+[0.]),
           "latent_rows":latent}
    return counts,audit


def lower_factor(log_m,rv,threshold):
    gap=np.maximum(0,-np.expm1(math.log(threshold)-log_m))
    each=np.divide(gap*gap,rv+gap*gap,out=np.zeros_like(gap),where=gap>0)
    return each.max(axis=1)


def upper_factor(log_m,rv,threshold):
    gap=np.maximum(0,np.expm1(math.log(threshold)-log_m))
    each=np.divide(rv,rv+gap*gap,out=np.ones_like(gap),where=gap>0)
    markov=np.exp(np.minimum(0,log_m-math.log(threshold)))
    return np.minimum(1,np.minimum(each,markov).sum(axis=1))


def bounded_lower(mean,bound,n,alpha=ALPHA):
    if mean<=0:return 0.
    x=min(1.,mean/bound);target=math.log(1/alpha)/n
    lo,hi=0.,x
    for _ in range(70):
        mid=(lo+hi)/2
        kl=x*math.log(x/mid)
        if x<1:kl+=(1-x)*math.log((1-x)/(1-mid))
        if kl>target:lo=mid
        else:hi=mid
    return bound*lo


def stats(z,bound=None,signed=False):
    total=float(z.sum());square=float(np.dot(z,z))
    out={"mean":float(z.mean()),
         "SE_descriptive":float(z.std(ddof=1)/math.sqrt(len(z))),
         "active":int(np.count_nonzero(z)),"min":float(z.min()),
         "max":float(z.max())}
    if signed:
        out["Hoeffding_one_sided_lower"]=float(
            z.mean()-4*bound*math.sqrt(math.log(1/ALPHA)/(2*len(z))))
        out["range_bound"]=[-2*bound,2*bound]
    else:
        out["ESS"]=total*total/square if square else 0.
        if bound is not None:
            assert z.min()>=0 and z.max()<=bound*(1+1e-9)+1e-10
            out["contribution_cap"]=bound
            out["Chernoff_one_sided_lower"]=bounded_lower(out["mean"],bound,len(z))
    return out


def estimate_memory(src,samples):
    j=len(src.lengths)
    # Release counts before forming coarse-mixture working arrays.
    phase=samples*8*j*R*8
    stage1=samples*8*src.q*8+2*phase+12*src.q*R*8+40_000_000
    stage2=8*phase+12*src.q*R*8+40_000_000
    return int(max(stage1,stage2))


def run_case(v,samples,seed):
    started=time.time();src,grid=configure_source(v);jcount=len(grid)
    rng=np.random.default_rng(seed)
    counts,audit=draw_fine_counts(src,rng,samples)
    flat=counts.reshape(samples*8,src.q)
    fine_z=np.empty((samples,8,jcount,R))
    fine_e=np.empty_like(fine_z)
    norms=[];order_discrepancies=[]
    for li in range(jcount):
        means,trans=src.query(li)
        norm=np.array([np.ascontiguousarray(means[:,r]).mean()for r in range(R)])
        norms.append(norm)
        order_discrepancies.append(src.n*np.max(np.abs(
            np.log(means.mean(axis=0))-np.log(norm))))
        fine_z[:,:,li,:]=(flat@np.log(means)).reshape(samples,8,R)
        nt=(flat@trans.astype(float)).reshape(samples,8,R)
        fine_e[:,:,li,:]=((src.n//8)*math.log1p(src.t0)
            +nt*(math.log1p(src.t1)-math.log1p(src.t0)))
        if li in [0,jcount//2,jcount-1]:
            print(json.dumps({"v":v,"completed_L":li+1,"J":jcount,
                              "seconds":time.time()-started}),flush=True)
    # All later queries use stored fine phase logs/exponents, not counts.
    del flat,counts
    norms=np.array(norms);log_norm=np.log(norms)
    c_z=math.exp(max(0.,float(src.n*log_norm.max())))
    model_z={};model_rv={};proposal_z=[]
    for g in GROUPS:
        coarse=fine_z.reshape(samples,g,8//g,jcount,R).sum(axis=2)
        exponent=fine_e.reshape(samples,g,8//g,jcount,R).sum(axis=2)
        phase_rv=np.expm1(exponent)
        log_block=lse(coarse,axis=3)-math.log(R)
        phase_weights=np.exp(coarse-lse(coarse,axis=3)[:,:,:,None])
        rv_block=np.sum(phase_weights*np.sqrt(phase_rv),axis=3)**2
        model_z[g]=log_block.sum(axis=1)
        model_rv[g]=np.expm1(np.log1p(rv_block).sum(axis=1))
        normalized=coarse-(src.n//g)*log_norm[None,None,:,:]
        proposal_z.append((lse(normalized,axis=3)-math.log(R)).sum(axis=1))
        del coarse,exponent,phase_rv,phase_weights,normalized
    log_tilt=lse(np.stack(proposal_z,axis=1).reshape(samples,-1)) \
              -math.log(len(GROUPS)*jcount)
    log_d=np.logaddexp(math.log(.5),math.log(.5)+log_tilt)
    weights=np.exp(-log_d);assert weights.max()<=2
    product_z=fine_z[:,:,:,0].sum(axis=1)
    product_rv=np.expm1(fine_e[:,:,:,0].sum(axis=1))
    model_z["product"]=product_z;model_rv["product"]=product_rv
    cap=2*len(GROUPS)*jcount*c_z
    fixed={};cut=samples//2
    for h in FIXED_LOG:
        t=math.exp(h)
        lower={g:lower_factor(model_z[g],model_rv[g],t)for g in GROUPS}
        upper={g:upper_factor(model_z[g],model_rv[g],t)for g in GROUPS}
        up=upper_factor(product_z,product_rv,t)
        entry={"lambda":t,"G":{},"between_G_lower_minus_upper":{}}
        for g in GROUPS:
            lo=stats(t*weights*lower[g],min(2*t,cap))
            uu=stats(t*weights*upper[g],2*t)
            uu.pop("Chernoff_one_sided_lower",None)
            difference=t*weights*(lower[g]-up)
            entry["G"][str(g)]={"true_tail_lower_integral":lo,
                "conditional_union_upper_integral":uu,
                "lower_minus_phase0_product_upper":
                    stats(difference,bound=t,signed=True)}
            assert np.max(lower[g]-upper[g])<=1e-12
            for other in GROUPS:
                if g!=other:
                    entry["between_G_lower_minus_upper"][f"{g}_minus_{other}"] \
                        =stats(t*weights*(lower[g]-upper[other]),bound=t,signed=True)
        entry["phase0_product_upper"]=stats(t*weights*up,2*t)
        entry["phase0_product_upper"].pop("Chernoff_one_sided_lower",None)
        if h==2:
            entry["primary_paired_D_G_minus_product"]={
                str(g):(t*weights*(lower[g]-up)).tolist()for g in GROUPS}
        fixed[str(h)]=entry
    holdout={}
    for g in GROUPS+["product"]:
        discovery=[]
        for t in THRESHOLDS:
            b=lower_factor(model_z[g],model_rv[g],float(t))
            discovery.append(float(t*np.mean(weights[:cut]*b[:cut])))
        ix=int(np.argmax(discovery));t=float(THRESHOLDS[ix])
        b=lower_factor(model_z[g],model_rv[g],t)
        bound=2*t if g=="product" else min(2*t,cap)
        row=stats(t*weights[cut:]*b[cut:],bound)
        row["lambda"]=t;row["discovery_coefficient"]=discovery[ix]
        row["scan_boundary_selected"]=ix in [0,len(THRESHOLDS)-1]
        if g!="product":
            up=upper_factor(product_z[cut:],product_rv[cut:],t)
            row["lower_minus_product_upper_holdout"]=stats(
                t*weights[cut:]*(b[cut:]-up),bound=t,signed=True)
        holdout[str(g)]=row
    meta=src.metadata()
    meta.update({"J":jcount,"query_first":"10/9","query_last":"11/9",
                 "natural_scale_count":src.k*v//9+1,
                 "union_old9_added":6,"G":GROUPS,"R":R})
    for key in ["natural_scale_first","natural_scale_last","gap_to_11_over_9"]:
        meta.pop(key,None)
    result={"metadata":meta,"seed":seed,"N":samples,"numpy":np.__version__,
        "python":sys.executable,"source_archive33_SHA256":ARCHIVE33_SHA256,
        "mass_mean":float(weights.mean()),
        "mass_SE_descriptive":float(weights.std(ddof=1)/math.sqrt(samples)),
        "weight_max":float(weights.max()),"C_Z":c_z,
        "max_raw_mean_mass_residual":float(np.abs(norms-1).max()),
        "max_n_log_normalizer_order_difference":float(max(order_discrepancies)),
        "proposal_audit":audit,"C_Z_corrected_lower_contribution_cap":cap,
        "fixed_log_thresholds":FIXED_LOG,"fixed":fixed,"holdout":holdout,
        "relative_variance_cap_selected_mean_max":{
            str(g):np.quantile(model_rv[g][np.arange(samples),
                np.argmax(model_z[g],axis=1)],[0,.5,.9,1]).tolist()
                for g in GROUPS},
        "memory_major_arrays_estimate_bytes":estimate_memory(src,samples),
        "runtime_seconds":time.time()-started}
    print(json.dumps({"DONE_v":v,"seconds":result["runtime_seconds"],
                      "mass":result["mass_mean"],"holdout":
                      {g:{"mean":x["mean"],"SE":x["SE_descriptive"],"ESS":x["ESS"]}
                       for g,x in holdout.items()}}),flush=True)
    return result


def main():
    import hashlib
    registration=json.loads((OUT/"registration.json").read_text())
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest() \
        ==registration["script_SHA256"]
    if (OUT/"RESULTS.json").exists():raise RuntimeError("Refusing repeated experiment.")
    control=small_control()
    (OUT/"small_control.json").write_text(json.dumps(control,separators=(",",":")))
    first=run_case(8,128,202610290341)
    cases=[first]
    source2,_=configure_source(16)
    ratio=(source2.q*len(source2.lengths))/(48960*233)
    predicted=first["runtime_seconds"]*ratio
    gate={"first_seconds":first["runtime_seconds"],
          "predicted_second_seconds":predicted,
          "second_memory_estimate":estimate_memory(source2,128)}
    proceed=(first["runtime_seconds"]<60 and predicted<60
             and estimate_memory(source2,128)<1_000_000_000)
    gate["proceed"]=proceed
    if proceed:cases.append(run_case(16,128,202610290342))
    out={"registration_SHA256":hashlib.sha256(
        (OUT/"registration.json").read_bytes()).hexdigest(),
        "numpy":np.__version__,"python":sys.executable,
        "scope":"true common-L block-mixture finite-grid conditional tail evidence;"
                " not original matched lower-bound family",
        "small_control":control,"cost_gate":gate,"cases":cases}
    (OUT/"RESULTS.json").write_text(json.dumps(out,separators=(",",":")))


if __name__=="__main__":
    main()

