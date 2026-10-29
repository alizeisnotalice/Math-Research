"""Common-D quantile-gated source tables. No R/J cache; ordinary FP, not interval arithmetic.

The analytic midpoint error and the empirical floating-point guard are exposed
separately. The latter is not a rigorous rounding certificate.
"""
from fractions import Fraction
import math
import numpy as np

# Copied verbatim from round034_block_mixture/block_probe.py.
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


def trapezoid_cdf(t, width_a, width_b):
    """CDF of sum/difference of centered uniforms; no distant four-square differences."""
    a, b = max(width_a, width_b), min(width_a, width_b)
    t = np.asarray(t, dtype=float)
    d, e = (a+b)/2, (a-b)/2
    out = np.zeros_like(t)
    first = (t > -d) & (t < -e)
    middle = (t >= -e) & (t <= e)
    last = (t > e) & (t < d)
    out[first] = (t[first]+d)**2/(2*a*b)
    out[middle] = .5+t[middle]/a
    out[last] = 1-(d-t[last])**2/(2*a*b)
    out[t >= d] = 1
    return out


class Source:
    def __init__(self, R=4, w=4, v=8):
        if R not in [4,8,16,32] or w != 4 or v != 8:
            raise ValueError("Require R=4/8/16/32, with frozen w4/v8.")
        self.R, self.w, self.v, self.p = R, w, v, .5
        self.K, self.labels, self.field = field_labels(8, 0x11d)
        self.k = self.K
        self.n = 116*v*(self.K//3)**2
        self.Q = 6*w*self.K*v
        self.q = self.Q
        self.delta = 1/(w*self.K*v)
        self.binwidth = 1/self.Q
        self.bins_center = (np.arange(self.Q)+.5)/self.Q
        self.M = math.ceil(math.sqrt(v))
        self.freq = np.arange(1,self.M+1,dtype=float)
        bcoef = 4*(729/29)*((np.sin(np.pi*self.freq*10/27)/(np.pi*self.freq))**4
                  +3*(np.sin(np.pi*self.freq/9)/(np.pi*self.freq))**4)
        self.c = np.sqrt(2*bcoef/3)
        self.cut = self.K//2
        self.a = self.cut/self.K
        self.pchi = self.a+self.p/(2*self.K)*float(
            np.sin(2*np.pi*self.freq*self.a)@self.c)
        self.class_mass = {True:self.pchi, False:1-self.pchi}
        harmonic = sum((Fraction(1,j)for j in range(1,self.M+1)),Fraction(0))
        # Exact rational coefficient bounds, not floating maximum searches.
        self.Dsup_bound = float(Fraction(20,3)*harmonic)
        self.Dprime_bound = float(Fraction(880,21)*self.M)
        self.background_lower = 1-self.p-self.p/self.K*self.Dsup_bound
        assert self.background_lower > 0
        self.atoms = np.empty((R,self.K))
        self.atom_offset_numerators = np.empty((R,self.K),dtype=np.int64)
        indices = np.arange(self.K)
        for table in range(R):
            # u=(label/4+(b+.5)/R) mod1 in plus cells, and
            #   (label/4+(.5-b)/R) mod1 in minus cells.
            # Odd numerator over 2R keeps the wrap and cell-edge cases exact.
            numerator=(self.labels*(R//2)+1
                       +np.where(indices<self.cut,2*table,-2*table))%(2*R)
            self.atom_offset_numerators[table]=numerator
            self.atoms[table]=(indices+numerator/(2*R))/self.K
            # In cell units, half-box is 1/(2wv). These integer comparisons
            # include exact touching at R32, without a fragile FP > test.
            assert np.all(numerator*self.w*self.v>=R)
            assert np.all((2*R-numerator)*self.w*self.v>=R)
        # A numerical guard only. The module never calls it interval certified.
        self.floating_point_guard_abs = 1e-13

    def _class_atoms(self, table, plus):
        return self.atoms[table,:self.cut] if plus else self.atoms[table,self.cut:]

    def _extended_atoms(self, table, plus):
        atoms = self._class_atoms(table,plus)
        return np.concatenate([atoms+s for s in [-1,0,1,2]])

    def _background_extra(self, center, ell, plus):
        """Exact real formula at query center; constant and D parts separate."""
        center = np.asarray(center,dtype=float)
        overlap = np.zeros_like(center)
        d_integral = np.zeros_like(center)
        a, z = (0.,self.a) if plus else (self.a,1.)
        for shift in [-1,0,1,2]:
            lower = np.maximum(center-ell/2,a+shift)
            upper = np.minimum(center+ell/2,z+shift)
            width = np.maximum(upper-lower,0)
            midpoint = (upper+lower)/2
            overlap += width
            for m,c in zip(self.freq,self.c):
                # 1/2 c [sin(2pi*m*upper)-sin(2pi*m*lower)].
                d_integral += c*np.cos(2*np.pi*m*midpoint)*np.sin(np.pi*m*width)
        return overlap, d_integral

    def background_breakpoints(self, ell, plus):
        a,z = (0.,self.a) if plus else (self.a,1.)
        return sorted({edge+shift+side*ell/2 for edge in [a,z]
                       for shift in [-1,0,1,2] for side in [-1,1]})

    def _background_bin_mean(self, centers, ell, plus):
        h = self.binwidth
        overlap, d_integral = self._background_extra(centers,ell,plus)
        # Most bins have no class breakpoint: overlap linear, D-integral C2.
        for point in self.background_breakpoints(ell,plus):
            indices = np.flatnonzero((centers-h/2 < point)&(point < centers+h/2))
            for index in indices:
                left,right = centers[index]-h/2,centers[index]+h/2
                cuts = [left]+[v for v in self.background_breakpoints(ell,plus)
                               if left<v<right]+[right]
                cuts = np.unique(cuts)
                mids=(cuts[:-1]+cuts[1:])/2
                o,d=self._background_extra(mids,ell,plus)
                widths=np.diff(cuts)/h
                overlap[index] = float(o@widths)
                d_integral[index] = float(d@widths)
        return (1-self.p)*overlap+self.p/self.K*d_integral

    def _peak_bin_mean_and_transition_count(self, centers, ell, table, plus):
        atoms=self._extended_atoms(table,plus)
        support=(self.delta+self.binwidth)/2
        if ell <= 2*support:
            raise ValueError("Use L=1 or L-1>delta+binwidth in this macro interface.")
        lo,hi=centers-ell/2,centers+ell/2
        full=(np.searchsorted(atoms,hi-support,side="right")
              -np.searchsorted(atoms,lo+support,side="left"))
        mean=full.astype(float)
        transition=np.zeros(len(centers),dtype=np.int64)
        for edge,is_left in [(lo,True),(hi,False)]:
            transition_count=(np.searchsorted(atoms,edge+support,side="right")
                              -np.searchsorted(atoms,edge-support,side="left"))
            assert transition_count.max()<=1
            transition+=transition_count
            # Exclude the endpoint already counted as a full atom; retain
            # the zero-contribution outer endpoint harmlessly.
            first=np.searchsorted(atoms,edge-support,side="left" if is_left else "right")
            last=np.searchsorted(atoms,edge+support,side="left" if is_left else "right")
            rows=np.flatnonzero(last-first)
            atom=atoms[first[rows]]
            # The other endpoint is outside compact support.
            cdf=trapezoid_cdf(edge[rows]-atom,self.delta,self.binwidth)
            mean[rows]+=1-cdf if is_left else cdf
        return mean,transition

    def evaluate(self,L,b):
        L=float(L);b=int(b)
        if not (0<=b<self.R and 1<=L<2):
            raise ValueError("Require 0<=b<R and 1<=L<2.")
        if L==1:
            plus=np.full(self.Q,self.pchi);minus=np.full(self.Q,1-self.pchi)
            zero=np.zeros(self.Q)
            return {"mean_plus":plus,"mean_minus":minus,
                "relative_variance_class_max":zero,"mean_error_abs":self.floating_point_guard_abs,
                "analytic_mean_error_abs":0.,"floating_point_guard_abs":self.floating_point_guard_abs,
                "normalizer":1.,"range_plus":zero,"range_minus":zero,
                "transition_count_plus":zero.astype(int),"transition_count_minus":zero.astype(int)}
        ell=L-1
        centers=self.bins_center+.5
        means={};ranges={};counts={}
        analytic_error=self.p/(self.K*L)*self.Dprime_bound*self.binwidth**2/12
        error=analytic_error+self.floating_point_guard_abs
        relatives=[]
        for plus in [True,False]:
            bg=self._background_bin_mean(centers,ell,plus)
            peak,transition=self._peak_bin_mean_and_transition_count(centers,ell,b,plus)
            mean=(self.class_mass[plus]+bg+self.p/self.K*peak)/L
            width=((1-self.p)*self.binwidth
                   +self.p/self.K*transition*min(1,self.binwidth/self.delta)
                   +self.p/self.K*2*self.Dsup_bound*self.binwidth)/L
            lower=mean-error
            assert lower.min()>0
            relative=(width/(2*lower))**2
            means[plus]=mean;ranges[plus]=width;counts[plus]=transition
            relatives.append(relative)
        return {"mean_plus":means[True],"mean_minus":means[False],
                "relative_variance_class_max":np.maximum(*relatives),
                "mean_error_abs":error,"analytic_mean_error_abs":analytic_error,
                "floating_point_guard_abs":self.floating_point_guard_abs,
                "normalizer":float(np.mean(means[True]+means[False])),
                "range_plus":ranges[True],"range_minus":ranges[False],
                "transition_count_plus":counts[True],"transition_count_minus":counts[False]}

    def point_class_kernel(self,x,L,b,plus):
        """Reference true-box kernel, independent of bin-mean CDF implementation."""
        x=np.asarray(x,dtype=float)
        ell=L-1;centers=x+.5
        overlap,d=self._background_extra(centers,ell,plus)
        atoms=self._extended_atoms(b,plus)
        distance=centers[:,None]-atoms[None,:]
        peak=(np.clip((distance+ell/2+self.delta/2)/self.delta,0,1)
              -np.clip((distance-ell/2+self.delta/2)/self.delta,0,1)).sum(axis=1)
        return (self.class_mass[plus]+(1-self.p)*overlap+self.p/self.K*d
                +self.p/self.K*peak)/L
