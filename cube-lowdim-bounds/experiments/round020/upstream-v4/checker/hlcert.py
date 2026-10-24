"""hl47a certificate verifier (search side; rigorous arb ball arithmetic via python-flint).

Kernel (diamond coordinates (u,v), r = |u|+|v|, D = {r <= 1}):
    K(u,v) = c * Psi(r) + S(u,v),  Psi(r) = (r^-alpha - R^-alpha)_+,
    S(u,v) = sum_{i,j=-N..N} a_ij B(u/h - i) B(v/h - j),  B = centred cubic B-spline, a_ij = a_|i||j| = a_ji.
Operator: A f(x) = sum_{j=1,2} int_R [f(x+t e_j) + f(x-t e_j) - 2 f(x)] |t|^{-1-alpha} dt  (Lean's jumpGen alpha).
Checks: (C1) K >= 0 on R^2\\{0}; (C2) K >= 1 on D\\{0}; (C3) (1/2) int K <= bound; (C4) A K >= 0 a.e. on R^2\\{0}.
Facts used (statements; proofs in PROGRESS.md / paper):
  F1  L_alpha B(s) = kappa sum_k w_k |s-k|^{3-alpha}, kappa = -2/(alpha(1-alpha)(2-alpha)(3-alpha)), w = (1,-4,6,-4,1),
      hence A S(u,v) = sc [ sum_m |u/h-m|^{p} D_m(v) + sum_m |v/h-m|^{p} D_m(u) ], p = 3-alpha, sc = h^-alpha kappa,
      D_m(y) = sum_j (sum_k w_k a_{m-k, j}) B(y/h - j).
  F2  A r^-alpha = C_alpha r^-2alpha off the axes, C_alpha = -8 Gamma(1-alpha) Gamma(2alpha) cos^2(pi alpha/2)/Gamma(alpha+1).
  F3  for r < R (u, v >= 0): A Psi = C_alpha r^-2alpha + 4 I(R-r) + 2 I(R+u-v) + 2 I(R+v-u), I decreasing, with
      I(d) = R^-alpha d^-alpha/alpha - J_alpha(d),  J_b(d) = int_d^inf (s+R-d)^-b s^{-1-alpha} ds
           = d^{-alpha-b} 2F1(b, alpha+b; alpha+b+1; -(R-d)/d)/(alpha+b),
      I' = -alpha J_{alpha+1},  I'' = alpha R^{-alpha-1} d^{-1-alpha} - alpha(alpha+1) J_{alpha+2},
      I''' = -alpha(1+alpha) R^{-alpha-1} d^{-2-alpha} + alpha(alpha+1) R^{-alpha-2} d^{-1-alpha} - alpha(alpha+1)(alpha+2) J_{alpha+3}.
  F4  for x = (u,v), u, v >= 0, r >= R: A Psi(x) >= 2 sum over the two lines of int_{-L}^{L} Psi(s, v1) ... (monotone bound by the
      far corner (u1, v1) of a cell), i.e. A Psi on a cell (all points with r >= R) >= A Psi(u1, v1).
"""
import sys, os, json, math, time
from fractions import Fraction
from flint import arb, acb, fmpq, ctx

ctx.prec = 96
W = [1, -4, 6, -4, 1]


def Q(x):
    f = Fraction(x)
    return fmpq(f.numerator, f.denominator)


def apow(x, p):
    """x^p for x > 0 (arb or exact), p arb."""
    x = arb(x)
    return (x.log() * p).exp()


def lo_of(x):
    """arb ball [lower(x), lower(x)] (exact lower endpoint)."""
    return arb(x.lower())


def lo_min(xs):
    """NaN-safe lower bound of min over the balls xs: the least exact lower endpoint; any non-finite ball gives
    [+/- inf] (no information, so every later '>= 0' test fails). Python's min() is not used on arbs because it can
    silently drop a NaN depending on the order."""
    out = None
    for x in xs:
        x = arb(x)
        if not x.is_finite():
            return arb('-inf')
        v = arb(x.lower())
        if out is None or v < out:
            out = v
    if out is None:
        raise ValueError('lo_min of an empty sequence')
    return out


def lo_max(xs):
    """lower bound of max over the balls xs from the finite ones (each finite lower endpoint is a valid lower bound of
    the max, so non-finite balls may be skipped); [+/- inf] if none is finite."""
    out = None
    for x in xs:
        x = arb(x)
        if not x.is_finite():
            continue
        v = arb(x.lower())
        if out is None or v > out:
            out = v
    return out if out is not None else arb('-inf')


def up_max(xs):
    """NaN-safe upper bound of max over the balls xs: the largest exact upper endpoint; any non-finite ball gives
    [+/- inf]."""
    out = None
    for x in xs:
        x = arb(x)
        if not x.is_finite():
            return arb('inf')
        v = arb(x.upper())
        if out is None or v > out:
            out = v
    if out is None:
        raise ValueError('up_max of an empty sequence')
    return out


def min_lo(a, b):
    """exact lower endpoint of min(a, b) as a point arb (sufficient for lower bounds); NaN-safe."""
    return lo_min((a, b))


def max_lo(a, b):
    """a point arb whose value is <= max over the balls (valid lower bound for max)."""
    return lo_max((a, b))


_BP = {-2: [8, 12, 6, 1], -1: [4, 0, -6, -3], 0: [4, 0, -6, 3], 1: [8, -12, 6, -1]}


def bern_of_piece(off, t0, t1):
    """B(x - off), local x in [t0,t1] subset [0,1], off in {-1,0,1,2}: exact Bernstein coefficients (fmpq), degree 3."""
    m = -off
    p = [fmpq(c, 6) for c in _BP[m]]
    cx = [fmpq(0)] * 4
    for k, pk in enumerate(p):
        for j in range(k + 1):
            cx[j] += pk * math.comb(k, j) * fmpq(-off) ** (k - j)
    cs = [fmpq(0)] * 4
    dd = t1 - t0
    for k, ck in enumerate(cx):
        for j in range(k + 1):
            cs[j] += ck * math.comb(k, j) * t0 ** (k - j) * dd ** j
    return [sum((fmpq(math.comb(i, j), math.comb(3, j)) * cs[j] for j in range(i + 1)), fmpq(0)) for i in range(4)]


def elev1(vec, newdeg):
    n = len(vec) - 1
    while n < newdeg:
        nv = [None] * (n + 2)
        nv[0] = vec[0]; nv[n + 1] = vec[n]
        for i in range(1, n + 1):
            a = arb(fmpq(i, n + 1))
            nv[i] = a * vec[i - 1] + (1 - a) * vec[i]
        vec = nv; n += 1
    return vec


def _elev2(B, d0, d1):
    """elevate a Bernstein array B[i][j] (deg n0 x n1) to (d0, d1)."""
    rows = [elev1(list(row), d1) for row in B]
    cols = [elev1([rows[i][j] for i in range(len(rows))], d0) for j in range(d1 + 1)]
    return [[cols[j][i] for j in range(d1 + 1)] for i in range(d0 + 1)]


class SplineComp:
    """one D4-symmetric tensor cubic B-spline component S(u,v) = sum a_ij B(u/h - i) B(v/h - j), |i|,|j| <= N."""

    def __init__(self, h, N, a, al):
        self.h_q = Q(h); self.h = arb(self.h_q); self.N = N = int(N)
        A = [[fmpq(0)] * (2 * N + 1) for _ in range(2 * N + 1)]
        for key, val in a.items():
            I_, J_ = map(int, key.split(','))
            v = Q(val)
            for i in {I_, -I_}:
                for j in {J_, -J_}:
                    A[i + N][j + N] = v; A[j + N][i + N] = v
        self.A = A
        self.D = {}
        for m in range(-N - 2, N + 3):
            row = []
            for j in range(-N, N + 1):
                s_ = fmpq(0)
                for kk in range(5):
                    i = m - (kk - 2)
                    if -N <= i <= N:
                        s_ += W[kk] * A[i + N][j + N]
                row.append(s_)
            if any(x != 0 for x in row):
                self.D[m] = row
        self.p = 3 - al
        self.kap = -2 / (al * (1 - al) * (2 - al) * (3 - al))
        self.sc = apow(self.h, -al) * self.kap
        self._col = {}; self._row = {}

class Kern:
    def __init__(self, spec):
        self.alpha_q = Q(spec['alpha']); self.R_q = Q(spec['R']); self.c_q = Q(spec['c'])
        self.al = arb(self.alpha_q); self.R = arb(self.R_q); self.c = arb(self.c_q)
        al = self.al
        self.comps = [SplineComp(spec['h'], spec['N'], spec['a'], al)]
        for ex in spec.get('extra', []):
            self.comps.append(SplineComp(ex['h'], ex['N'], ex['a'], al))
        c0 = self.comps[0]
        self.h_q, self.h, self.N, self.A, self.D = c0.h_q, c0.h, c0.N, c0.A, c0.D     # primary component (compat)
        self.p = 3 - al
        self.kap = -2 / (al * (1 - al) * (2 - al) * (3 - al))
        pi = arb.pi()
        self.Calpha = -8 * (1 - al).gamma() * (2 * al).gamma() * (pi * al / 2).cos() ** 2 / (al + 1).gamma()

    # ---------------------------------------------------------------- singular function pieces
    def J(self, b, d):
        al, R = self.al, self.R
        d = arb(d)
        z = -(R - d) / d
        return apow(d, -al - b) * z.hypgeom_2f1(b, al + b, al + b + 1) / (al + b)

    def Jbound(self, b, d0):
        """sup_{d >= d0} J_b(d): for d <= R, J_b(d) <= d^{-alpha-b}/(alpha+b) (since s+R-d >= s); for d >= R,
        J_b(d) <= R^-b d^-alpha/alpha (since s+R-d >= R for s >= d). Both bounds decrease in d."""
        al, R = self.al, self.R
        d0 = arb(d0)
        if d0 >= R:
            return up_max((apow(R, -b) * apow(d0, -al) / al,))
        return up_max((apow(d0, -al - b) / (al + b), apow(R, -b - al) / al))

    def I0(self, d):
        al, R = self.al, self.R
        return apow(R, -al) * apow(d, -al) / al - self.J(al, d)

    def I1(self, d):
        return -self.al * self.J(self.al + 1, d)

    def I2(self, d):
        al, R = self.al, self.R
        return al * apow(R, -al - 1) * apow(d, -1 - al) - al * (al + 1) * self.J(al + 2, d)

    def I3bound(self, d0):
        al, R = self.al, self.R
        d0 = arb(d0)
        return (al * (1 + al) * apow(R, -al - 1) * apow(d0, -2 - al) + al * (al + 1) * apow(R, -al - 2) * apow(d0, -1 - al)
                + al * (al + 1) * (al + 2) * self.Jbound(al + 3, d0))

    # ---------------------------------------------------------------- column / row caches for the spline part
    def col_data(self, s0, s1, comp=None):
        """per u-interval [s0,s1] (units of h, inside one knot interval [k,k+1]):
        for every knot m with D_m != 0: Bernstein-in-tau (deg 2) lower/upper data of t(s-m) = |s-m|^p.
        Returns (k, dict m -> ('adj', tan0, tan1, ch0, ch1) or ('far', f0, b1, f1, t3rem))."""
        comp = comp or self.comps[0]
        key = (s0, s1)
        if key in comp._col:
            return comp._col[key]
        p = self.p
        k = int(math.floor(float(s0) + 1e-12))
        if not (fmpq(k) <= s0 and s1 <= fmpq(k + 1)):
            raise ValueError('column interval must lie in one knot interval')
        ds = s1 - s0; dsA = arb(ds); smid = (s0 + s1) / 2
        out = {}
        for m in comp.D:
            x0, x1, xm = s0 - m, s1 - m, smid - m
            axm = abs(xm); sg = 1 if xm > 0 else -1
            tm = apow(axm, p); t1 = p * apow(axm, p - 1) * sg
            if m == k or m == k + 1:
                ch0 = apow(abs(x0), p) if x0 != 0 else arb(0)
                ch1 = apow(abs(x1), p) if x1 != 0 else arb(0)
                out[m] = ('adj', tm - t1 * dsA / 2, tm + t1 * dsA / 2, ch0, ch1)
            else:
                t2 = p * (p - 1) * apow(axm, p - 2)
                Aa, Bb, Cc = tm, t1 * dsA, t2 * dsA * dsA / 2
                f0 = Aa - Bb / 2 + Cc / 4; f1 = Aa + Bb / 2 + Cc / 4; b1 = 2 * Aa - (f0 + f1) / 2
                xmin = min(abs(x0), abs(x1))
                t3 = p * (p - 1) * abs(p - 2) * apow(xmin, p - 3)
                out[m] = ('far', f0, b1, f1, t3 * (dsA / 2) ** 3 / 6)
        res = (k, out)
        comp._col[key] = res
        return res

    def row_data(self, q0, q1, comp=None):
        """per v-interval [q0,q1] (units of h, inside [l,l+1]): Bernstein (deg 3 in sigma) of sc*D_m(v) for all m."""
        comp = comp or self.comps[0]
        key = (q0, q1)
        if key in comp._row:
            return comp._row[key]
        N = comp.N
        l = int(math.floor(float(q0) + 1e-12))
        Bv = {l + off: bern_of_piece(off, q0 - l, q1 - l) for off in (-1, 0, 1, 2)}
        out = {}
        for m, row in comp.D.items():
            dv = [fmpq(0)] * 4
            nz = False
            for j, bc in Bv.items():
                if -N <= j <= N and row[j + N] != 0:
                    nz = True
                    for qq in range(4):
                        dv[qq] += row[j + N] * bc[qq]
            if nz and any(x != 0 for x in dv):
                out[m] = [comp.sc * arb(x) for x in dv]
        comp._row[key] = out
        return out

    def _part(self, cd, rd):
        """lower Bernstein array [tau deg 2][sigma deg 3] and remainder for sum_m t(s-m) * (sc D_m)(sigma)."""
        k, tdat = cd
        res = [[arb(0)] * 4 for _ in range(3)]
        rem = arb(0)
        for m, dco in rd.items():
            td = tdat.get(m)
            if td is None:
                continue
            if td[0] == 'adj':
                _, tan0, tan1, ch0, ch1 = td
                for qq in range(4):
                    d = dco[qq]
                    if d >= 0:
                        L0, L1 = d * tan0, d * tan1
                    elif d <= 0:
                        L0, L1 = d * ch0, d * ch1
                    else:
                        L0 = min_lo(d * tan0, d * ch0); L1 = min_lo(d * tan1, d * ch1)
                    # a linear function in tau with end values >= L0, L1: Bernstein deg 2 = (L0, (L0+L1)/2, L1)
                    res[0][qq] += L0; res[1][qq] += (L0 + L1) / 2; res[2][qq] += L1
            else:
                _, f0, b1, f1, r3 = td
                for qq in range(4):
                    d = dco[qq]
                    res[0][qq] += f0 * d; res[1][qq] += b1 * d; res[2][qq] += f1 * d
                dmax = up_max(abs(x) for x in dco)
                rem += dmax * r3
        return res, rem

    def spline_array(self, u0, u1, v0, v1):
        """Bernstein (3,3) array (tau = u-local, sigma = v-local) of a lower polynomial for A S on the cell, and remainder,
        summed over all spline components (the cell must lie inside one knot cell of every component)."""
        Z = [[arb(0)] * 4 for _ in range(4)]
        rem = arb(0)
        for comp in self.comps:
            h = comp.h_q
            s0, s1, q0, q1 = u0 / h, u1 / h, v0 / h, v1 / h
            ru, remu = self._part(self.col_data(s0, s1, comp), self.row_data(q0, q1, comp))
            rv, remv = self._part(self.col_data(q0, q1, comp), self.row_data(s0, s1, comp))
            Zu = _elev2(ru, 3, 3)
            Zv = _elev2(rv, 3, 3)
            for i in range(4):
                for j in range(4):
                    Z[i][j] = Z[i][j] + Zu[i][j] + Zv[j][i]
            rem = rem + remu + remv
        return Z, rem

    def spline_lower(self, u0, u1, v0, v1):
        Z, rem = self.spline_array(u0, u1, v0, v1)
        return lo_min(z for row in Z for z in row) - rem

    def A_lower(self, u0, u1, v0, v1):
        """rigorous lower bound of A K on the cell: spline array + singular Taylor array combined before the min."""
        Z, rem = self.spline_array(u0, u1, v0, v1)
        sing = self.sing_parts(u0, u1, v0, v1)        # ('poly', B[3][3] arb (already times c), rem) or ('const', lb)
        if sing[0] == 'poly':
            _, Bs, rs = sing
            lo = lo_min(Z[i][j] + Bs[i][j] for i in range(4) for j in range(4))
            return lo - rem - rs
        lo = lo_min(z for row in Z for z in row)
        return lo - rem + sing[1]

    # ---------------------------------------------------------------- singular part on a cell
    def sing_lower(self, u0, u1, v0, v1):
        sp_ = self.sing_parts(u0, u1, v0, v1)
        if sp_[0] == 'poly':
            return lo_min(x for row in sp_[1] for x in row) - sp_[2]
        return sp_[1]

    def sing_mono(self, u0, u1, v0, v1):
        al, R = self.al, self.R
        rmin, rmax = u0 + v0, u1 + v1
        cands = []
        if rmin < self.R_q:
            rr = min(rmax, self.R_q)
            cands.append(self.Calpha * apow(arb(rr), -2 * al) + 4 * self.I0(R - arb(rmin))
                         + 2 * self.I0(R + arb(u1 - v0)) + 2 * self.I0(R + arb(v1 - u0)))
        if rmax > self.R_q:
            cands.append(self.Apsi_outside_point(u1, v1))
        lo = cands[0]
        for x in cands[1:]:
            lo = min_lo(lo, x)
        return self.c * lo

    # ---------------------------------------------------------------- outside (r >= R) Taylor-1 model
    def _Hint(self, a, b, kind):
        """rigorous integrals for H(a,b) = 2 int_0^L psi(s,b) W(a,s) ds (L = R - b), psi = (s+b)^-al - R^-al,
        W = (a-s)^{-1-al} + (a+s)^{-1-al}.  kind: 'H', 'Ha', 'Hb', 'Haa', 'Hab', 'Hbb_int' (integral part of H_bb)."""
        al = self.al
        L = self.R_q - b
        if L <= 0:
            return arb(0)
        aa, bb, alc = acb(a), acb(b), acb(al)
        Ra = acb(apow(self.R, -al))
        pw = lambda z, e: ((z.log()) * e).exp()
        if kind == 'H':
            f = lambda t, an: (pw(t + bb, -alc) - Ra) * (pw(aa - t, -1 - alc) + pw(aa + t, -1 - alc))
        elif kind == 'Ha':
            f = lambda t, an: (pw(t + bb, -alc) - Ra) * (-(1 + alc)) * (pw(aa - t, -2 - alc) + pw(aa + t, -2 - alc))
        elif kind == 'Hb':
            f = lambda t, an: (-alc) * pw(t + bb, -alc - 1) * (pw(aa - t, -1 - alc) + pw(aa + t, -1 - alc))
        elif kind == 'Haa':
            f = lambda t, an: (pw(t + bb, -alc) - Ra) * (1 + alc) * (2 + alc) * (pw(aa - t, -3 - alc) + pw(aa + t, -3 - alc))
        elif kind == 'Hab':
            f = lambda t, an: (-alc) * pw(t + bb, -alc - 1) * (-(1 + alc)) * (pw(aa - t, -2 - alc) + pw(aa + t, -2 - alc))
        elif kind == 'Hbb_int':
            f = lambda t, an: alc * (alc + 1) * pw(t + bb, -alc - 2) * (pw(aa - t, -1 - alc) + pw(aa + t, -1 - alc))
        val = acb.integral(f, acb(0), acb(L))
        return 2 * val.real

    def _outside_taylor(self, u0, u1, v0, v1):
        """linear lower model of A Psi on a cell with u0 + v0 > R (all points outside). Returns (B[3][3] arb, rem) or (None, None)."""
        al, R = self.al, self.R
        uc, vc = (u0 + u1) / 2, (v0 + v1) / 2
        du, dv = (u1 - u0) / 2, (v1 - v0) / 2
        F0 = self._Hint(uc, vc, 'H') + self._Hint(vc, uc, 'H')
        Fu = self._Hint(uc, vc, 'Ha') + self._Hint(vc, uc, 'Hb')
        Fv = self._Hint(uc, vc, 'Hb') + self._Hint(vc, uc, 'Ha')
        # Hessian bounds over the cell (monotone integrands: worst at the near corner (u0, v0))
        def W_L(a, b):
            L = self.R_q - b
            if L <= 0:
                return arb(0)
            return apow(arb(a - L), -1 - al) + apow(arb(a + L), -1 - al)
        def Hbb_bound(a, b):
            if self.R_q - b <= 0:
                return arb(0)
            return 2 * al * apow(R, -al - 1) * W_L(a, b) + abs(self._Hint(a, b, 'Hbb_int'))
        Muu = abs(self._Hint(u0, v0, 'Haa')) + Hbb_bound(v0, u0)
        Mvv = Hbb_bound(u0, v0) + abs(self._Hint(v0, u0, 'Haa'))
        Muv = abs(self._Hint(u0, v0, 'Hab')) + abs(self._Hint(v0, u0, 'Hab'))
        dua, dva = arb(du), arb(dv)
        rem = (Muu * dua ** 2 + 2 * Muv * dua * dva + Mvv * dva ** 2) / 2
        # linear model values at the corners -> Bernstein (1,1) -> elevate (3,3)
        B1 = [[F0 + Fu * (-dua if i == 0 else dua) + Fv * (-dva if j == 0 else dva) for j in range(2)] for i in range(2)]
        if not (rem.is_finite() and all(x.is_finite() for row in B1 for x in row)):
            return None, None
        return _elev2(B1, 3, 3), rem

    def sing_parts(self, u0, u1, v0, v1):
        mono = self.sing_mono(u0, u1, v0, v1)
        if u0 + v0 > self.R_q and u0 + v0 - self.R_q < fmpq(1, 2):
            B, rem = self._outside_taylor(u0, u1, v0, v1)
            if B is not None:
                Bc = [[self.c * x for x in row] for row in B]
                return ('poly', Bc, self.c * rem)
        if u1 + v1 < self.R_q and u0 + v0 > 0:
            B, rem = self._sing_inside_taylor(u0, u1, v0, v1)
            if B is not None:
                Bc = [[self.c * x for x in row] for row in B]
                remc = self.c * rem
                tay_lo = lo_min(x for row in Bc for x in row) - remc
                if tay_lo.lower() >= mono.lower() - 1e-3 * abs(float(mono.mid())) or True:
                    # use the polynomial form (it captures gradient cancellation with the spline part); fall back to the
                    # constant if the polynomial bound is much worse than the monotone constant
                    if tay_lo.lower() > mono.lower() - 10 * abs(float(mono.mid())) - 10:
                        return ('poly', Bc, remc)
        return ('const', mono)

    def _sing_inside_taylor(self, u0, u1, v0, v1):
        al, R = self.al, self.R
        if u0 + v0 == 0:
            return None, None
        uc, vc = (u0 + u1) / 2, (v0 + v1) / 2
        du, dv = (u1 - u0) / 2, (v1 - v0) / 2
        zc = arb(uc + vc); wc = arb(uc - vc)
        dR = R - zc
        g0 = self.Calpha * apow(zc, -2 * al) + 4 * self.I0(dR)
        g1 = -2 * al * self.Calpha * apow(zc, -2 * al - 1) - 4 * self.I1(dR)
        g2 = 2 * al * (2 * al + 1) * self.Calpha * apow(zc, -2 * al - 2) + 4 * self.I2(dR)
        g3b = 2 * al * (2 * al + 1) * (2 * al + 2) * self.Calpha * apow(arb(u0 + v0), -2 * al - 3) \
            + 4 * self.I3bound(R - arb(u1 + v1))
        Gp, Gm = R + wc, R - wc
        G0 = 2 * self.I0(Gp) + 2 * self.I0(Gm)
        G1 = 2 * self.I1(Gp) - 2 * self.I1(Gm)
        G2 = 2 * self.I2(Gp) + 2 * self.I2(Gm)
        wmax = max(abs(u1 - v0), abs(u0 - v1))
        G3b = 4 * self.I3bound(R - arb(wmax))
        dua, dva = arb(du), arb(dv)

        def F(x, y):
            return g0 + g1 * (x + y) + g2 * (x + y) ** 2 / 2 + G0 + G1 * (x - y) + G2 * (x - y) ** 2 / 2
        xs = [-dua, arb(0), dua]; ys = [-dva, arb(0), dva]
        V = [[F(x, y) for y in ys] for x in xs]
        tb = lambda f0, fm, f1: [f0, 2 * fm - (f0 + f1) / 2, f1]
        rowsB = [tb(*V[i]) for i in range(3)]                      # rowsB[i][j]: x-node i, sigma-Bernstein j
        B2 = [[None] * 3 for _ in range(3)]
        for j in range(3):
            col = tb(rowsB[0][j], rowsB[1][j], rowsB[2][j])        # tau-Bernstein for sigma-index j
            for i in range(3):
                B2[i][j] = col[i]
        B3 = _elev2(B2, 3, 3)
        rem = (g3b + G3b) * arb(du + dv) ** 3 / 6
        if not (rem.is_finite() and all(x.is_finite() for row in B3 for x in row)):
            return None, None
        return B3, rem

    def Apsi_outside_point(self, u, v):
        """A Psi at (u, v), u + v > R, as an arb (F4 lower bound expression evaluated exactly)."""
        al, R = self.al, self.R
        tot = arb(0)
        for (a, b) in ((u, v), (v, u)):
            L = self.R_q - b
            if L <= 0:
                continue
            aa, bb, alc = acb(a), acb(b), acb(al)
            Ra = acb(apow(R, -al))

            def f(s, analytic):
                return ((((s + bb).log()) * (-alc)).exp() - Ra) * \
                    ((((aa - s).log()) * (-1 - alc)).exp() + (((aa + s).log()) * (-1 - alc)).exp())
            val = acb.integral(f, acb(0), acb(L))
            tot += 2 * val.real
        return tot

    # ---------------------------------------------------------------- K
    def S_power(self, u0, u1, v0, v1):
        """exact power coefficients C[x][y] (fmpq) of S (all components) in local (tau, sigma) on a cell."""
        C = [[fmpq(0)] * 4 for _ in range(4)]
        for comp in self.comps:
            Cc = self._S_power_comp(comp, u0, u1, v0, v1)
            for x in range(4):
                for y in range(4):
                    C[x][y] += Cc[x][y]
        return C

    def _S_power_comp(self, comp, u0, u1, v0, v1):
        h, N = comp.h_q, comp.N
        s0, s1, q0, q1 = u0 / h, u1 / h, v0 / h, v1 / h
        k = int(math.floor(float(s0) + 1e-12)); l = int(math.floor(float(q0) + 1e-12))
        if not (fmpq(k) <= s0 and s1 <= fmpq(k + 1) and fmpq(l) <= q0 and q1 <= fmpq(l + 1)):
            raise ValueError('cell must lie in one knot cell')

        def pw(off, t0, t1):
            m = -off
            p = [fmpq(c, 6) for c in _BP[m]]
            cx = [fmpq(0)] * 4
            for kk, pk in enumerate(p):
                for j in range(kk + 1):
                    cx[j] += pk * math.comb(kk, j) * fmpq(-off) ** (kk - j)
            cs = [fmpq(0)] * 4
            dd = t1 - t0
            for kk, ck in enumerate(cx):
                for j in range(kk + 1):
                    cs[j] += ck * math.comb(kk, j) * t0 ** (kk - j) * dd ** j
            return cs
        Pu = {k + off: pw(off, s0 - k, s1 - k) for off in (-1, 0, 1, 2)}
        Pv = {l + off: pw(off, q0 - l, q1 - l) for off in (-1, 0, 1, 2)}
        C = [[fmpq(0)] * 4 for _ in range(4)]
        for i, pu in Pu.items():
            if not -N <= i <= N:
                continue
            for j, pv in Pv.items():
                if not -N <= j <= N:
                    continue
                aij = comp.A[i + N][j + N]
                if aij == 0:
                    continue
                for x in range(4):
                    for y in range(4):
                        C[x][y] += aij * pu[x] * pv[y]
        return C

    def T3_power(self, u0, u1, v0, v1):
        """power coefficients (arb) in (tau, sigma) of the cubic Taylor polynomial of r^-alpha - R^-alpha at the centre;
        it is <= Psi on the cell when u0+v0 > 0 (4th derivative of r^-alpha > 0; (y)_+ >= y)."""
        al, R = self.al, self.R
        uc, vc = (u0 + u1) / 2, (v0 + v1) / 2
        du, dv = (u1 - u0) / 2, (v1 - v0) / 2
        rc = arb(uc + vc)
        t = [apow(rc, -al) - apow(R, -al), -al * apow(rc, -al - 1), al * (al + 1) * apow(rc, -al - 2) / 2,
             -al * (al + 1) * (al + 2) * apow(rc, -al - 3) / 6]
        base = {(0, 0): arb(-du - dv), (1, 0): arb(2 * du), (0, 1): arb(2 * dv)}
        PS = [[arb(0)] * 4 for _ in range(4)]
        powk = {(0, 0): arb(1)}
        for n in range(4):
            for (a_, b_), v_ in powk.items():
                PS[a_][b_] += t[n] * v_
            newp = {}
            for (a_, b_), v_ in powk.items():
                for (c_, d_), w_ in base.items():
                    key = (a_ + c_, b_ + d_)
                    newp[key] = newp.get(key, arb(0)) + v_ * w_
            powk = newp
        return PS

    @staticmethod
    def power_to_tensor_bern(C):
        """3x3-degree power coefficients -> tensor Bernstein coefficients."""
        def conv(vec):
            return [sum((arb(fmpq(math.comb(i, j), math.comb(3, j))) * vec[j] for j in range(i + 1)), arb(0)) for i in range(4)]
        rows = [conv([arb(C[x][y]) for y in range(4)]) for x in range(4)]
        cols = [conv([rows[x][y] for x in range(4)]) for y in range(4)]
        return [[cols[y][x] for y in range(4)] for x in range(4)]

    @staticmethod
    def power_to_tri_bb(C):
        """power coefficients C[i][j] -> degree-6 BB coefficients on the triangle tau, sigma >= 0, tau+sigma <= 1."""
        n = 6
        from math import factorial
        bb = {}
        for i in range(4):
            for j in range(4):
                cij = C[i][j]
                if isinstance(cij, fmpq) and cij == 0:
                    continue
                m = n - i - j
                for a in range(m + 1):
                    for b in range(m - a + 1):
                        c = m - a - b
                        coef = arb(fmpq(factorial(m), factorial(a) * factorial(b) * factorial(c)))
                        key = (a, b + i, c + j)
                        bb[key] = bb.get(key, arb(0)) + arb(cij) * coef
        out = []
        for a in range(n + 1):
            for b in range(n - a + 1):
                c = n - a - b
                out.append(bb.get((a, b, c), arb(0)) / arb(fmpq(factorial(n), factorial(a) * factorial(b) * factorial(c))))
        return out

    @staticmethod
    def power_to_tensor_bern_exact(C):
        def conv(vec):
            return [sum((fmpq(math.comb(i, j), math.comb(3, j)) * vec[j] for j in range(i + 1)), fmpq(0)) for i in range(4)]
        rows = [conv([C[x][y] for y in range(4)]) for x in range(4)]
        cols = [conv([rows[x][y] for x in range(4)]) for y in range(4)]
        return [[cols[y][x] for y in range(4)] for x in range(4)]

    @staticmethod
    def power_to_tri_bb_exact(C):
        n = 6
        from math import factorial
        bb = {}
        for i in range(4):
            for j in range(4):
                cij = C[i][j]
                if cij == 0:
                    continue
                m = n - i - j
                for a in range(m + 1):
                    for b in range(m - a + 1):
                        c = m - a - b
                        key = (a, b + i, c + j)
                        bb[key] = bb.get(key, fmpq(0)) + cij * fmpq(factorial(m), factorial(a) * factorial(b) * factorial(c))
        out = []
        for a in range(n + 1):
            for b in range(n - a + 1):
                c = n - a - b
                out.append(bb.get((a, b, c), fmpq(0)) / fmpq(factorial(n), factorial(a) * factorial(b) * factorial(c)))
        return out

    def K_bounds(self, u0, u1, v0, v1, triangle=False):
        """rigorous lower bound of K on the closed cell (and on its lower-left triangle if triangle=True).
        Spline part: exact (fmpq) Bernstein/BB coefficients; singular part: Bernstein/BB of its cubic Taylor lower bound (arb)."""
        C = self.S_power(u0, u1, v0, v1)
        rmin, rmax = u0 + v0, u1 + v1
        Sx = self.power_to_tri_bb_exact(C) if triangle else [x for row in self.power_to_tensor_bern_exact(C) for x in row]
        if rmin >= self.R_q:
            return arb(min(Sx))
        if rmin == 0:
            psi = apow(arb(rmax), -self.al) - apow(self.R, -self.al) if rmax < self.R_q else arb(0)
            return arb(min(Sx)) + self.c * psi
        T = self.T3_power(u0, u1, v0, v1)
        Tx = self.power_to_tri_bb(T) if triangle else [x for row in self.power_to_tensor_bern(T) for x in row]
        cands = [lo_min(arb(a) + self.c * b for a, b in zip(Sx, Tx))]
        # tangent line of the convex function r -> r^-alpha - R^-alpha at r_t = min(rmax, R): a lower bound of Psi on the
        # cell that is exact at the far corner (useful where K touches 0 at r = R)
        al, R = self.al, self.R
        rt = min(rmax, self.R_q)
        p0 = arb(0) if rt == self.R_q else apow(arb(rt), -al) - apow(R, -al)     # exact 0 at r_t = R
        p1 = -al * apow(arb(rt), -al - 1)
        du, dv = (u1 - u0) / 2, (v1 - v0) / 2
        rc = (u0 + u1 + v0 + v1) / 2
        # a linear function's Bernstein (tensor deg 3) / triangle-BB (deg 6) coefficients are its values at the domain points;
        # with exact rational factors (r - rt) the coefficient at the far corner is exactly p0
        def lin_at(t_, s_):
            fac = (rc - du - dv) + 2 * du * t_ + 2 * dv * s_ - rt         # exact fmpq
            return p0 + (p1 * arb(fac) if fac != 0 else arb(0))
        if triangle:
            Lx = [lin_at(fmpq(b, 6), fmpq(c, 6)) for a in range(7) for b in range(7 - a) for c in [6 - a - b]]
        else:
            Lx = [lin_at(fmpq(i, 3), fmpq(j, 3)) for i in range(4) for j in range(4)]
        cands.append(lo_min(arb(a) + self.c * b for a, b in zip(Sx, Lx)))
        if rmax >= self.R_q:
            # the cell reaches r = R: Psi >= 0 is also valid
            cands.append(arb(min(Sx)))
        return lo_max(cands)

    def K_lower_cell(self, u0, u1, v0, v1):
        return self.K_bounds(u0, u1, v0, v1, triangle=False)

    def mass_upper(self):
        al = self.al
        spl = arb(0)
        for comp in self.comps:
            tot = fmpq(0)
            for row in comp.A:
                for x in row:
                    tot += x
            spl += arb(comp.h_q ** 2 * tot)
        return (self.c * 2 * al * apow(self.R, 2 - al) / (2 - al) + spl) / 2
