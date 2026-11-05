#!/usr/bin/env python3
"""Audit (47a) core: independent evaluation of the certificate kernel K and of the jump operator A applied to K, with
rigorous ball enclosures (python-flint arb).  Written from the certificate definitions and their proofs only (not from the
first program's code).

K(u,v) = c Psi(|u|+|v|) + S(u,v),  Psi(r) = (r^-a - R^-a)_+,  S = sum_{i,j} a_ij B(u/h - i) B(v/h - j)   (DEFS.md).
A f(x) = sum_j int [f(x+t e_j) + f(x-t e_j) - 2 f(x)] |t|^{-1-a} dt.

A S (via F1, verified separately by rigorous quadrature):
  A S(u,v) = sc sum_m |u/h - m|^e D_m(v) + sc sum_m |v/h - m|^e D_m(u),   e = 3 - a,  sc = h^-a kappa,
  kappa = -2/(a(1-a)(2-a)(3-a)),  D_m(t) = sum_j delta_{m,j} B(t/h - j),  delta_{m,j} = sum_{k=-2}^{2} w_k a_{m-k,j},
  w = (1,-4,6,-4,1)   (a_ij = a_ji gives the same D_m in both directions).
A Psi (F2/F3, verified by quadrature): for u, v >= 0, r = u + v < R:
  A Psi = C_a r^-2a + 4 I(R-r) + 2 I(R+u-v) + 2 I(R+v-u),  I(d) = int_d^inf [R^-a - (s+R-d)^-a] s^{-1-a} ds;
  for r > R:  A Psi = 2 J(u,v) + 2 J(v,u),  J(x,y) = int_0^{R-y} ((t+y)^-a - R^-a) [(x-t)^{-1-a} + (x+t)^{-1-a}] dt (0 if y >= R),
  evaluated in closed form through Euler-integral / incomplete-beta antiderivatives (2F1), re-derived in this file.
"""
import json
from fractions import Fraction as Fr

from flint import arb, acb, ctx, fmpq

ctx.prec = 128
ctx.cap = 16

W4 = {-2: 1, -1: -4, 0: 6, 1: -4, 2: 1}


def A(x):
    x = Fr(x)
    return arb(fmpq(x.numerator, x.denominator))


def pw(x, e):
    """x^e for an arb x > 0 (ball must not contain 0)."""
    return (e * x.log()).exp()


def apow(x, e):
    """|x|^e for real arb x (handles exact zero); for balls containing 0 returns a ball [0, max^e] (e > 0 only)."""
    if not x.is_finite():
        raise ValueError("apow: non-finite argument")
    if x.is_zero():
        return arb(0)
    ax = abs(x)
    if ax > 0:
        return pw(ax, e)
    # ball contains 0: e must be > 0
    hi = ax.upper()
    return arb(0).union(pw(hi, e)) if hi > 0 else arb(0)


def hyp2f1(a, b, c, z, **flags):
    return acb(z).hypgeom_2f1(acb(a), acb(b), acb(c), **flags).real


# ---------------- local uniform cubic B-spline basis on [0,1] ----------------
# piece of B(x + p - i) for i = p-1, p, p+1, p+2 at local x in [0,1]
NPOLY = [  # coefficients of 1, x, x^2, x^3 (times 1/6)
    (Fr(1, 6), Fr(-3, 6), Fr(3, 6), Fr(-1, 6)),     # (1-x)^3/6
    (Fr(4, 6), Fr(0), Fr(-6, 6), Fr(3, 6)),         # (3x^3 - 6x^2 + 4)/6
    (Fr(1, 6), Fr(3, 6), Fr(3, 6), Fr(-3, 6)),      # (-3x^3 + 3x^2 + 3x + 1)/6
    (Fr(0), Fr(0), Fr(0), Fr(1, 6)),                # x^3/6
]
# Bernstein (Bezier) coefficients of the four basis pieces on [0,1]
NBEZ = [
    (Fr(1, 6), Fr(0), Fr(0), Fr(0)),
    (Fr(4, 6), Fr(4, 6), Fr(2, 6), Fr(1, 6)),
    (Fr(1, 6), Fr(2, 6), Fr(4, 6), Fr(4, 6)),
    (Fr(0), Fr(0), Fr(0), Fr(1, 6)),
]


def basis(x):
    """the four local basis values at x (Fraction -> Fractions, arb -> arbs)."""
    out = []
    for (c0, c1, c2, c3) in NPOLY:
        if isinstance(x, arb):
            out.append(A(c0) + x * (A(c1) + x * (A(c2) + x * A(c3))))
        else:
            out.append(c0 + x * (c1 + x * (c2 + x * c3)))
    return out


class Spec:
    def __init__(self, path):
        d = json.load(open(path))
        # fail closed on fields outside the DEFS.md format (e.g. extra spline components would change K)
        unknown = set(d) - {"alpha", "R", "h", "N", "c", "a"} - {"lp_val", "source", "bound", "note"}
        if unknown:
            raise ValueError("spec has fields outside the DEFS.md format: %s" % sorted(unknown))
        self.alpha = Fr(d["alpha"]); self.R = Fr(d["R"]); self.h = Fr(d["h"]); self.N = int(d["N"])
        self.c = Fr(d["c"])
        # range of validity of the closed forms and enclosures used here (1 < alpha < 2; c >= 0 for the tangent bounds)
        if not (1 < self.alpha < 2 and self.R > 0 and self.h > 0 and self.c >= 0):
            raise ValueError("spec parameters outside the range handled by this checker")
        self.acls = {}
        for k, v in d["a"].items():
            I, J = (int(z) for z in k.split(","))
            if not (0 <= J <= I <= self.N):
                raise ValueError("spec class key %s is not of the form I,J with 0 <= J <= I <= N" % k)
            self.acls[(I, J)] = Fr(v)
        self.claim = d.get("bound")
        # arb constants
        a = A(self.alpha)
        self.a_ = a
        self.e = 3 - a
        self.kappa = -2 / (a * (1 - a) * (2 - a) * (3 - a))
        self.hA = A(self.h)
        self.sc = pw(self.hA, -a) * self.kappa
        self.RA = A(self.R)
        self.Ra = pw(self.RA, -a)
        self.cA = A(self.c)
        self.Ca = -8 * (1 - a).gamma() * (2 * a).gamma() * (arb.pi() * a / 2).cos() ** 2 / (a + 1).gamma()
        # 4th differences delta[m][j]
        self.maxidx = max(I for (I, J), v in self.acls.items() if v != 0)
        M = self.maxidx + 2
        self.delta = {}
        for m in range(-M, M + 1):
            row = {}
            for j in range(-self.maxidx, self.maxidx + 1):
                s = sum(W4[k] * self.coef(m - k, j) for k in range(-2, 3))
                if s != 0:
                    row[j] = s
            if row:
                self.delta[m] = row
        self.mlist = sorted(self.delta)

    def coef(self, i, j):
        I, J = max(abs(i), abs(j)), min(abs(i), abs(j))
        return self.acls.get((I, J), Fr(0))

    # ---------- the spline S ----------
    def S_bezier(self, p, q):
        """4x4 exact Bezier net of S on the knot cell [p h,(p+1) h] x [q h,(q+1) h] in local coordinates."""
        blk = [[self.coef(p - 1 + i, q - 1 + j) for j in range(4)] for i in range(4)]
        out = [[Fr(0)] * 4 for _ in range(4)]
        for k in range(4):
            for l in range(4):
                s = Fr(0)
                for i in range(4):
                    if NBEZ[i][k] == 0:
                        continue
                    for j in range(4):
                        if blk[i][j] and NBEZ[j][l]:
                            s += blk[i][j] * NBEZ[i][k] * NBEZ[j][l]
                out[k][l] = s
        return out

    def S_zero_cell(self, p, q):
        return all(self.coef(p - 1 + i, q - 1 + j) == 0 for i in range(4) for j in range(4))

    def S_point(self, u, v):
        """exact S(u,v) for rational u, v."""
        h = self.h
        p = int((u / h) // 1); q = int((v / h) // 1)
        x = u / h - p; y = v / h - q
        bx = basis(x); by = basis(y)
        return sum(self.coef(p - 1 + i, q - 1 + j) * bx[i] * by[j] for i in range(4) for j in range(4))

    # ---------- D_m on a knot interval ----------
    def Dcoef(self, m, q):
        row = self.delta.get(m, {})
        return [row.get(q - 1 + jj, Fr(0)) for jj in range(4)]

    # ---------- A S at a point (rational u, v) ----------
    def AS_point(self, u, v, skip=()):
        """returns (A1S, A2S) as arbs at the rational point (u, v); 'skip' = set of ('u', m) / ('v', m) terms omitted."""
        h = self.h
        p = int((u / h) // 1); q = int((v / h) // 1)
        x = u / h - p; y = v / h - q
        by = basis(y); bx = basis(x)
        A1 = arb(0); A2 = arb(0)
        e = self.e
        for m in self.mlist:
            if ('u', m) not in skip:
                dc = self.Dcoef(m, q)
                Dv = sum(dc[k] * by[k] for k in range(4))
                if Dv != 0:
                    A1 += apow(A(u / h - m), e) * A(Dv)
            if ('v', m) not in skip:
                dc = self.Dcoef(m, p)
                Du = sum(dc[k] * bx[k] for k in range(4))
                if Du != 0:
                    A2 += apow(A(v / h - m), e) * A(Du)
        return self.sc * A1, self.sc * A2

    # ---------- Psi and A Psi ----------
    def Psi(self, r):
        """Psi(r) for an arb r > 0."""
        if r < self.RA:
            return pw(r, -self.a_) - self.Ra
        if r >= self.RA:
            return arb(0)
        return arb(0).union(pw(r, -self.a_) - self.Ra)

    def I_(self, d):
        a = self.a_; Rr = self.RA
        z = -(Rr - d) / d
        return self.Ra * pw(d, -a) / a - pw(d, -2 * a) * hyp2f1(a, 2 * a, 2 * a + 1, z, bc=True) / (2 * a)

    def APsi_inside(self, U, V):
        """closed form for u, v >= 0, u + v < R (arbs, points or balls)."""
        a = self.a_
        r = U + V
        return self.Ca * pw(r, -2 * a) + 4 * self.I_(self.RA - r) + 2 * self.I_(self.RA + U - V) + 2 * self.I_(self.RA + V - U)

    # outside pieces: J(x, y) = (i) + (ii) + (iii)
    def _Ji(self, x, y, q):
        """int_y^R s^-a (x + y - s)^{-q} ds  (x + y > R): incomplete beta via 2F1 antiderivative."""
        a = self.a_
        rr = x + y
        ap = 1 - a
        def F(t):
            return pw(t, ap) / ap * hyp2f1(ap, q, ap + 1, t, ac=True)
        return pw(rr, 1 - a - q) * (F(self.RA / rr) - F(y / rr))

    def _Jii(self, x, y, q, p=None):
        """int_y^R s^-p (s + x - y)^{-q} ds via s = 1/w:  int_{1/R}^{1/y} w^{p+q-2} (1 + c w)^{-q} dw, c = x - y."""
        a = self.a_
        if p is None:
            p = a
        c = x - y
        sexp = p + q - 1               # w^{sexp-1}
        def G(w):
            return pw(w, sexp) / sexp * hyp2f1(q, sexp, sexp + 1, -c * w, bc=True)
        return G(1 / y) - G(1 / self.RA)

    def J(self, x, y):
        """J(x,y) for x + y > R, y > 0 (arbs)."""
        if y >= self.RA:
            return arb(0)
        a = self.a_
        q = 1 + a
        L = self.RA - y
        i_ = self._Ji(x, y, q)
        ii = self._Jii(x, y, q)
        iii = -self.Ra * (pw(x - L, -a) - pw(x + L, -a)) / a
        return i_ + ii + iii

    def APsi_outside(self, U, V):
        return 2 * self.J(U, V) + 2 * self.J(V, U)

    def APsi_point(self, u, v):
        """A Psi at a rational point with u, v >= 0, u + v != R; +inf on the axes outside."""
        U, V = A(u), A(v)
        if u + v < self.R:
            return self.APsi_inside(U, V)
        if u + v > self.R:
            if v == 0 or u == 0:
                return arb('inf')
            return self.APsi_outside(U, V)
        return arb('inf')

    def APsi_lower_const(self, u0, v0, u1, v1):
        """monotone lower bound of A Psi over the sub-cell [u0,u1] x [v0,v1] (u, v >= 0), any position."""
        Rf = self.R
        bounds = []
        if u0 + v0 < Rf:                         # inside part nonempty
            a = self.a_
            rmax = min(u1 + v1, Rf)
            lb = self.Ca * pw(A(rmax), -2 * a) + 4 * self.I_(self.RA - A(u0 + v0)) + \
                2 * self.I_(self.RA + A(u1 - v0)) + 2 * self.I_(self.RA - A(u0 - v1))
            bounds.append(lb)
        if u1 + v1 > Rf:                         # outside part nonempty: decreasing in u and v
            bounds.append(self.APsi_outside(A(u1), A(v1)))
        if not all(b.is_finite() for b in bounds):
            return arb('-inf')                   # fail closed (NaN-safe)
        lo = bounds[0]
        for b in bounds[1:]:
            lo = lo.union(b)                     # lower endpoint of the union = min of the lower endpoints
        return arb(lo.lower())
