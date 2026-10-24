#!/usr/bin/env python3
"""Audit (47a): independent checker for the claims C1-C4 of a kernel spec in the DEFS.md format.

usage: hl_check.py SPEC --claim BOUND [--procs P] [--maxdepth D] [--kdepth D] [--only c1,c2,c3,c4] [--cells p:q,...]
                   [--tight] [--out REPORT.json]

C3  (1/2) int K = (1/2)[c 2 a R^{2-a}/(2-a) + h^2 sum_{i,j} a_ij]  (arb upper end) <= BOUND.
C1  K >= 0 and C2  K >= 1 on D: on each knot cell S is a bicubic (exact Bernstein net on sub-rectangles, Fractions);
    Psi(u+v) (u, v >= 0) is convex in r = u + v, so c Psi >= c [Psi(rt) + Psi'(rt)(u + v - rt)] for any 0 < rt < R
    (and >= 0); the tangent is linear in (u, v) and is added to the net before taking the minimum.  In a knot
    cell cut by u + v = 1, D is split into rectangles and at most one triangle; a triangle P0 P1 P2 is handled through
    x = P0 + s[(1-t)(P1-P0) + t(P2-P0)] (S o map: exact degree (6,6) Bernstein net; u + v affine in (s, s t)).
    Quadrisection where a bound fails.  Every knot cell meeting D is checked for C2 (no skipping).
C4  A K >= 0 on every knot cell of the fundamental domain 0 <= q <= p where K is not identically 0 (where K = 0, A K >= 0
    follows from C1).  Each knot cell is cut by the graded grid GRID in both directions (finer toward the knot lines);
    on a sub-cell with local variables (xi, eta) in [-1, 1]^2:
      A_1 S = sc sum_m phi(s - m) D_m(y):  phi(t) = |t|^e Taylor-expanded to degree 3 in xi around the centre
        (Lagrange remainder |(e)_4| dist^{e-4} rho^4/24), or, for the two terms whose knot line touches the sub-cell,
        replaced by its chord (error <= (2 rho)^e tau*(1 - 1/e)/2 about the shifted chord); D_m is an exact cubic in eta.
      A_2 S = the same with u and v exchanged (a_ij = a_ji).
      c A Psi: first the monotone constant lower bound (F3 inside; F4 outside, A Psi decreasing in u and v there);
        if that fails, a Taylor model:
        inside (0 < u+v < R on the sub-cell): F(z) + G(w), z = u+v, w = u-v, to degree 3 (I^(k) in closed form via 2F1;
          |I^(k)| decreasing in d bounds the 4th-order remainders);
        outside (u+v > R, u, v > 0, each of u, v entirely below or above R): 2 J(u,v) + 2 J(v,u) to degree 2 at the
          centre (J and its first and second derivatives in closed form through M_{p,q} = int_0^{R-y} (t+y)^-p
          [(x-t)^-q + (x+t)^-q] dt, 2F1), third derivatives bounded by monotonicity in x and y.
    The sum of the polynomial parts is a polynomial of degree (3,3) in (xi, eta); its Bernstein coefficients on [-1,1]^2
    bound it below; the remainders are subtracted.  Sub-cells that fail are quadrisected up to --maxdepth.
    --tight: always evaluate the Taylor models too (for informative minima; the verdict does not depend on it).
"""
import argparse
import hashlib
import json
import os
import sys
import time
from fractions import Fraction as Fr
from math import comb, factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, ctx  # noqa: E402
from hl_core import Spec, A, pw, hyp2f1, NPOLY  # noqa: E402

ctx.prec = 128

GRID = [Fr(0), Fr(1, 256), Fr(1, 128), Fr(1, 64), Fr(1, 32), Fr(1, 16), Fr(1, 8), Fr(1, 4), Fr(3, 8), Fr(1, 2),
        Fr(5, 8), Fr(3, 4), Fr(7, 8), Fr(15, 16), Fr(31, 32), Fr(63, 64), Fr(127, 128), Fr(255, 256), Fr(1)]


# ---------------------------------------------------------------- provenance (per process)
def _sha_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def provenance():
    """sha256 of the files this process actually loaded (this script, hl_core, the python-flint modules used) and
    versions; recorded per run and per worker."""
    import flint
    import hl_core
    import flint.types.arb as _fa
    import flint.types.acb as _fc
    import flint.types.arb_mat as _fm
    import flint.types.fmpq as _fq
    d = {"hl_check.py": _sha_file(os.path.abspath(__file__)), "hl_check_path": os.path.realpath(os.path.abspath(__file__)),
         "hl_core.py": _sha_file(hl_core.__file__), "hl_core_path": os.path.realpath(hl_core.__file__),
         "flint_version": flint.__version__, "python": sys.version.split()[0]}
    for m in (_fa, _fc, _fm, _fq):
        d[os.path.basename(m.__file__)] = _sha_file(m.__file__)
    return d


def prov_digest(d):
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()


# ---------------------------------------------------------------- polynomial helpers
def power_to_bern_matrix(n, a, b):
    """matrix M (Fractions) with bern = M @ power: coefficients of sum_k p_k s^k on s in [a, b] (s = a + (b-a) t) in the
    degree-n Bernstein basis in t."""
    d = b - a
    # s^k = sum_j C(k,j) a^(k-j) d^j t^j ;  t^j = sum_{i>=j} C(i,j)/C(n,j) B_i
    M = [[Fr(0)] * (n + 1) for _ in range(n + 1)]
    for k in range(n + 1):
        for j in range(k + 1):
            cj = comb(k, j) * a ** (k - j) * d ** j
            for i in range(j, n + 1):
                M[i][k] += cj * Fr(comb(i, j), comb(n, j))
    return M


MC3 = arb_mat([[A(z) for z in row] for row in power_to_bern_matrix(3, Fr(-1), Fr(1))])
MC3T = MC3.transpose()


def poly_shift_scale(coefs, c, r):
    """coefficients in xi of sum_k coefs[k] x^k with x = c + r xi (Fractions)."""
    n = len(coefs)
    out = [Fr(0)] * n
    for k in range(n):
        if coefs[k] == 0:
            continue
        for j in range(k + 1):
            out[j] += coefs[k] * comb(k, j) * c ** (k - j) * r ** j
    return out


NEG_INF = arb('-inf')


def fin(x):
    """x for a finite arb, -inf otherwise (lower bounds: non-finite values go to the fail-closed side)."""
    return x if x.is_finite() else NEG_INF


def bern_min(P):
    """lower end of the minimum Bernstein coefficient of the 4x4 power-basis polynomial P(xi, eta) on [-1,1]^2
    (-inf if any coefficient is not finite)."""
    B = MC3 * P * MC3T
    lo = None
    for i in range(4):
        for j in range(4):
            v = B[i, j].lower()
            if not v.is_finite():
                return NEG_INF
            if lo is None or v < lo:
                lo = v
    return lo


# ---------------------------------------------------------------- C3
def check_c3(sp, claim):
    a = sp.a_
    tot = Fr(0)
    for (I, J), v in sp.acls.items():
        if I == 0 and J == 0:
            mult = 1
        elif J == 0 or I == J:
            mult = 4
        else:
            mult = 8
        tot += mult * v
    mass_psi = sp.cA * 2 * a * pw(sp.RA, 2 - a) / (2 - a)
    half = (mass_psi + A(sp.h ** 2 * tot)) / 2
    ok = half.upper() <= A(Fr(claim))
    return ok, half, tot


# ---------------------------------------------------------------- C1 / C2
def bern1d_basis(x0, x1):
    """Bx[i][k]: degree-3 Bernstein coefficients on [x0, x1] of the local basis piece N_i."""
    M = power_to_bern_matrix(3, x0, x1)
    return [[sum(M[k][d] * NPOLY[i][d] for d in range(4)) for k in range(4)] for i in range(4)]


class KChecker:
    def __init__(self, sp):
        self.sp = sp
        self.tri_cache = {}

    def tangent_options(self, rlo, rhi):
        """list of (T0, T1): c Psi(r) >= T0 + T1 r on the region r in [rlo, rhi] (arbs)."""
        sp = self.sp
        opts = [(arb(0), arb(0))]                    # Psi >= 0
        cands = set()
        for rt in (rlo, (rlo + rhi) / 2, rhi):
            if 0 < rt < sp.R:
                cands.add(rt)
        for rt in sorted(cands):
            R_t = A(rt)
            v = pw(R_t, -sp.a_) - sp.Ra
            d = -sp.a_ * pw(R_t, -sp.a_ - 1)
            opts.append((sp.cA * (v - d * R_t), sp.cA * d))
        return opts

    def rect_lower(self, p, q, x0, x1, y0, y1):
        sp = self.sp
        h = sp.h
        blk = [[sp.coef(p - 1 + i, q - 1 + j) for j in range(4)] for i in range(4)]
        Bx = bern1d_basis(x0, x1)
        By = bern1d_basis(y0, y1)
        net = [[sum(blk[i][j] * Bx[i][k] * By[j][l] for i in range(4) for j in range(4) if blk[i][j])
                for l in range(4)] for k in range(4)]
        rlo = h * (p + q + x0 + y0)
        rhi = h * (p + q + x1 + y1)
        best = None
        for (T0, T1) in self.tangent_options(rlo, rhi):
            lo = None
            for k in range(4):
                for l in range(4):
                    r = h * (p + q + x0 + (x1 - x0) * Fr(k, 3) + y0 + (y1 - y0) * Fr(l, 3))
                    val = fin((A(net[k][l]) + T0 + T1 * A(r)).lower())
                    if lo is None or val < lo:
                        lo = val
            if best is None or lo > best:
                best = lo
        return best

    def tri_lower(self, p, q, tri, s0, s1, t0, t1):
        """lower bound of K on Phi([s0,s1]x[t0,t1]), Phi(s,t) = P0 + s[(1-t)(P1-P0) + t(P2-P0)] in the local coordinates
        of knot cell (p,q) (Phi maps [0,1]^2 onto the triangle P0 P1 P2).  S o Phi is a polynomial of degree (6,6) in
        (s,t) (exact Fractions); u + v is affine in (s, s t), so the tangent bounds are bilinear in (s,t)."""
        sp = self.sp
        h = sp.h
        P0, P1, P2 = tri
        x0, y0 = P0
        a1, a2 = P1[0] - P0[0], P1[1] - P0[1]
        b1, b2 = P2[0] - P1[0], P2[1] - P1[1]          # x = x0 + s (a1 + b1 t), y = y0 + s (a2 + b2 t)
        key = (p, q, tri)
        P = self.tri_cache.get(key)
        if P is None:
            blk = [[sp.coef(p - 1 + i, q - 1 + j) for j in range(4)] for i in range(4)]
            sa = [[sum(blk[i][j] * NPOLY[i][a] * NPOLY[j][b] for i in range(4) for j in range(4) if blk[i][j])
                   for b in range(4)] for a in range(4)]
            P = [[Fr(0)] * 7 for _ in range(7)]
            for ea in range(4):
                for eb in range(4):
                    if sa[ea][eb] == 0:
                        continue
                    for i in range(ea + 1):
                        for j in range(eb + 1):
                            cf = sa[ea][eb] * comb(ea, i) * comb(eb, j) * x0 ** (ea - i) * y0 ** (eb - j)
                            if cf == 0:
                                continue
                            # s^(i+j) (a1 + b1 t)^i (a2 + b2 t)^j
                            pt = [Fr(1)]
                            for (c0, c1, k) in ((a1, b1, i), (a2, b2, j)):
                                for _ in range(k):
                                    nxt = [Fr(0)] * (len(pt) + 1)
                                    for d, z in enumerate(pt):
                                        nxt[d] += z * c0
                                        nxt[d + 1] += z * c1
                                    pt = nxt
                            for d, z in enumerate(pt):
                                P[i + j][d] += cf * z
            self.tri_cache[key] = P
        Ms = power_to_bern_matrix(6, s0, s1)
        Mt = power_to_bern_matrix(6, t0, t1)
        tmp = [[sum(Ms[i][k] * P[k][l] for k in range(7)) for l in range(7)] for i in range(7)]
        net = [[sum(tmp[i][l] * Mt[j][l] for l in range(7)) for j in range(7)] for i in range(7)]

        def rval(s_, t_):
            return h * (p + q + x0 + y0 + s_ * (a1 + a2 + (b1 + b2) * t_))
        rc = [rval(s_, t_) for s_ in (s0, s1) for t_ in (t0, t1)]         # bilinear: extremes at the corners
        best = None
        for (T0, T1) in self.tangent_options(min(rc), max(rc)):
            lo = None
            for i in range(7):
                si = s0 + (s1 - s0) * Fr(i, 6)
                for j in range(7):
                    tj = t0 + (t1 - t0) * Fr(j, 6)
                    val = fin((A(net[i][j]) + T0 + T1 * A(rval(si, tj))).lower())
                    if lo is None or val < lo:
                        lo = val
            if best is None or lo > best:
                best = lo
        return best

    def check_piece(self, p, q, target, piece, maxdepth=10):
        """piece = ('rect', x0, x1, y0, y1) or ('tri', (P0, P1, P2)) in local coordinates of knot cell (p,q); returns
        (ok, min lower bound over the accepted leaves, number of leaves, failures list)."""
        if piece[0] == "rect":
            stack = [(piece[1], piece[2], piece[3], piece[4], 0)]
        else:
            stack = [(Fr(0), Fr(1), Fr(0), Fr(1), 0)]
        leaves = 0
        mn = None
        fails = []
        while stack:
            a0, a1, b0, b1, dep = stack.pop()
            if piece[0] == "rect":
                lo = self.rect_lower(p, q, a0, a1, b0, b1)
            else:
                lo = self.tri_lower(p, q, piece[1], a0, a1, b0, b1)
            if lo >= target:
                leaves += 1
                mn = lo if mn is None or lo < mn else mn
                continue
            if dep >= maxdepth:
                fails.append((p, q, piece[0], a0, a1, b0, b1, float(lo.mid())))
                continue
            am = (a0 + a1) / 2
            bm = (b0 + b1) / 2
            stack += [(a0, am, b0, bm, dep + 1), (am, a1, b0, bm, dep + 1), (a0, am, bm, b1, dep + 1),
                      (am, a1, bm, b1, dep + 1)]
        return (not fails), mn, leaves, fails


def diamond_pieces(sp, p, q):
    """D cap knot cell (p,q) (first quadrant) in local coordinates, {x + y <= L}, L = 1/h - p - q, as rectangles and
    triangles (empty if L <= 0)."""
    L = 1 / sp.h - p - q
    if L >= 2:
        return [("rect", Fr(0), Fr(1), Fr(0), Fr(1))]
    if L <= 0:
        return []
    if L <= 1:
        return [("tri", ((Fr(0), Fr(0)), (L, Fr(0)), (Fr(0), L)))]
    m = L - 1
    return [("rect", Fr(0), m, Fr(0), Fr(1)), ("rect", m, Fr(1), Fr(0), m), ("tri", ((m, m), (Fr(1), m), (m, Fr(1))))]


def check_c12(sp, maxdepth=10):
    kc = KChecker(sp)
    out = {"c1_cells": 0, "c2_cells": 0, "c1_fail": [], "c2_fail": [], "c1_min": None, "c2_min": None,
           "c1_leaves": 0, "c2_leaves": 0}
    full = ("rect", Fr(0), Fr(1), Fr(0), Fr(1))
    P = max(sp.maxidx + 3, int(sp.R / sp.h) + 2)
    for p in range(0, P + 1):                  # C1: every knot cell of 0 <= q <= p where K is not identically 0
        for q in range(0, p + 1):
            if sp.S_zero_cell(p, q) and (p + q) * sp.h >= sp.R:
                continue
            ok, mn, lv, f = kc.check_piece(p, q, 0, full, maxdepth=maxdepth)
            out["c1_cells"] += 1
            out["c1_leaves"] += lv
            out["c1_fail"] += f
            if mn is not None and (out["c1_min"] is None or mn < out["c1_min"]):
                out["c1_min"] = mn
    P2 = int(1 / sp.h) + 1
    for p in range(0, P2 + 1):                 # C2: every knot cell of 0 <= q <= p meeting D (no skipping)
        for q in range(0, p + 1):
            pieces = diamond_pieces(sp, p, q)
            if not pieces:
                continue
            out["c2_cells"] += 1
            for pc in pieces:
                ok, mn, lv, f = kc.check_piece(p, q, 1, pc, maxdepth=maxdepth)
                out["c2_leaves"] += lv
                out["c2_fail"] += f
                if mn is not None and (out["c2_min"] is None or mn < out["c2_min"]):
                    out["c2_min"] = mn
    return out


# ---------------------------------------------------------------- C4
class AKChecker:
    def __init__(self, sp):
        self.sp = sp
        e = sp.e
        self.e = e
        self.ef = [arb(1), e, e * (e - 1), e * (e - 1) * (e - 2), e * (e - 1) * (e - 2) * (e - 3)]   # falling factorials
        self.e4 = abs(self.ef[4])
        tau = pw(1 / e, 1 / (e - 1))
        self.gstar = tau * (1 - 1 / e)                # max_{[0,1]} (t - t^e)
        self.mlist = sp.mlist
        self.nm = len(self.mlist)
        self.cache_phi = {}
        self.cache_D = {}
        self.cache_I = {}
        self.tight = False             # True: always try the Taylor models too (tighter reported minima, slower)

    # ---- 1D tables
    def phi_table(self, idx, x0, x1):
        """Taylor/chord models of phi(idx + x - m) = |idx + x - m|^e on x in [x0, x1] for every m in mlist:
        (n x 4 arb_mat of xi-coefficients, n x 1 arb_mat of remainder bounds)."""
        key = (idx, x0, x1)
        if key in self.cache_phi:
            return self.cache_phi[key]
        e = self.e
        rho = (x1 - x0) / 2
        xc = (x0 + x1) / 2
        rA = A(rho)
        C = arb_mat(self.nm, 4)
        Rm = arb_mat(self.nm, 1)
        for row, m in enumerate(self.mlist):
            tc = idx + xc - m
            T = abs(tc)
            if T > rho:
                sg = 1 if tc > 0 else -1
                TA = A(T)
                P0 = pw(TA, e)
                x = rA / TA
                C[row, 0] = P0
                C[row, 1] = sg * self.ef[1] * P0 * x
                C[row, 2] = self.ef[2] / 2 * P0 * x * x
                C[row, 3] = sg * self.ef[3] / 6 * P0 * x * x * x
                Rm[row, 0] = arb((self.e4 / 24 * pw(A(T - rho), e - 4) * rA ** 4).upper())
            elif T == rho:
                W = pw(2 * rA, e)
                g = W * self.gstar
                C[row, 0] = W / 2 - g / 2
                C[row, 1] = W / 2 if tc > 0 else -W / 2
                Rm[row, 0] = arb((g / 2).upper())
            else:
                raise ValueError("sub-interval straddles a knot")
        self.cache_phi[key] = (C, Rm)
        return C, Rm

    def D_table(self, idx, y0, y1):
        """eta-coefficients of D_m(y) on y in [y0, y1] of knot interval idx, and bounds sum |coef| >= max |D_m|."""
        key = (idx, y0, y1)
        if key in self.cache_D:
            return self.cache_D[key]
        sp = self.sp
        rho = (y1 - y0) / 2
        yc = (y0 + y1) / 2
        T = [poly_shift_scale(list(NPOLY[jj]), yc, rho) for jj in range(4)]
        C = arb_mat(self.nm, 4)
        Dm = arb_mat(self.nm, 1)
        for row, m in enumerate(self.mlist):
            dc = sp.Dcoef(m, idx)
            if not any(dc):
                continue
            co = [sum(dc[jj] * T[jj][d] for jj in range(4)) for d in range(4)]
            for d in range(4):
                C[row, d] = A(co[d])
            Dm[row, 0] = A(sum(abs(z) for z in co))
        self.cache_D[key] = (C, Dm)
        return C, Dm

    def AS_model(self, p, q, x0, x1, y0, y1):
        """(4x4 arb_mat P[k][l] of xi^k eta^l, remainder bound) for A S on the sub-cell."""
        sp = self.sp
        Px, Rx = self.phi_table(p, x0, x1)
        Dy, My = self.D_table(q, y0, y1)
        Py, Ry = self.phi_table(q, y0, y1)
        Dx, Mx = self.D_table(p, x0, x1)
        P = (Px.transpose() * Dy + Dx.transpose() * Py) * sp.sc
        rem = ((Rx.transpose() * My)[0, 0] + (Mx.transpose() * Ry)[0, 0]) * sp.sc
        return P, rem

    # ---- A Psi pieces
    def Ik(self, k, d):
        """I^(k)(d) for exact rational d > 0 (arb)."""
        key = (k, d)
        if key in self.cache_I:
            return self.cache_I[key]
        sp = self.sp
        a = sp.a_
        D = A(d)
        z = -(sp.RA - D) / D
        poch = arb(1)
        for t in range(1, k + 1):
            poch *= (t + a)
        val = (-1) ** k * poch * (sp.Ra * pw(D, -a - k) / (a + k)
                                  - pw(D, -2 * a - k) * hyp2f1(a, 2 * a + k, 2 * a + k + 1, z, bc=True) / (2 * a + k))
        self.cache_I[key] = val
        return val

    def APsi_taylor(self, p, q, x0, x1, y0, y1):
        """(4x4 polynomial of A Psi in (xi, eta), remainder) for a sub-cell inside 0 < u+v < R."""
        sp = self.sp
        h = sp.h
        a = sp.a_
        Rf = sp.R
        xc = (x0 + x1) / 2; yc = (y0 + y1) / 2
        rx = (x1 - x0) / 2; ry = (y1 - y0) / 2
        zc = h * (p + q + xc + yc); z0 = h * (p + q + x0 + y0); z1 = h * (p + q + x1 + y1)
        wc = h * (p - q + xc - yc); w0 = h * (p - q + x0 - y1); w1 = h * (p - q + x1 - y0)
        if not (0 < z0 and z1 < Rf and Rf + w0 > 0 and Rf - w1 > 0):
            raise ValueError("APsi_taylor called outside 0 < u+v < R")
        az = A(h * rx); bz = A(h * ry)
        # F(z) = C z^{-2a} + 4 I(R - z);  G(w) = 2 I(R + w) + 2 I(R - w)
        Fk = []
        Gk = []
        ff = arb(1)
        ZC = A(zc)
        for k in range(4):
            Fk.append((sp.Ca * ff * pw(ZC, -2 * a - k) + 4 * (-1) ** k * self.Ik(k, Rf - zc)) / factorial(k))
            Gk.append((2 * self.Ik(k, Rf + wc) + 2 * (-1) ** k * self.Ik(k, Rf - wc)) / factorial(k))
            ff *= (-2 * a - k)
        F4 = abs(sp.Ca) * abs(ff) * pw(A(z0), -2 * a - 4) + 4 * abs(self.Ik(4, Rf - z1))
        G4 = 2 * abs(self.Ik(4, Rf + w0)) + 2 * abs(self.Ik(4, Rf - w1))
        rem = (F4 + G4) * (az + bz) ** 4 / 24
        P = arb_mat(4, 4)
        for k in range(4):
            for i in range(k + 1):
                j = k - i
                cb = comb(k, i)
                P[i, j] += cb * az ** i * bz ** j * (Fk[k] + Gk[k] * (-1) ** j)
        return P, arb(rem.upper())

    # ---- A Psi outside (u + v > R):  A Psi = 2 J(u,v) + 2 J(v,u),  J(x,y) = int_0^L Psi(t+y) k_{1+a}(x,t) dt, L = R - y,
    #      k_q(x,t) = (x-t)^-q + (x+t)^-q;  M_{p,q}(x,y) = int_0^L (t+y)^-p k_q(x,t) dt  (>= 0, decreasing in x and y)
    def Mpq(self, X, Y, p, q):
        """M_{p,q}(x,y) for arbs x, y with 0 < y < R < x + y."""
        sp = self.sp
        Rr = sp.RA
        rr = X + Y
        if not (rr > Rr and Y > 0 and Y < Rr):
            raise ValueError("Mpq outside its domain")
        ap = 1 - p

        def F(t):          # antiderivative of t^-p (1-t)^-q
            return pw(t, ap) / ap * hyp2f1(ap, q, ap + 1, t, ac=True)
        i_ = pw(rr, 1 - p - q) * (F(Rr / rr) - F(Y / rr))           # int_y^R s^-p (x+y-s)^-q ds
        cc = X - Y
        se = p + q - 1

        def G(w):          # antiderivative of w^{p+q-2} (1 + cc w)^-q
            return pw(w, se) / se * hyp2f1(q, se, se + 1, -cc * w, bc=True)
        ii = G(1 / Y) - G(1 / Rr)                                   # int_y^R s^-p (s+x-y)^-q ds
        return i_ + ii

    def kq(self, X, Lr, q):
        return pw(X - Lr, -q) + pw(X + Lr, -q)

    def Jq(self, X, Y, q):
        """int_0^L Psi(t+y) k_q(x,t) dt = M_{a,q} - R^-a ((x-L)^{1-q} - (x+L)^{1-q})/(q-1)."""
        sp = self.sp
        Lr = sp.RA - Y
        return self.Mpq(X, Y, sp.a_, q) - sp.Ra * (pw(X - Lr, 1 - q) - pw(X + Lr, 1 - q)) / (q - 1)

    def J_derivs(self, X, Y):
        """J, J_x, J_y, J_xx, J_xy, J_yy at (x, y) (0 < y < R < x + y)."""
        sp = self.sp
        a = sp.a_
        Lr = sp.RA - Y
        J0 = self.Jq(X, Y, 1 + a)
        Jx = -(1 + a) * self.Jq(X, Y, 2 + a)
        Jxx = (1 + a) * (2 + a) * self.Jq(X, Y, 3 + a)
        Jy = -a * self.Mpq(X, Y, a + 1, 1 + a)
        Jxy = a * (1 + a) * self.Mpq(X, Y, a + 1, 2 + a)
        Jyy = a * pw(sp.RA, -a - 1) * self.kq(X, Lr, 1 + a) + a * (a + 1) * self.Mpq(X, Y, a + 2, 1 + a)
        return J0, Jx, Jy, Jxx, Jxy, Jyy

    def J3_bounds(self, x0, y0, y1):
        """upper bounds of |J_xxx|, |J_xxy|, |J_xyy|, |J_yyy| over x >= x0, y in [y0, y1] (x0 + y0 > R, y1 < R)."""
        sp = self.sp
        a = sp.a_
        Rr = sp.RA
        X0 = A(x0); Y0 = A(y0); Y1 = A(y1)
        dm = X0 + Y0 - Rr            # min of x - L
        dp = X0 + Rr - Y1            # min of x + L
        R1 = pw(Rr, -a - 1)
        R2 = pw(Rr, -a - 2)
        bxxx = (1 + a) * (2 + a) * (3 + a) * self.Jq(X0, Y0, 4 + a)
        bxxy = a * (1 + a) * (2 + a) * self.Mpq(X0, Y0, a + 1, 3 + a)
        bxyy = a * (1 + a) * R1 * (pw(dm, -2 - a) + pw(dp, -2 - a)) + a * (1 + a) ** 2 * self.Mpq(X0, Y0, a + 2, 2 + a)
        byyy = a * (1 + a) * R1 * (pw(dp, -2 - a) + pw(dm, -2 - a)) + a * (a + 1) * R2 * (pw(dm, -1 - a) + pw(dp, -1 - a)) \
            + a * (a + 1) * (a + 2) * self.Mpq(X0, Y0, a + 3, 1 + a)
        return [arb(abs(z).upper()) for z in (bxxx, bxxy, bxyy, byyy)]

    def APsi_out_taylor(self, p, q, x0, x1, y0, y1):
        """quadratic Taylor model (4x4 in xi, eta) + cubic remainder of A Psi on a sub-cell with u0 + v0 > R,
        u0, v0 > 0, and each of u, v either < R on the whole sub-cell or >= R on the whole sub-cell."""
        sp = self.sp
        h = sp.h
        Rf = sp.R
        u0 = h * (p + x0); u1 = h * (p + x1); v0 = h * (q + y0); v1 = h * (q + y1)
        uc = (u0 + u1) / 2; vc = (v0 + v1) / 2
        ru = A((u1 - u0) / 2); rv = A((v1 - v0) / 2)
        f = [arb(0)] * 6          # value, f_u, f_v, f_uu, f_uv, f_vv
        M = [arb(0)] * 4          # bounds f_uuu, f_uuv, f_uvv, f_vvv
        U = A(uc); V = A(vc)
        if v1 < Rf:               # 2 J(u, v): x = u, y = v
            J0, Jx, Jy, Jxx, Jxy, Jyy = self.J_derivs(U, V)
            f = [f[0] + 2 * J0, f[1] + 2 * Jx, f[2] + 2 * Jy, f[3] + 2 * Jxx, f[4] + 2 * Jxy, f[5] + 2 * Jyy]
            bxxx, bxxy, bxyy, byyy = self.J3_bounds(u0, v0, v1)
            M = [M[0] + 2 * bxxx, M[1] + 2 * bxxy, M[2] + 2 * bxyy, M[3] + 2 * byyy]
        if u1 < Rf:               # 2 J(v, u): x = v, y = u
            J0, Jx, Jy, Jxx, Jxy, Jyy = self.J_derivs(V, U)
            f = [f[0] + 2 * J0, f[1] + 2 * Jy, f[2] + 2 * Jx, f[3] + 2 * Jyy, f[4] + 2 * Jxy, f[5] + 2 * Jxx]
            bxxx, bxxy, bxyy, byyy = self.J3_bounds(v0, u0, u1)
            M = [M[0] + 2 * byyy, M[1] + 2 * bxyy, M[2] + 2 * bxxy, M[3] + 2 * bxxx]
        P = arb_mat(4, 4)
        P[0, 0] = f[0]
        P[1, 0] = f[1] * ru
        P[0, 1] = f[2] * rv
        P[2, 0] = f[3] * ru * ru / 2
        P[1, 1] = f[4] * ru * rv
        P[0, 2] = f[5] * rv * rv / 2
        rem = (M[0] * ru ** 3 + 3 * M[1] * ru * ru * rv + 3 * M[2] * ru * rv * rv + M[3] * rv ** 3) / 6
        return P, arb(rem.upper())

    def subcell_lower(self, p, q, x0, x1, y0, y1):
        sp = self.sp
        h = sp.h
        PS, remS = self.AS_model(p, q, x0, x1, y0, y1)
        u0 = h * (p + x0); u1 = h * (p + x1); v0 = h * (q + y0); v1 = h * (q + y1)
        const = bern_min(PS) - remS + sp.cA * sp.APsi_lower_const(u0, v0, u1, v1)
        best = fin(const.lower())
        mode = "const"
        if (best < 0 or self.tight) and u0 + v0 > 0 and u1 + v1 < sp.R:
            PP, remP = self.APsi_taylor(p, q, x0, x1, y0, y1)
            tay = fin((bern_min(PS + PP * sp.cA) - remS - sp.cA * remP).lower())
            if tay > best:
                best = tay
                mode = "taylor"
        elif (best < 0 or self.tight) and u0 + v0 > sp.R and u0 > 0 and v0 > 0 and (u1 < sp.R or u0 >= sp.R) and \
                (v1 < sp.R or v0 >= sp.R):
            PP, remP = self.APsi_out_taylor(p, q, x0, x1, y0, y1)
            tay = fin((bern_min(PS + PP * sp.cA) - remS - sp.cA * remP).lower())
            if tay > best:
                best = tay
                mode = "taylor"
        return best, mode

    def check_knot_cell(self, p, q, maxdepth):
        t0 = time.time()
        fails = []
        mn = None
        leaves = 0
        ntay = 0
        stack = []
        for i in range(len(GRID) - 1):
            for j in range(len(GRID) - 1):
                stack.append((GRID[i], GRID[i + 1], GRID[j], GRID[j + 1], 0))
        maxd = 0
        while stack:
            x0, x1, y0, y1, dep = stack.pop()
            lo, mode = self.subcell_lower(p, q, x0, x1, y0, y1)
            if lo > 0:
                leaves += 1
                ntay += (mode == "taylor")
                maxd = max(maxd, dep)
                if mn is None or lo < mn:
                    mn = lo
                continue
            if dep >= maxdepth:
                fails.append((str(x0), str(x1), str(y0), str(y1), float(lo.mid())))
                continue
            xm = (x0 + x1) / 2; ym = (y0 + y1) / 2
            stack += [(x0, xm, y0, ym, dep + 1), (xm, x1, y0, ym, dep + 1), (x0, xm, ym, y1, dep + 1),
                      (xm, x1, ym, y1, dep + 1)]
        return {"p": p, "q": q, "ok": not fails, "min": (float(mn.mid()) if mn is not None else None),
                "min_str": (mn.str(10) if mn is not None else None), "leaves": leaves, "taylor": ntay,
                "maxdepth": maxd, "fails": fails[:20], "nfail": len(fails), "sec": time.time() - t0}


_W = {}


def _init(spec_path, tight=False):
    sp = Spec(spec_path)
    _W["ak"] = AKChecker(sp)
    _W["ak"].tight = tight
    _W["prov"] = prov_digest(provenance())
    _W["spec_sha"] = _sha_file(spec_path)


def _work(args):
    p, q, maxdepth = args
    r = _W["ak"].check_knot_cell(p, q, maxdepth)
    r["prov"] = _W["prov"]
    r["spec_sha"] = _W["spec_sha"]
    return r


def c4_cells(sp):
    cells = []
    P = max(sp.maxidx + 3, int(sp.R / sp.h) + 2)
    for p in range(0, P + 1):
        for q in range(0, p + 1):
            if sp.S_zero_cell(p, q) and (p + q) * sp.h >= sp.R:
                continue
            cells.append((p, q))
    return cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("--claim", required=True)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--maxdepth", type=int, default=8)
    ap.add_argument("--kdepth", type=int, default=10)
    ap.add_argument("--only", default="c1,c2,c3,c4")
    ap.add_argument("--cells", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--tight", action="store_true")
    args = ap.parse_args()
    sp = Spec(args.spec)
    only = set(args.only.split(","))
    prov = provenance()
    report = {"spec": args.spec, "spec_sha256": _sha_file(args.spec), "claim": args.claim, "provenance": prov,
              "args": vars(args)}
    print("spec sha256 %s ; hl_check %s ; hl_core %s ; flint %s" % (report["spec_sha256"][:16], prov["hl_check.py"][:16],
          prov["hl_core.py"][:16], prov["flint_version"]), flush=True)
    t0 = time.time()
    if "c3" in only:
        ok, half, tot = check_c3(sp, args.claim)
        print("C3 mass: (1/2) int K in %s ; claim %s ; %s" % (half.str(15), args.claim, "PASS" if ok else "FAIL"),
              flush=True)
        report["c3"] = {"ok": ok, "half_mass": half.str(20), "sum_a": str(tot)}
    if "c1" in only or "c2" in only:
        o = check_c12(sp, maxdepth=args.kdepth)
        print("C1 K >= 0: %d knot cells, %d leaves, min lower %s, failures %d" %
              (o["c1_cells"], o["c1_leaves"], o["c1_min"].str(8) if o["c1_min"] is not None else None, len(o["c1_fail"])),
              flush=True)
        print("C2 K >= 1 on D: %d knot cells, %d leaves, min lower %s, failures %d" %
              (o["c2_cells"], o["c2_leaves"], o["c2_min"].str(8) if o["c2_min"] is not None else None, len(o["c2_fail"])),
              flush=True)
        report["c1"] = {"ok": not o["c1_fail"], "cells": o["c1_cells"], "leaves": o["c1_leaves"],
                        "min": o["c1_min"].str(10) if o["c1_min"] is not None else None,
                        "fails": [[str(z) for z in f] for f in o["c1_fail"][:50]]}
        report["c2"] = {"ok": not o["c2_fail"], "cells": o["c2_cells"], "leaves": o["c2_leaves"],
                        "min": o["c2_min"].str(10) if o["c2_min"] is not None else None,
                        "fails": [[str(z) for z in f] for f in o["c2_fail"][:50]]}
        print("  (%.1f s)" % (time.time() - t0), flush=True)
    if "c4" in only:
        if args.cells:
            cells = [tuple(int(z) for z in c.split(":")) for c in args.cells.split(",")]
        else:
            cells = c4_cells(sp)
        tasks = [(p, q, args.maxdepth) for (p, q) in cells]
        res = []
        if args.procs > 1:
            import multiprocessing as mp
            with mp.Pool(args.procs, initializer=_init, initargs=(args.spec, args.tight)) as pool:
                for r in pool.imap_unordered(_work, tasks):
                    res.append(r)
                    print("cell (%d,%d): %s min %s leaves %d taylor %d depth %d fails %d (%.1f s)" %
                          (r["p"], r["q"], "ok" if r["ok"] else "FAIL", r["min_str"], r["leaves"], r["taylor"],
                           r["maxdepth"], r["nfail"], r["sec"]), flush=True)
        else:
            _init(args.spec, args.tight)
            for t in tasks:
                r = _work(t)
                res.append(r)
                print("cell (%d,%d): %s min %s leaves %d taylor %d depth %d fails %d (%.1f s)" %
                      (r["p"], r["q"], "ok" if r["ok"] else "FAIL", r["min_str"], r["leaves"], r["taylor"],
                       r["maxdepth"], r["nfail"], r["sec"]), flush=True)
        dg = prov_digest(prov)
        prov_ok = all(r["prov"] == dg and r["spec_sha"] == report["spec_sha256"] for r in res)
        if not prov_ok:
            print("PROVENANCE MISMATCH between the main process and a worker", flush=True)
        nf = sum(r["nfail"] for r in res)
        mins = [r["min"] for r in res if r["min"] is not None]
        print("C4 A K >= 0: %d knot cells, %d leaves, min lower %.6g, failing sub-cells %d" %
              (len(res), sum(r["leaves"] for r in res), min(mins) if mins else float("nan"), nf), flush=True)
        report["c4"] = {"ok": nf == 0 and len(res) == len(tasks) and prov_ok, "prov_consistent": prov_ok, "cells": len(res),
                        "leaves": sum(r["leaves"] for r in res), "min": min(mins) if mins else None,
                        "results": res}
    report["sec"] = time.time() - t0
    if args.out:
        with open(args.out, "w") as f:
            json.dump(report, f, indent=1)
    allok = all(report[k]["ok"] for k in ("c1", "c2", "c3", "c4") if k in report)
    print("OVERALL (%s): %s" % (",".join(sorted(only)), "PASS" if allok else "FAIL"), flush=True)
    if not allok:
        sys.exit(1)


if __name__ == "__main__":
    main()
