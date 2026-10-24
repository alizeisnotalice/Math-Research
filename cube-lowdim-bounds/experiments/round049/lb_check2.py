#!/usr/bin/env python3
"""Second checker for periodic point-mass lower bounds of c_2 (weak (1,1) constant of the centred maximal operator over
closed squares in R^2), written from NOTE_lower_bound.md (Lemmas 1 and 2) before reading the first checker's code.
Same mathematics as lb_check.py; this version prunes the enumeration so that 50-atom cells with a large reach run fast.

Family of sets (NOTE): S = the atoms in a box [x_a, x_b] x [y_a, y_b] whose sides pass through atom coordinates and
have widths <= reach; rho_S = (k - 1)/D, k = isqrt(floor(D^2 mu(S)/(4 lam))) (strictly below sqrt(mu(S)/(4 lam)));
R_S(rho) = intersection over p in S of Q(p, rho) = [max x - rho, min x + rho] x [max y - rho, min y + rho] (from the atoms
of S, not the box corners). U = union of all lattice translates; A = |U cap C|; bound = lam A / mu(C) (Lemma 2).

Enumeration (exact, with pruning that provably loses no rectangle meeting C):
 * every S of the family is the set of atoms in its own bounding box, whose widths are <= reach; so it suffices to run
   over boxes whose left side x_a and bottom side y_a are atom coordinates of S (S's own min x and min y);
 * Rbound := sqrt(Mmax/(4 lam)) with Mmax = mu(C) (ceil(reach/P_x) + 1)(ceil(reach/P_y) + 1) >= the mass of any box of
   widths <= reach, so rho_S < Rbound for every S of the family; a rectangle meeting C has min x(S) + rho > 0 and
   min y(S) + rho > 0, max x(S) - rho < P_x, max y(S) - rho < P_y, hence x_a = min x(S) in (-Rbound, P_x + Rbound - 0)
   ... (precisely: x_a > -Rbound and x_a <= max x(S) < P_x + Rbound), and the same for y_a;
 * R_S nonempty requires max(spread_x, spread_y) <= 2 rho_S, so mu(S) >= lam max(spread)^2; along increasing y_b the
   spread_y only grows, so the y_b loop stops once spread_y^2 lam > mass of the whole column within height reach;
 * candidates passing a float pre-test with a relative margin of 1e-9 are recomputed exactly (Fractions / integers);
   candidates failing it by more than the margin are certainly empty.
Union area: exact integers on the common grid 1/G (G = lcm of all denominators and D); 2-D difference array on the
compressed coordinates, coverage by prefix sums, area = sum over rows of row width x (covered length from run ends and
exact cumulative column widths).
Exit 0 iff the input is valid and lam A / mu(C) >= claimed (if a claim is given).
usage: lb_check2.py CERT.json D [CLAIMED]"""
import bisect
import json
import math
import sys
from fractions import Fraction as F
from math import gcd, isqrt

import numpy as np


def lcm(a, b):
    return a // gcd(a, b) * b


def fail(msg):
    print('FAIL:', msg, flush=True)
    print('RESULT: FAIL')
    sys.exit(1)


def frac(v, name):
    if isinstance(v, bool) or not isinstance(v, (str, int)):
        fail(f'{name}: expected a rational string or an integer')
    try:
        return F(v)
    except (ValueError, ZeroDivisionError):
        fail(f'{name}: not a rational: {v!r}')


def union_area(rects):
    """exact area of a union of integer rectangles (x0, x1, y0, y1), x0 < x1, y0 < y1."""
    if not rects:
        return 0
    xs = sorted({r[0] for r in rects} | {r[1] for r in rects})
    ys = sorted({r[2] for r in rects} | {r[3] for r in rects})
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    nx, ny = len(xs) - 1, len(ys) - 1
    ycum = [0] * (ny + 1)
    for j in range(ny):
        ycum[j + 1] = ycum[j] + (ys[j + 1] - ys[j])
    area = 0
    strip = max(1, min(nx, 30_000_000 // (ny + 1)))
    rs = [(xi[r[0]], xi[r[1]], yi[r[2]], yi[r[3]]) for r in rects]
    for s0 in range(0, nx, strip):
        s1 = min(nx, s0 + strip)
        diff = np.zeros((s1 - s0 + 1, ny + 1), dtype=np.int32)
        for (a, b, c, d) in rs:
            if b <= s0 or a >= s1:
                continue
            a2, b2 = max(a, s0) - s0, min(b, s1) - s0
            diff[a2, c] += 1; diff[a2, d] -= 1; diff[b2, c] -= 1; diff[b2, d] += 1
        cov = np.cumsum(np.cumsum(diff, axis=0), axis=1)[: s1 - s0, :ny] > 0
        for i in range(s1 - s0):
            row = cov[i]
            if not row.any():
                continue
            e = np.diff(np.concatenate(([0], row.view(np.int8), [0])))
            starts = np.nonzero(e == 1)[0]; ends = np.nonzero(e == -1)[0]
            L = sum(ycum[int(t)] - ycum[int(s)] for s, t in zip(starts, ends))
            area += (xs[s0 + i + 1] - xs[s0 + i]) * L
    return area


def main():
    if len(sys.argv) not in (3, 4):
        print(__doc__); sys.exit(2)
    cert = json.load(open(sys.argv[1]))
    D = int(sys.argv[2])
    claimed = frac(sys.argv[3], 'claimed') if len(sys.argv) == 4 else None
    for key in ('lam', 'Px', 'Py', 'pts', 'reach'):
        if key not in cert:
            fail(f'missing key {key}')
    lam, Px, Py, reach = (frac(cert[k], k) for k in ('lam', 'Px', 'Py', 'reach'))
    if not (D >= 1 and lam > 0 and Px > 0 and Py > 0 and reach >= 0):
        fail('need D >= 1, lam > 0, P_x > 0, P_y > 0, reach >= 0')
    atoms = []
    for i, a in enumerate(cert['pts']):
        if not (isinstance(a, list) and len(a) == 3):
            fail(f'atom {i}: expected [x, y, m]')
        x, y, m = (frac(v, f'atom {i}') for v in a)
        if not (0 <= x < Px and 0 <= y < Py):
            fail(f'atom {i} outside the cell [0,P_x) x [0,P_y)')
        if not m > 0:
            fail(f'atom {i}: mass <= 0')
        atoms.append((x, y, m))
    if not atoms:
        fail('no atoms')
    if len({(x, y) for x, y, _ in atoms}) != len(atoms):
        fail('two atoms at the same point')
    muC = sum(m for _, _, m in atoms)
    G = D
    for v in [Px, Py, reach] + [a[0] for a in atoms] + [a[1] for a in atoms]:
        G = lcm(G, v.denominator)
    gx, gy, gr = int(Px * G), int(Py * G), int(reach * G)
    Mmax = muC * (math.ceil(reach / Px) + 1) * (math.ceil(reach / Py) + 1)
    Rb = math.sqrt(float(Mmax / (4 * lam))) * (1 + 1e-9) + 1e-9          # float upper bound of sqrt(Mmax/(4 lam))
    gRb = int(math.ceil(Rb * G)) + 1
    # atoms needed: x in (-Rb, Px + Rb + reach], y likewise
    imin, imax = -math.ceil(Rb / float(Px)) - 1, math.ceil((float(Px) + Rb + float(reach)) / float(Px)) + 1
    jmin, jmax = -math.ceil(Rb / float(Py)) - 1, math.ceil((float(Py) + Rb + float(reach)) / float(Py)) + 1
    W = []
    for ci in range(imin, imax + 1):
        for cj in range(jmin, jmax + 1):
            for (x, y, m) in atoms:
                X, Y = int(x * G) + ci * gx, int(y * G) + cj * gy
                if -gRb < X <= gx + gRb + gr and -gRb < Y <= gy + gRb + gr:
                    W.append((X, Y, m))
    W.sort()
    Xs = [p[0] for p in W]
    xvals = sorted(set(Xs))
    four_lam = 4 * lam
    lamf = float(lam); Gf = float(G)
    rects = []
    ncand = nexact = 0
    for xa in xvals:
        if not (-gRb < xa < gx + gRb):
            continue
        i0 = bisect.bisect_left(Xs, xa)
        for xb in xvals[bisect.bisect_left(xvals, xa):]:
            if xb - xa > gr:
                break
            i1 = bisect.bisect_right(Xs, xb)
            col = sorted(W[i0:i1], key=lambda p: p[1])            # column atoms by y
            ys_ = [p[1] for p in col]
            n = len(col)
            masses = [p[2] for p in col]
            mf = np.array([float(m) for m in masses])
            cum = np.concatenate(([0.0], np.cumsum(mf)))
            wx2 = ((xb - xa) / Gf) ** 2
            for s in range(n):
                ya = ys_[s]
                if not (-gRb < ya < gy + gRb):
                    continue
                if s > 0 and ys_[s - 1] == ya:
                    continue                                       # same y_a: handled at its first index
                # column mass within height reach from ya
                tcap = bisect.bisect_right(ys_, ya + gr)
                mcol = cum[tcap] - cum[s]
                if mcol * (1 + 1e-9) < lamf * wx2:
                    continue
                mnx = mxx = None
                msum = 0.0
                for t in range(s, tcap):
                    p = col[t]
                    msum += mf[t]
                    mnx = p[0] if mnx is None or p[0] < mnx else mnx
                    mxx = p[0] if mxx is None or p[0] > mxx else mxx
                    if t + 1 < tcap and ys_[t + 1] == ys_[t]:
                        continue                                   # include all atoms at the same y before testing
                    wy = (p[1] - ya) / Gf
                    if lamf * wy * wy > mcol * (1 + 1e-9):
                        break
                    spread = max((mxx - mnx) / Gf, wy)
                    ncand += 1
                    if msum * (1 + 1e-9) < lamf * spread * spread:
                        continue                                   # certainly empty
                    # exact: S = col[s..t]
                    S = col[s:t + 1]
                    muS = sum(q[2] for q in S)
                    qv = muS / four_lam * D * D
                    k = isqrt(qv.numerator // qv.denominator)
                    rho = F(k - 1, D)
                    nexact += 1
                    if rho <= 0:
                        continue
                    r = int(rho * G)
                    x0, x1 = mxx - r, mnx + r
                    y0, y1 = p[1] - r, ya + r
                    if x0 >= x1 or y0 >= y1:
                        continue
                    cx0, cx1, cy0, cy1 = max(x0, 0), min(x1, gx), max(y0, 0), min(y1, gy)
                    if cx0 < cx1 and cy0 < cy1:
                        rects.append((cx0, cx1, cy0, cy1))
    rects = sorted(set(rects))
    import os
    if os.environ.get('LB_DUMP'):
        json.dump({'G': str(G), 'rects': [[str(v) for v in r] for r in rects]}, open(os.environ['LB_DUMP'], 'w'))
    A = F(union_area(rects), G * G)
    value = lam * A / muC
    print(f'atoms/cell {len(atoms)}, window atoms {len(W)}, mu(C) = {float(muC):.10f}, lam = {lam}, D = {D}, reach = {reach}')
    print(f'Rbound {Rb:.4f}; candidates {ncand}, exact evaluations {nexact}, distinct rectangles meeting C {len(rects)}')
    print(f'A = {float(A):.15f} of cell area {float(Px * Py):.15f}')
    print(f'lam A / mu(C) = {value.numerator}/{value.denominator}')
    print(f'              ~= {float(value):.13f}')
    if claimed is not None:
        print(f'claimed       = {claimed.numerator}/{claimed.denominator} ({"equal" if value == claimed else ("value > claim" if value > claimed else "value < claim")})')
        ok = value >= claimed
        print('RESULT:', 'PASS' if ok else 'FAIL')
        sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
