"""Exact lower bound for the weak-type ratio of a periodic point-mass configuration: integer-arithmetic version of
check_exact.py (same mathematics, much faster).

All coordinates are integers over a common denominator Dx, all masses integers over Dm. For each box with bounds at atom
coordinates (sides <= reach, bounds within rho_max of the cell), S = atoms in the box, mu(S) = Ms/Dm, and
rho_S = (isqrt(floor(Ms Dr^2 / (4 lam Dm))) - 1) / Dr, strictly below sqrt(mu(S)/(4 lam)). The rectangle
intersection_{p in S} Q(p, rho_S) is clipped to the cell [0,Px) x [0,Py); the area of the union is computed exactly
(integers over (lcm(Dx, Dr))^2); the output lam * A / mu(cell) is a lower bound for c_2 (NOTE_lower_bound.md, Lemma 2).
lam must be 1 here (configs are normalized so).
Usage: python -I check_exact2.py CONFIG.json [--Dr 2^40]
"""
import sys, json, math, argparse
from fractions import Fraction as F


def lcm(a, b):
    return a * b // math.gcd(a, b)


def union_area_int(R):
    """R: list of (x0, x1, y0, y1) integers; returns the integer area of the union."""
    if not R:
        return 0
    ys = sorted({r[2] for r in R} | {r[3] for r in R})
    yi = {v: i for i, v in enumerate(ys)}
    ev = []
    for (x0, x1, y0, y1) in R:
        ev.append((x0, 1, yi[y0], yi[y1])); ev.append((x1, -1, yi[y0], yi[y1]))
    ev.sort()
    n = len(ys) - 1
    cnt = [0] * (4 * n + 4); ln = [0] * (4 * n + 4)

    def upd(node, l, r, a, b, d):
        if b <= l or r <= a:
            return
        if a <= l and r <= b:
            cnt[node] += d
        else:
            m = (l + r) // 2
            upd(2 * node, l, m, a, b, d); upd(2 * node + 1, m, r, a, b, d)
        if cnt[node] > 0:
            ln[node] = ys[r] - ys[l]
        elif r - l == 1:
            ln[node] = 0
        else:
            ln[node] = ln[2 * node] + ln[2 * node + 1]

    area = 0; px = ev[0][0]
    for (x, d, a, b) in ev:
        area += ln[1] * (x - px)
        px = x
        if n > 0:
            upd(1, 0, n, a, b, d)
    return area


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('config'); ap.add_argument('--Dr', type=int, default=2 ** 40)
    a = ap.parse_args()
    cfg = json.load(open(a.config))
    if F(cfg['lam']) != 1:
        raise ValueError('lam must be 1')
    Px, Py = F(cfg['Px']), F(cfg['Py'])
    pts = [(F(x), F(y), F(m)) for (x, y, m) in cfg['pts']]
    for (x, y, m) in pts:
        if not (0 <= x < Px and 0 <= y < Py and m > 0):
            raise ValueError('atoms must lie in [0,Px) x [0,Py) with positive mass')
    reach = F(cfg['reach'])
    Dr = a.Dr
    Dx = 1
    for v in [Px, Py, reach] + [p[0] for p in pts] + [p[1] for p in pts]:
        Dx = lcm(Dx, v.denominator)
    Dm = 1
    for p in pts:
        Dm = lcm(Dm, p[2].denominator)
    L = lcm(Dx, Dr)                       # common denominator of every rectangle coordinate
    sx, sr = L // Dx, L // Dr
    PxI, PyI, reachI = int(Px * Dx), int(Py * Dx), int(reach * Dx)
    P = [(int(x * Dx), int(y * Dx), int(m * Dm)) for (x, y, m) in pts]
    mcellI = sum(p[2] for p in P)
    # rho_max bound: a box of side <= reach meets at most (reach/P + 2)^2 cells
    nb = (int(reach / Px) + 2) * (int(reach / Py) + 2)
    rho_max = F(math.isqrt(nb * mcellI * Dr * Dr // (4 * Dm)) + 1, Dr)
    rmI = int(math.ceil(rho_max * Dx)) + 1
    ext = reachI + rmI
    mx = ext // PxI + 2; my = ext // PyI + 2
    W = [(x + i * PxI, y + j * PyI, m) for i in range(-mx, mx + 1) for j in range(-my, my + 1) for (x, y, m) in P]
    # keep only atoms that can be box bounds or box members: x in [-rmI - reachI, PxI + rmI + reachI] etc.
    W = [w for w in W if -rmI - reachI <= w[0] <= PxI + rmI + reachI and -rmI - reachI <= w[1] <= PyI + rmI + reachI]
    X = sorted({w[0] for w in W}); Y = sorted({w[1] for w in W})
    xi = {v: i for i, v in enumerate(X)}; yi = {v: i for i, v in enumerate(Y)}
    nx, ny = len(X), len(Y)
    G = [[0] * (ny + 1) for _ in range(nx + 1)]
    for (x, y, m) in W:
        G[xi[x] + 1][yi[y] + 1] += m
    for i in range(1, nx + 1):
        row, prev = G[i], G[i - 1]
        acc = 0
        for j in range(1, ny + 1):
            acc += row[j]
            row[j] = prev[j] + acc
    R = []
    four_dm = 4 * Dm
    for ia in range(nx):
        xa = X[ia]
        if xa < -rmI:
            continue
        for ib in range(ia, nx):
            xb = X[ib]
            if xb - xa > reachI or xb > PxI + rmI:
                break
            Gb, Ga = G[ib + 1], G[ia]
            for ja in range(ny):
                ya = Y[ja]
                if ya < -rmI:
                    continue
                for jb in range(ja, ny):
                    yb = Y[jb]
                    if yb - ya > reachI or yb > PyI + rmI:
                        break
                    Ms = Gb[jb + 1] - Ga[jb + 1] - Gb[ja] + Ga[ja]
                    if Ms <= 0:
                        continue
                    k = math.isqrt(Ms * Dr * Dr // four_dm) - 1
                    if k <= 0:
                        continue
                    rho = k * sr                                  # rho_S in units of 1/L
                    x0 = max(xb * sx - rho, 0); x1 = min(xa * sx + rho, PxI * sx)
                    if x1 <= x0:
                        continue
                    y0 = max(yb * sx - rho, 0); y1 = min(ya * sx + rho, PyI * sx)
                    if y1 <= y0:
                        continue
                    R.append((x0, x1, y0, y1))
    Aint = union_area_int(R)
    A = F(Aint, L * L)
    ratio = A / F(mcellI, Dm)
    print('atoms per cell', len(P), 'window atoms', len(W), 'rectangles', len(R))
    print('area in cell %.15f of %.15f' % (float(A), float(Px * Py)))
    print('EXACT ratio lower bound = %s' % (ratio if len(str(ratio)) < 300 else '(%d-digit fraction)' % len(str(ratio))))
    print('ratio >= %s' % (str(math.floor(ratio * 10 ** 10) / 10 ** 10)))
    json.dump(dict(config=a.config, ratio_num=str(ratio.numerator), ratio_den=str(ratio.denominator),
                   ratio_floor_1e10=str(math.floor(ratio * 10 ** 10)), rectangles=len(R), Dr=Dr),
              open(a.config.replace('.json', '') + '.check2.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
