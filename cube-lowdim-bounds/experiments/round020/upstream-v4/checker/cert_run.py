"""Driver for the hl47a certificate checks (search-side verifier). Usage:
   python -I cert_run.py SPEC.json [--Q 8] [--maxdepth 4] [--procs 24] [--only K|A|M]
Checks (C1) K >= 0, (C2) K >= 1 on D, (C3) mass, (C4) A K >= 0 on the fundamental domain 0 <= v <= u (D4 symmetry).
"""
import sys, os, json, math, time, argparse
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, fmpq, ctx
import hlcert
from hlcert import Kern, Q

ctx.prec = 96
SPEC = None
KER = None


def _init(spec):
    global KER
    ctx.prec = 96
    KER = Kern(spec)


# ------------------------------------------------------------------------------------------ K checks
def tri_bb_min(Pw):
    """Pw: power coefficients c[i][j] (fmpq) of a bicubic in (tau, sigma). Return min BB coefficient of degree 6 on the
    triangle tau, sigma >= 0, tau + sigma <= 1 (exact)."""
    n = 6
    # tau^i sigma^j (lam0+lam1+lam2)^{n-i-j}, lam1 = tau, lam2 = sigma: coefficient of lam0^a lam1^b lam2^c
    from math import comb, factorial
    bb = {}
    for i in range(4):
        for j in range(4):
            cij = Pw[i][j]
            if cij == 0:
                continue
            m = n - i - j
            # (lam0+lam1+lam2)^m = sum_{a+b+c=m} m!/(a!b!c!) lam0^a lam1^b lam2^c
            for a in range(m + 1):
                for b in range(m - a + 1):
                    c = m - a - b
                    coef = fmpq(factorial(m), factorial(a) * factorial(b) * factorial(c))
                    key = (a, b + i, c + j)
                    bb[key] = bb.get(key, fmpq(0)) + cij * coef
    # BB coefficient = (monomial coefficient) / (n!/(a!b!c!))
    out = []
    for a in range(n + 1):
        for b in range(n - a + 1):
            c = n - a - b
            v = bb.get((a, b, c), fmpq(0))
            out.append(v / fmpq(math.factorial(n), math.factorial(a) * math.factorial(b) * math.factorial(c)))
    return min(out)


def bicubic_power(u0, u1, v0, v1):
    """power coefficients in local (tau, sigma) of S on the cell (exact)."""
    K = KER
    h, N = K.h_q, K.N
    s0, s1, q0, q1 = u0 / h, u1 / h, v0 / h, v1 / h
    k = int(math.floor(float(s0) + 1e-12)); l = int(math.floor(float(q0) + 1e-12))
    # power coefficients of B(x - off) in local x of the knot interval, then substitute x = t0 + dt * tau
    def pw(off, t0, t1):
        m = -off
        p = [fmpq(c, 6) for c in hlcert._BP[m]]
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
            aij = K.A[i + N][j + N]
            if aij == 0:
                continue
            for x in range(4):
                for y in range(4):
                    C[x][y] += aij * pu[x] * pv[y]
    return C


def K_cell_job(cell):
    """returns list of failures for the knot cell (with recursive subdivision)."""
    (u0, u1, v0, v1, depth, maxdepth) = cell
    K = KER
    one = fmpq(1)
    rmin, rmax = u0 + v0, u1 + v1
    if rmax <= one:
        lb = K.K_bounds(u0, u1, v0, v1); ok = lb >= 1; target = '1'
    elif rmin >= one:
        lb = K.K_bounds(u0, u1, v0, v1); ok = lb >= 0; target = '0'
    else:
        lb0 = K.K_bounds(u0, u1, v0, v1)
        ok = False; lb = lb0; target = '1 on D-part'
        if u0 + v0 + (u1 - u0) == one and (u1 - u0) == (v1 - v0) and lb0 >= 0:
            lbt = K.K_bounds(u0, u1, v0, v1, triangle=True)
            ok = lbt >= 1; lb = lbt
    if ok:
        return []
    if depth < maxdepth:
        um, vm = (u0 + u1) / 2, (v0 + v1) / 2
        fails = []
        for (a, b, c_, d) in ((u0, um, v0, vm), (um, u1, v0, vm), (u0, um, vm, v1), (um, u1, vm, v1)):
            fails += K_cell_job((a, b, c_, d, depth + 1, maxdepth))
        return fails
    return [(str(u0), str(u1), str(v0), str(v1), str(lb), target)]


# ------------------------------------------------------------------------------------------ A K checks
def A_cell(u0, u1, v0, v1):
    K = KER
    return K.A_lower(u0, u1, v0, v1)


def cell_has_support(u0, u1, v0, v1):
    """False if K vanishes identically on the cell (then A K >= 0 there automatically, given K >= 0)."""
    K = KER
    if u0 + v0 < K.R_q:
        return True
    for comp in K.comps:
        h, N = comp.h_q, comp.N
        i_lo = int(math.floor(float(u0 / h))) - 1; i_hi = int(math.floor(float(u1 / h))) + 2
        j_lo = int(math.floor(float(v0 / h))) - 1; j_hi = int(math.floor(float(v1 / h))) + 2
        for i in range(max(i_lo, -N), min(i_hi, N) + 1):
            for j in range(max(j_lo, -N), min(j_hi, N) + 1):
                if comp.A[i + N][j + N] != 0:
                    return True
    return False


def A_col_job(args):
    """all cells of one column (u-interval) of size delta, with recursive refinement on failure."""
    (m, delta, maxdepth, nmax) = args
    u0, u1 = m * delta, (m + 1) * delta
    fails = []; nchecked = 0; minlb = None
    for n in range(0, m + 1):
        v0, v1 = n * delta, (n + 1) * delta
        if not cell_has_support(u0, u1, v0, v1):
            continue
        res, cnt, mlb = _A_rec(u0, u1, v0, v1, 0, maxdepth)
        fails += res; nchecked += cnt
        minlb = mlb if minlb is None else min(minlb, mlb)
    return fails, nchecked, minlb


def _A_rec(u0, u1, v0, v1, depth, maxdepth):
    lb = A_cell(u0, u1, v0, v1)
    lbl = float(lb.lower())
    if lb.lower() >= 0:
        return [], 1, lbl
    if depth < maxdepth:
        um, vm = (u0 + u1) / 2, (v0 + v1) / 2
        fails = []; cnt = 1; mlb = None
        for (a, b, c_, d) in ((u0, um, v0, vm), (um, u1, v0, vm), (u0, um, vm, v1), (um, u1, vm, v1)):
            if c_ > b:          # entirely above the diagonal v = u: outside the fundamental domain
                continue
            r, c2, m2 = _A_rec(a, b, c_, d, depth + 1, maxdepth)
            fails += r; cnt += c2
            mlb = m2 if mlb is None else min(mlb, m2)
        return fails, cnt, (mlb if mlb is not None else lbl)
    return [(str(u0), str(u1), str(v0), str(v1), lbl)], 1, lbl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('spec')
    ap.add_argument('--Q', type=int, default=8)
    ap.add_argument('--maxdepth', type=int, default=3)
    ap.add_argument('--Kdepth', type=int, default=6)
    ap.add_argument('--procs', type=int, default=24)
    ap.add_argument('--only', type=str, default='MKA')
    ap.add_argument('--out', type=str, default='')
    a = ap.parse_args()
    spec = json.load(open(a.spec))
    _init(spec)
    K = KER
    t0 = time.time()
    report = {}
    if 'M' in a.only:
        m = K.mass_upper()
        report['mass_half'] = str(m)
        print('(C3) (1/2) int K in', m, ' upper', float(m.upper()), flush=True)
    # support box over all spline components; base cell size = finest knot spacing
    U = fmpq(0)
    hmin = None
    for comp_spec in [spec] + spec.get('extra', []):
        hc = Q(comp_spec['h'])
        hmin = hc if hmin is None else min(hmin, hc)
        if comp_spec['a']:
            imax = max(max(abs(int(k.split(',')[0])), abs(int(k.split(',')[1]))) for k in comp_spec['a'])
            U = max(U, (imax + 2) * hc)
    h = hmin
    if not all((Q(cs['h']) / h).q == 1 for cs in [spec] + spec.get('extra', [])):
        raise ValueError('knot spacings must be multiples of the finest')
    N = K.N
    report['U'] = str(U)
    if 'K' in a.only:
        cells = []
        n_k = int(math.ceil(float(U / h)))
        for i in range(n_k):
            for j in range(i + 1):
                cells.append((i * h, (i + 1) * h, j * h, (j + 1) * h, 0, a.Kdepth))
        with Pool(a.procs, initializer=_init, initargs=(spec,)) as pool:
            res = pool.map(K_cell_job, cells, chunksize=4)
        kf = [f for r in res for f in r]
        report['K_fail'] = kf
        print(f'(C1/C2) K checks: {len(cells)} knot cells, failures {len(kf)}  t={time.time()-t0:.1f}', flush=True)
        for f in kf[:10]:
            print('   ', f)
    if 'A' in a.only:
        delta = h / a.Q
        M = int(math.ceil(float(U / delta)))
        jobs = [(m, delta, a.maxdepth, M) for m in range(M)]
        with Pool(a.procs, initializer=_init, initargs=(spec,)) as pool:
            res = pool.map(A_col_job, jobs, chunksize=1)
        af = [f for r in res for f in r[0]]
        nch = sum(r[1] for r in res)
        mins = [r[2] for r in res if r[2] is not None]
        report['A_fail'] = af; report['A_cells'] = nch
        print(f'(C4) A K checks: delta = h/{a.Q}, cells evaluated {nch}, failures {len(af)}, '
              f'min cell lower bound {min(mins) if mins else None:.4g}  t={time.time()-t0:.1f}', flush=True)
        af.sort(key=lambda f: f[4])
        for f in af[:15]:
            print('   ', f)
    if a.out:
        json.dump(report, open(a.out, 'w'), indent=1)
    nfail = len(report.get('K_fail', [])) + len(report.get('A_fail', []))
    print('PASS' if nfail == 0 else f'FAIL ({nfail} failing cells)', flush=True)
    if nfail:
        sys.exit(1)


if __name__ == '__main__':
    main()
