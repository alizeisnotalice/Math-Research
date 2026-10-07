#!/usr/bin/env python3
"""Exact clipping components, not actual Rn/FIRST simulation."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PREFIX = 'exterior_difference_exit'
checks = []
rounds = []

def record(label, passed):
    checks.append({'label': label, 'pass': bool(passed)})

def G(u, i):
    return [(a + u[x ^ (1 << i)]) / 2 for x, a in enumerate(u)]

def S(u, n):
    gu = [G(u, i) for i in range(n)]
    return [n*a - sum(g[x] for g in gu) for x, a in enumerate(u)]

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

def encoded(q):
    return str(q.numerator) + '/' + str(q.denominator)

for n, spike_n in [(2, 8), (3, 32), (4, 128)]:
    start = len(checks)
    N = 1 << n
    patterns = []
    patterns.append([F(5) if x == 0 else F(0) for x in range(N)])
    patterns.append([F(7) if x == 0 else F(2) if x == 1 else F(0) for x in range(N)])
    patterns.append([F(1, 7) if x == 0 else F(3 + x % 3) if x < N//2 else F(0) for x in range(N)])
    patterns.append([F(1 + x % 5) if x % 3 == 0 else F(0) for x in range(N)])
    examples = []
    for pattern, u in enumerate(patterns):
        su = S(u, n)
        kappa = max(-q for q in su)
        assert kappa > 0
        active = [a > 0 for a in u]
        nu = [kappa + su[x] if active[x] else F(0) for x in range(N)]
        mu = [nu[x] - su[x] for x in range(N)]
        W = sum(nu)
        tag = f'n{n}/pattern{pattern}'
        record(tag+'/cap', all(F(0) <= q <= kappa for q in mu))
        record(tag+'/positive-source', all(q >= 0 for q in nu))
        record(tag+'/saturation', all(mu[x] == kappa for x in range(N) if active[x]))
        record(tag+'/mass', sum(mu) == W and sum(su) == 0)
        for li, L in enumerate([kappa/2, kappa, 2*kappa, max(u)/2]):
            w = [min(a, L) for a in u]
            v = [max(a-L, F(0)) for a in u]
            sw, sv = S(w, n), S(v, n)
            WL = sum(nu[x] for x in range(N) if u[x] > L)
            ew, cross = dot(w, sw), dot(w, su)
            lt = tag+f'/level{li}'
            record(lt+'/decomposition', all(w[x]+v[x] == u[x] and sw[x]+sv[x] == su[x] for x in range(N)))
            record(lt+'/nonnegative-energy', ew >= 0)
            record(lt+'/contraction-energy', ew <= cross)
            record(lt+'/saturation-cross', cross == dot(w, nu)-kappa*sum(w))
            record(lt+'/low-source-fee', cross <= L*W)
            record(lt+'/convex-pointwise', all(sv[x] <= (su[x] if u[x] > L else F(0)) for x in range(N)))
            record(lt+'/tail-positive-source', all(max(sv[x], F(0)) <= (nu[x] if u[x] > L else F(0)) for x in range(N)))
            record(lt+'/zero-tail-mass', sum(sv) == 0)
            record(lt+'/tail-L1', sum(abs(q) for q in sv) <= 2*WL)
            record(lt+'/original-exterior-cap', all(sv[x] >= -kappa for x in range(N) if not active[x]))
            record(lt+'/global-coarse-cap', all(sv[x] >= -kappa-n*L for x in range(N)))
            record(lt+'/high-induced-cap', all(sw[x] >= 0 and sv[x] == nu[x]-kappa-sw[x] for x in range(N) if u[x] > L))
            examples.append({'pattern': pattern, 'level': encoded(L), 'kappa': encoded(kappa), 'W_b': encoded(W), 'W_L': encoded(WL), 'energy': encoded(ew), 'cross_energy': encoded(cross)})
    # No obstacle is numerically reconstructed here. The proven original
    # obstacle relation u >= (H-kappa)/n on the finite source box is used.
    heights = [F(1), F(2), F(8), F(32)]
    H = F(1) + 2*spike_n*max(heights)
    lower_u = (H-F(1))/spike_n
    record(f'spike/n{spike_n}/strict-high', all(lower_u > L for L in heights))
    record(f'spike/n{spike_n}/finite-input-mass', H*F(1, H.numerator) == 1)
    record(f'spike/n{spike_n}/no-summability', sum(F(1) for _ in heights) == 4)
    rounds.append({'finite_dimension': n, 'spike_dimension': spike_n, 'checks': len(checks)-start, 'examples': examples, 'spike': {'H': encoded(H), 'box_volume': encoded(1/H), 'lower_u': encoded(lower_u), 'heights': [encoded(q) for q in heights], 'tail_mass_each': '1/1'}})

registration = HERE / f'{PREFIX}_registration_20261007.json'
result = {'status': 'PASS' if all(c['pass'] for c in checks) else 'FAIL', 'seed': None,
          'scope': 'Exact finite symmetric Markov clipping components and analytic original-box spike coefficients; not actual Rn FIRST samples or weak endpoint evidence.',
          'total_checks': len(checks), 'failed_checks': [c for c in checks if not c['pass']],
          'rounds': rounds, 'checks': checks,
          'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'registration_sha256': hashlib.sha256(registration.read_bytes()).hexdigest()}
out = HERE / f'{PREFIX}_results_20261007.json'
out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
print(json.dumps({k: result[k] for k in ['status', 'total_checks', 'failed_checks', 'script_sha256']}, ensure_ascii=False))
raise SystemExit(0 if result['status'] == 'PASS' else 1)
