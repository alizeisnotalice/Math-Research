from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import json, hashlib, math

ROOT = Path(__file__).parent
PREFIX = 'original_ordered_martingale_dilation'
DATE = '20261007'
checks = 0

def check(condition):
    global checks
    assert condition
    checks += 1

def mean(v):
    return sum(v, F(0)) / len(v)

def conditional(v, j):
    width = 2 ** j
    out = []
    for start in range(0, len(v), width):
        avg = mean(v[start:start+width])
        out.extend([avg] * width)
    return out

def packed(v):
    return [str(x) for x in v]

rows = []
for n, J in [(4, 3), (16, 6), (64, 9)]:
    before = checks
    qs = [F(x) for x in [0, 1, 3, 31, 255]]
    rs = [F(0), F(1, n+1), F(1, 3), F(2, 3), F(1)]
    symbol_checks = []
    for r, s in combinations(rs, 2):
        A, B = 1-r, 1-s
        for q in qs:
            for k in range(1, 7):
                # (-1)^(k-1) times kth derivative of the log increment.
                positive_derivative = F(math.factorial(k-1), 2) * (
                    (A/(1+A*q))**k - (B/(1+B*q))**k)
                check(positive_derivative >= 0)
            kr = (1+(1-r)*q)/(1+q)
            ks = (1+(1-s)*q)/(1+q)
            check(0 < ks/kr <= 1)
        symbol_checks.append({'r':str(r), 's':str(s)})
    for r, s, t in combinations(rs, 3):
        for q in qs:
            k = lambda v: (1+(1-v)*q)/(1+q)
            check((k(s)/k(r))*(k(t)/k(s)) == k(t)/k(r))
    tensor_q = [qs[i % len(qs)] for i in range(n)]
    for r in rs:
        factors = [(1+(1-r)*q)/(1+q) for q in tensor_q]
        K = math.prod(factors)
        direct = sum((-tensor_q[i]/(1+tensor_q[i])) *
                     math.prod(factors[:i]+factors[i+1:]) for i in range(n))
        via_log = -K * sum((q/(1+(1-r)*q) for q in tensor_q), F(0))
        check(direct == via_log)

    count = 2 ** (J+1)
    # Algebra-only finite filtration; not an original G_c spatial discretization.
    f = [F(1+x.bit_count()) + F(x % 11, 3) if x % 2 else F(0)
         for x in range(count)]
    labels = [1+(x//2) % J if x % 2 == 0 else 0 for x in range(count)]
    a = {j:[F(labels[x] == j) for x in range(count)] for j in range(1,J+1)}
    M = {j:conditional(f,j) for j in range(1,J+1)}
    pi = {j:conditional(a[j],j) for j in range(1,J+1)}
    b = {j:[F(0)]*count for j in range(1,J+1)}
    survive = [F(1)]*count
    for j in range(J, 0, -1):
        b[j] = [pi[j][x]*survive[x] for x in range(count)]
        survive = [survive[x]*(1-pi[j][x]) for x in range(count)]
        check(conditional(b[j],j) == b[j])
        check(mean([b[j][x]*M[j][x] for x in range(count)]) ==
              mean([b[j][x]*f[x] for x in range(count)]))
    total_pi = [sum((pi[j][x] for j in range(1,J+1)),F(0)) for x in range(count)]
    total_b = [sum((b[j][x] for j in range(1,J+1)),F(0)) for x in range(count)]
    for x in range(count):
        check(sum(a[j][x] for j in range(1,J+1)) <= 1)
        check(0 <= total_b[x] <= 1)
        check(total_b[x] == 1-survive[x])
        check(total_pi[x]-total_b[x] >= 0)
    T = sum(mean([a[j][x]*M[j][x] for x in range(count)]) for j in range(1,J+1))
    check(T == mean([f[x]*total_pi[x] for x in range(count)]))
    cov_sum = F(0)
    for j in range(1,J+1):
        af = [a[j][x]*f[x] for x in range(count)]
        check(all(x == 0 for x in af))
        conditional_af = conditional(af,j)
        cov = [conditional_af[x]-pi[j][x]*M[j][x] for x in range(count)]
        cov_sum += mean(cov)
    check(T == -cov_sum)
    W = mean(f)
    stopped = mean([f[x]*total_b[x] for x in range(count)])
    deficit = mean([f[x]*(total_pi[x]-1+survive[x]) for x in range(count)])
    pair = mean([f[x]*sum((pi[i][x]*pi[j][x] for i,j in combinations(range(1,J+1),2)),F(0)) for x in range(count)])
    check(0 <= stopped <= W)
    check(T == stopped+deficit)
    check(0 <= deficit <= pair)
    rows.append({'tensor_n':n,'filtration_J':J,'terminal_atoms':count,
                 'checks':checks-before,'spectral_pair_nodes':symbol_checks,
                 'W':str(W),'selected_traffic':str(T),'stopped_traffic':str(stopped),
                 'exact_deficit':str(deficit),'pair_overlap_upper':str(pair),
                 'complete_terminal_f':packed(f),'receiver_labels':labels,
                 'M':{str(j):packed(M[j]) for j in M},
                 'pi':{str(j):packed(pi[j]) for j in pi},
                 'b':{str(j):packed(b[j]) for j in b}})

reg = ROOT / f'{PREFIX}_registration_{DATE}.json'
payload = {'status':'PASS','checks':checks,'rounds':rows,
           'scope':'Exact original-symbol algebra and universal finite-filtration identities only; finite filtration is not an original-kernel, obstacle, or actual geom input.',
           'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),
           'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'arithmetic':'Fraction exact; no tolerance, Monte Carlo, or confidence statements'}
dest = ROOT / f'{PREFIX}_results_{DATE}.json'
assert not dest.exists(), 'Never overwrite a previous result'
dest.write_text(json.dumps(payload,indent=2))
print(json.dumps({'status':'PASS','checks':checks,'rounds':[{'n':r['tensor_n'],'J':r['filtration_J'],'W':r['W'],'selected_traffic':r['selected_traffic'],'stopped_traffic':r['stopped_traffic'],'exact_deficit':r['exact_deficit']} for r in rows]}))
