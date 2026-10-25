#!/usr/bin/env python3
"""Finite 560-triple, two-linear-factor completion batch; no 5-basis scan."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import sys
import time

import flint
from flint import arb, ctx

HERE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def det3(a):
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])


def adj3(a):
    ans = [[0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rows = [[a[u][v] for v in range(3) if v != i] for u in range(3) if u != j]
            ans[i][j] = (-1)**(i+j)*(rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0])
    return ans


def det_fraction(rows):
    a = [list(map(F, row)) for row in rows]
    result = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F(0)
        if j != i:
            a[j], a[i] = a[i], a[j]
            result = -result
        pivot = a[i][i]
        result *= pivot
        for j in range(i+1, len(a)):
            ratio = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
            a[j][i] = 0
    return result


def ar(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def parse(s):
    center, radius = map(F, s[1:-1].split(" +/- "))
    return center-radius, center+radius


def hist_key(coefficients):
    total = sum(coefficients)
    assert total > 0
    return tuple(sorted(Counter(F(x, total) for x in coefficients if x).items()))


def main():
    began = time.monotonic()
    registration = json.loads((HERE/'registration.json').read_text())
    for row in registration['inputs']:
        assert sha(Path(row['source'])) == row['sha256']
    old_path = next(Path(r['source']) for r in registration['inputs'] if r['source'].endswith('upper_I020/results.json'))
    old = json.loads(old_path.read_text())
    old_J = parse(old['arb_384bit']['fractional_Holder_J_upper'])
    old_I = parse(old['arb_384bit']['fractional_Holder_complete_I_upper'])
    V = [(1,)+s for s in product((-1,1),repeat=4)]
    VV = [[dot(v,w) for w in V] for v in V]
    W, denominator = 16**5, 20*16**5
    type_groups = {}
    scope_occupancy = [0]*16
    triple_types = Counter()
    total_C = 0
    total_quad_C = 0
    completion_example = None
    geometry_hash = hashlib.sha256()
    group_count = 0

    for J in combinations(range(16),3):
        gram = [[VV[i][j] for j in J] for i in J]
        D = det3(gram)
        assert D > 0
        adj = adj3(gram)
        assert [[sum(gram[i][k]*adj[k][j] for k in range(3)) for j in range(3)] for i in range(3)] == [[D*int(i==j) for j in range(3)] for i in range(3)]
        z = [[VV[j][i] for j in J] for i in range(16)]
        za = [[sum(z[i][k]*adj[k][j] for k in range(3)) for j in range(3)] for i in range(16)]
        A = [[D*VV[i][j]-dot(za[i],z[j]) for j in range(16)] for i in range(16)]
        assert all(A[i][i] >= 0 for i in range(16))
        nonzero = {i for i in range(16) if A[i][i]}
        assert not (set(J) & nonzero)
        weights = [[0]*16 for _ in range(16)]
        for i in range(16):
            for k in range(i+1,16):
                numerator = A[i][i]*A[k][k]-A[i][k]**2
                assert numerator >= 0 and numerator % D == 0
                weights[i][k] = weights[k][i] = numerator//D
        C_J = sum(weights[i][k] for i in range(16) for k in range(i+1,16))
        assert C_J == D*16**2
        total_quad_C += C_J
        classes = []
        unused = set(nonzero)
        while unused:
            i = min(unused)
            line = {k for k in unused if weights[i][k] == 0}
            assert all(weights[u][v] == 0 for u in line for v in line)
            assert all(weights[u][v] > 0 for u in line for v in nonzero-line)
            classes.append(line)
            unused -= line
        scope = set(J) | nonzero
        triple_types[(D,len(nonzero),len(classes))] += 1
        local_C = 0
        for line in classes:
            i = min(line)
            outside = nonzero-line
            sum_diag = sum(A[u][u] for u in line)
            sum_L = sum(weights[i][k] for k in outside)
            Cg = sum(weights[u][k] for u in line for k in outside)
            assert Cg == F(sum_diag,A[i][i])*sum_L
            assert Cg > 0
            # Verify the rank-one coefficient factorization along the entire line.
            assert all(weights[u][k]*A[i][i] == A[u][u]*weights[i][k] for u in line for k in outside)
            X_hist = hist_key([A[u][u] for u in line])
            L_hist = hist_key([weights[i][k] for k in outside])
            pair = tuple(sorted((X_hist,L_hist)))
            key = (len(scope),pair)
            if key not in type_groups:
                type_groups[key] = {'scope_size':len(scope),'probability_histograms':pair,'sum_C':0,'count':0,
                                    'representative_J':J,'representative_line':sorted(line)}
            row = type_groups[key]
            row['sum_C'] += Cg
            row['count'] += 1
            local_C += Cg
            total_C += Cg
            group_count += 1
            for u in scope:
                scope_occupancy[u] += Cg
            geometry_hash.update(json.dumps([J,sorted(line),Cg,str(X_hist),str(L_hist)],separators=(',',':')).encode())
        assert local_C == 2*C_J
        if completion_example is None:
            completion_example = (J,D,weights)

    assert total_quad_C == 10*W and total_C == denominator
    assert len(set(scope_occupancy)) == 1
    rho = F(scope_occupancy[0],denominator)
    assert 0 < rho <= 1
    assert sum(F(g['sum_C'],denominator) for g in type_groups.values()) == 1
    assert sum(F(g['sum_C']*g['scope_size'],16*denominator) for g in type_groups.values()) == rho
    J0,D0,w0 = completion_example
    for i,k in combinations(range(16),2):
        assert det_fraction([V[j] for j in J0]+[V[i],V[k]])**2 == w0[i][k]

    beta,b,q = F(3,4),F(3,20),F(10,7)
    ctx.prec = 384
    R_product = arb(1)
    types_out = []
    for key,row in sorted(type_groups.items()):
        theta = F(row['sum_C'],denominator)
        R_pair = []
        for hist in row['probability_histograms']:
            R_pair.append(sum((count*ar(prob)**ar(beta) for prob,count in hist),arb(0)))
        Rg = R_pair[0]*R_pair[1]
        R_product *= Rg**ar(-theta*b/beta)
        types_out.append({'scope_size':row['scope_size'],'group_count':row['count'],
                          'sum_C':row['sum_C'],'theta':str(theta),
                          'probability_histograms':[{str(prob):count for prob,count in hist} for hist in row['probability_histograms']],
                          'representative_J':row['representative_J'],'representative_line':row['representative_line'],
                          'R_product_ball':Rg.str(80)})
    moment = ar(1+b*rho/beta).gamma()/ar(1+b*rho).gamma()
    Jnew = 2*ar(F(7,4)).gamma()*moment**ar(F(5)/rho)*R_product
    c = 3/(ar(2)**ar(F(15,2))*arb.pi()**ar(F(5,2)))
    Inew = c**ar(q)*arb.pi()**2/60*Jnew**ar(q)
    upper = F(str(Inew.upper().fmpq()))
    cap = F((upper*10**13).__ceil__(),10**13)
    assert Inew < ar(cap)
    checks = {'J_improves_locked020':bool(Jnew<ar(old_J[0])),
              'I_improves_locked020':bool(Inew<ar(old_I[0])),
              'I_below_target_11e_minus6':bool(Inew<ar(F(11,10**6)))}
    out = {'status':'finite_three_subset_two_level_batch_complete','new_C5_upper_bound':False,
           'script_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),
           'source_inputs':registration['inputs'],
           'geometry':{'triple_count':560,'two_level_group_count':group_count,
                       'completion_weight_sum':total_quad_C,'group_weight_sum':total_C,
                       'expected_20W':denominator,'coordinate_scope_rho':str(rho),
                       'triple_types':[{ 'gram_determinant':D,'nonzero_quotient_directions':m,
                                        'quotient_line_classes':l,'count':count}
                                       for (D,m,l),count in sorted(triple_types.items())],
                       'probability_type_count':len(type_groups),'probability_types':types_out,
                       'geometry_stream_sha256':geometry_hash.hexdigest(),
                       'independent_5x5_completion_example':{'J':J0,'Gram_det':D0,'pair_checks':120}},
           'parameters':{'beta':str(beta),'b':str(b),'rho':str(rho),'q':str(q)},
           'arb_384bit':{'two_linear_factor_scale_gain':R_product.str(85),
                        'ES_minus_b_rho':moment.str(85),'J_upper':Jnew.str(85),
                        'complete_I_upper':Inew.str(85),'new_over_locked020_I_lower':(Inew/ar(old_I[0])).str(85)},
           'strict_rational_I_upper':str(cap),'checks':checks,
           'runtime':{'python':sys.executable,'flint':flint.__version__,'precision_bits':ctx.prec,
                      'seconds':time.monotonic()-began,'basis_5subset_scan':False,'spatial_or_probability_samples':0},
           'scope':'New full-time/full-face auxiliary I mechanism pending independent audit; no type-weight optimization and no new g/C5 claim'}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'rho':str(rho),'groups':group_count,
                      'probability_types':len(type_groups),'J':Jnew.str(28),'I':Inew.str(28),
                      'rational_I_upper':str(cap),'checks':checks,'seconds':out['runtime']['seconds']},indent=2))


if __name__=='__main__':
    main()
