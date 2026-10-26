#!/usr/bin/env python3
"""Finite domain/symmetry/scalar interface audit; never evaluates a heat/Green kernel."""
import time
START = time.monotonic()
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product
from math import comb
from pathlib import Path
import hashlib
import json
import os
import sys

HERE = Path(__file__).resolve().parent
SCALE = 1 << 128


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while block := f.read(8*1024*1024):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(path)


def root(q, k):
    """Direct rational comparison, greedily setting binary root digits."""
    assert q > 0
    n, d = q.numerator*SCALE**k, q.denominator
    top = max(0, (n.bit_length()-d.bit_length())//k+2)
    a = 0
    for bit in range(top, -1, -1):
        candidate = a+(1 << bit)
        if candidate**k*d <= n:
            a = candidate
    low, high = Q(a, SCALE), Q(a+1, SCALE)
    assert low**k <= q < high**k
    return low, high


def decimal_outward(q, places, upper):
    scale = 10**places
    n = q.numerator*scale//q.denominator + int(upper)
    digits = str(n).zfill(places+1)
    return digits[:-places]+'.'+digits[-places:]


def main():
    reg = json.loads((HERE/'registration.json').read_text())
    assert reg['status'] == 'registered_before_single_finite_interface_batch'
    assert reg['source_sha256'] == sha(__file__)
    assert not (HERE/'results.json').exists(), 'No repeat of a completed or failed batch'
    out = {'pid': os.getpid(), 'python': sys.executable, 'status': 'running',
           'registration_sha256': sha(HERE/'registration.json'), 'source_sha256': sha(__file__)}
    checkpoint = {'pid': os.getpid(), 'status': 'input_hashes', 'completed_roots': 0,
                  'completed_ordered_grid_vertices': 0, 'completed_generators': 0}
    write(HERE/'checkpoint.json', checkpoint)
    print(json.dumps(out), flush=True)
    try:
        for item in reg['inputs']:
            assert sha(item['path']) == item['sha256'], item['path']
        # Only the owner manifest and scalar auxiliary input are parsed. Never parse final values.
        owner = json.loads(Path(reg['roles']['owner_manifest']).read_text())
        assert owner['roots'] == 52360 and owner['passed'] == 97960
        assert owner['split'] == 45600 and owner['queue'] == owner['inflight'] == 0
        assert owner['passed'] == owner['roots'] + owner['split']
        expected_roots = sorted(tuple(reversed(a)) for a in combinations_with_replacement(range(32),4))
        assert len(expected_roots) == comb(35,4) == 52360
        assert len(set(expected_roots)) == len(expected_roots)
        for n, r in enumerate(expected_roots, 1):
            assert 31 >= r[0] >= r[1] >= r[2] >= r[3] >= 0
            if n % 8192 == 0:
                checkpoint.update(status='finite_root_domain', completed_roots=n)
                write(HERE/'checkpoint.json', checkpoint)
        roots_digest = hashlib.sha256(json.dumps(expected_roots,separators=(',',':')).encode()).hexdigest()
        rootset = set(expected_roots)
        count = 0
        for a in combinations_with_replacement(range(33),4):
            vertex = tuple(reversed(a))
            chosen = tuple(min(x,31) for x in vertex)
            assert chosen in rootset
            assert all(i <= x <= i+1 for i,x in zip(chosen,vertex))
            count += 1
            if count % 8192 == 0:
                checkpoint.update(status='ordered_grid_vertex_endpoints', completed_roots=len(expected_roots),
                                  completed_ordered_grid_vertices=count)
                write(HERE/'checkpoint.json', checkpoint)
        assert count == comb(36,4) == 58905
        V = [(1,)+s for s in product((-1,1),repeat=4)]
        def canonical(v):
            return tuple(v[0]*x for x in v)
        generators = []
        for j in range(5):
            images = [tuple(-x if i == j else x for i,x in enumerate(v)) for v in V]
            generators.append((f'coordinate_sign_flip_{j+1}',images))
        for j in range(4):
            images = []
            for v in V:
                w = list(v); w[j],w[j+1] = w[j+1],w[j]
                images.append(tuple(w))
            generators.append((f'adjacent_coordinate_swap_{j+1}_{j+2}',images))
        symmetry = []
        for k,(name,images) in enumerate(generators,1):
            reps = [canonical(v) for v in images]
            assert set(reps) == set(V) and len(set(reps)) == 16
            symmetry.append({'generator':name, 'representative_permutation':[V.index(v) for v in reps],
                             'representative_signs':[v[0] for v in images]})
            checkpoint.update(status='signed_permutation_generators', completed_generators=k)
            write(HERE/'checkpoint.json', checkpoint)
        budget = json.loads(Path(reg['roles']['orbit_budget']).read_text())
        I = Q(budget['new_I_strict_simple_rational_cap'])
        assert I == Q(71555899,10**12)
        g = Q(1,10000); h = g/5; N = 32*(g-h)
        assert h == Q(1,50000) and N == Q(8,3125)
        cost4 = 3125*I**4/g**5
        assert cost4 > 32**4
        cost = root(cost4,4)
        threshold_I4 = Q(32**4)*g**5/3125
        threshold_I = root(threshold_I4,4)
        threshold_g5 = 3125*I**4/Q(32**4)
        threshold_g = root(threshold_g5,5)
        assert threshold_I[1] < I and g < threshold_g[0]
        oldI = Q(36112162239,500000000000000)
        assert I < oldI
        # A concrete smaller I budget would be sufficient only after its own full-integral proof.
        proposed_I = Q(4,100000)
        proposed_cost4 = 3125*proposed_I**4/g**5
        assert proposed_cost4 < 32**4
        proposed_cost = root(proposed_cost4,4)
        intervals = {'paired_cost':cost,'I_threshold_for_cost32':threshold_I,
                     'g_threshold_for_current_Icap':threshold_g,'conditional_cost_if_I_le_4e_minus5':proposed_cost}
        out.update({'status':'conditional_numeric_gate', 'interface_algebra_checks':'passed',
                    'owner_counts_are_unaccepted_claims_only': {'roots':owner['roots'],'passed':owner['passed'],
                                                             'split':owner['split'],'queue':0,'inflight':0},
                    'expected_root_count':len(expected_roots),'expected_root_canonical_sha256':roots_digest,
                    'actual_final_root_set_parsed_or_certified':False,
                    'ordered_grid_vertices_checked':count,'symmetry_generators':symmetry,
                    'Icap':str(I),'g_candidate':str(g),'h':str(h),'N':str(N),'paired_cost_fourth_power':str(cost4),
                    'cost_fourth_power_minus_32_fourth_power':str(cost4-32**4),
                    'threshold_I_fourth_power':str(threshold_I4),'threshold_g_fifth_power':str(threshold_g5),
                    'intervals':{name:{'lower':str(x[0]),'upper':str(x[1]),
                        'decimal_outward':[decimal_outward(x[0],20,False),decimal_outward(x[1],20,True)]}
                        for name,x in intervals.items()},
                    'new_kernel_or_special_function_calls':0,'numerical_gate081_passed_here':False,
                    'global_g_certified_here':False,'new_formal_C5_bound':False,'formal_C5_upper':32})
        for item in reg['inputs']:
            assert sha(item['path']) == item['sha256'], item['path']
        checkpoint.update(status='complete_conditional_numeric_gate',completed_roots=len(expected_roots),
                          completed_ordered_grid_vertices=count, completed_generators=len(generators))
    except Exception as error:
        out.update(status='stopped_failure_no_rerun',error=repr(error))
        checkpoint.update(status='failed',error=repr(error))
        raise
    finally:
        out['seconds_before_final_write'] = time.monotonic()-START
        write(HERE/'results.json',out)
        write(HERE/'checkpoint.json',checkpoint)
    print(json.dumps({'status':out['status'],'pid':out['pid'],'seconds':out['seconds_before_final_write'],
                      'intervals':out['intervals']}),flush=True)


if __name__ == '__main__':
    main()
