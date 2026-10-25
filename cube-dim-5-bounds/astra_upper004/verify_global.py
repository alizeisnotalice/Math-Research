"""Exact symmetric space polynomials and a full positive-face Bernstein certificate."""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations, product
from math import comb, factorial
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parent
RAW = Path('/tmp/astra_upper004_graph_moments.json')


def orbit(partition):
    return set(permutations(partition))


def read_family(raw, denominators):
    out = []
    for poly in raw:
        result = {}
        for key, row in poly.items():
            exponent = tuple(map(int, key.split(',')))
            value = sum((F(int(n), denominators[i]) for i, n in enumerate(row)), F(0))
            if value:
                result[exponent] = value
        out.append(result)
    return out


def corner_value(poly):
    return sum((value * len(orbit(key)) for key, value in poly.items()), F(0))


def bernstein_face(power):
    # Substitute X1=1. Symmetry and homogeneity of G reduce all cube points to this face.
    face = {}
    for partition, value in power.items():
        for exponent in orbit(partition):
            key = exponent[1:]
            face[key] = face.get(key, F(0)) + value
    face = {key: val for key, val in face.items() if val}
    degrees = tuple(max(k[i] for k in face) for i in range(4))
    indices = list(product(*(range(n + 1) for n in degrees)))
    tensor = {k: face.get(k, F(0)) for k in indices}
    for axis, degree in enumerate(degrees):
        nxt = {}
        for key in indices:
            total = F(0)
            for j in range(key[axis] + 1):
                prior = list(key)
                prior[axis] = j
                total += tensor[tuple(prior)] * F(comb(key[axis], j), comb(degree, j))
            nxt[key] = total
        tensor = nxt
    minimum_key = min(tensor, key=tensor.get)
    return face, degrees, tensor, minimum_key


def main():
    began = time.monotonic()
    raw = json.loads(RAW.read_text())
    denominator = list(map(int, raw['denominators']))
    assert denominator[0] == 1
    for n in range(1, 18):
        assert denominator[n] == denominator[n - 1] * (4 * n - 1)
    families = {name: read_family(raw[name], denominator) for name in ['d1', 'd2', 'b']}
    checkpoint_path = ROOT / 'corner_checkpoint.json'
    checkpoint = json.loads(checkpoint_path.read_text())
    names = {'d1': 'E_d_W_power', 'd2': 'E_d_squared_W_power', 'b': 'E_b_corner_W_power'}
    for name, family in families.items():
        for n, poly in enumerate(family):
            assert all(sum(key) == n + (name == 'b') for key in poly)
            assert corner_value(poly) == F(checkpoint['moments'][names[name]][n])
    assert families['d1'][0][(0, 0, 0, 0, 0)] == F(648, 1463)
    assert families['d2'][0][(0, 0, 0, 0, 0)] == F(126404864, 474619145)
    assert families['b'][0] == {(1, 0, 0, 0, 0): F(216, 385)}
    A, B, c = (F(checkpoint['constants'][name]) for name in ['A_lower', 'B_upper', 'c_lower'])
    results = []
    for candidate in checkpoint['candidates']:
        if not candidate['passes_necessary_corner_screen']:
            continue
        eta, a = F(candidate['eta']), F(candidate['a_selected'])
        power = {}
        for n in range(8):
            factor = F(11, 4) * a ** 7 * comb(7, n) * eta ** n
            for key, value in families['d1'][n].items():
                power[key] = power.get(key, F(0)) + factor * A * value
            for key, value in families['d2'][n].items():
                power[key] = power.get(key, F(0)) - factor * B * value
        for n in range(12):
            factor = -F(7, 4) * a ** 11 * comb(11, n) * eta ** n
            for key, value in families['b'][n].items():
                power[key] = power.get(key, F(0)) + factor * value
        power = {key: value for key, value in power.items() if value}
        assert corner_value(power) == F(candidate['corner_polynomial_value'])
        face, degrees, tensor, minimum_key = bernstein_face(power)
        minimum = tensor[minimum_key]
        lower = c * minimum
        simple = F(lower.numerator * 10 ** 8 // lower.denominator, 10 ** 8)
        # Every recorded coefficient is an exact rational global lower enclosure.
        assert all(value >= minimum for value in tensor.values())
        record = {'eta': str(eta), 'a': str(a), 'symmetric_power_coefficients': {','.join(map(str, k)): str(v) for k, v in sorted(power.items())},
                  'face_Bernstein_degrees': degrees, 'face_Bernstein_coefficient_count': len(tensor),
                  'minimum_Bernstein_index': minimum_key, 'minimum_Bernstein_value': str(minimum),
                  'global_g_exact_rational': str(lower), 'global_g_simple_rational': str(simple),
                  'strictly_improves_1_over_27000': simple > F(1, 27000),
                  'positive_global_certificate': minimum > 0,
                  'global_g_float_display_only': float(lower),
                  'Bernstein_coefficients_sha256': hashlib.sha256(json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True).encode()).hexdigest()}
        results.append(record)
    moment_summary = {'families': {name: [{','.join(map(str,k)): str(v) for k,v in sorted(poly.items())} for poly in family] for name,family in families.items()},
                      'crosscheck': 'Every corner value matches independent three-group Cauchy-Binet/Stieltjes checkpoint exactly.',
                      'canonical_graph_moments': raw['canonical_moment_count'], 'raw_graphs_canonicalized': raw['canonical_cache_count']}
    (ROOT / 'space_moments.json').write_text(json.dumps(moment_summary, indent=2) + '\n')
    out = {'status': 'exact_full_face_Bernstein_certificates', 'new_C5_upper_bound': False, 'c_lower': str(c),
           'candidates': results,
           'raw_graph_moment_input_sha256': hashlib.sha256(RAW.read_bytes()).hexdigest(),
           'space_moments_sha256': hashlib.sha256((ROOT / 'space_moments.json').read_bytes()).hexdigest(),
           'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT / 'graph_moments.cpp', ROOT / 'corner_checkpoint.py', checkpoint_path, ROOT / 'registration.md', Path(__file__)]},
           'runtime': {'python': sys.executable, 'seconds': time.monotonic() - began}}
    (ROOT / 'verification.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'candidates': [{key: value for key,value in row.items() if key != 'symmetric_power_coefficients'} for row in results], 'runtime': out['runtime']}, indent=2))


if __name__ == '__main__':
    main()
