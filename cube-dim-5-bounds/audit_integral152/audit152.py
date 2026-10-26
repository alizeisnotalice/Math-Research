"""Prepared finite73-factor/64-root independent checker. Not yet executed.

Time weights use integrated barycentric polynomials; packet sums use common
dyadic integer denominators. No producer module, Fraction, or power root.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import time

from rational152 import Q, ceil_div, dyadic_bits, on_grid

HERE = Path(__file__).resolve().parent
D = 1 << 192
N = 1 << 160


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    path = Path(path)
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def write(path, obj):
    with Path(path).open('x') as handle:
        json.dump(obj, handle, sort_keys=True, separators=(',', ':'))
        handle.write('\n')


def audit(plan, binding, output):
    start = time.monotonic()
    for path, h in plan['input_locks'].items():
        assert sha(path) == h, path
    for path, h in binding['candidate_locks'].items():
        assert plan['input_locks'].get(path, h) == h, path
        assert sha(path) == h, path
    paths = plan['source_paths']
    source_check = read(paths['source_check'])
    assert source_check['status'] == 'selected64_original_active_chord_root_geometry_packet_identity_bound'
    assert source_check['selected'] == 64 and source_check['all_unexpanded_roots']
    assert source_check['all_strict112_chord_arm_active'] and source_check['all_U_at_least_1over65536']
    assert source_check['new_factor_mean_H_R_or_old_power_calls'] == 0
    pre, selected = read(paths['preflight']), read(paths['selected64'])
    inherited = read(paths['audit121_results'])
    assert inherited['status'] == inherited['final']['status'] == 'pass'
    assert inherited['final']['final_budget']['UH'] == pre['saved_complete_UH']
    geos, old, catalog = read(paths['geometries']), read(paths['old_factors']), read(paths['old_ratio'])
    author = Path(binding['output_directory'])
    af, an, ar = read(author/'factors.json'), read(author/'nodes.json'), read(author/'results.json')
    state, reg = read(author/'state.json'), read(author/'registration.json')
    for name in ['factors.json', 'nodes.json', 'results.json', 'state.json', 'registration.json']:
        assert str(author/name) in binding['candidate_locks']
    assert str(binding['producer_final_manifest_path']) in binding['candidate_locks']
    identity = {'plan_sha256': plan['producer_plan_sha256'],
                'selected64_sha256': plan['producer_selected64_sha256'],
                'consumer_sha256': plan['producer_consumer_sha256']}
    assert state['status'] == 'complete' and state['identity'] == identity
    assert state['invocations'] == binding['actual_producer_invocations']
    assert reg['pid'] == state['invocations'][0]
    assert ar['pid'] == state['invocations'][-1] == binding['actual_producer_pid']
    assert all(reg[k] == v for k, v in identity.items())
    assert state['factors'] == af and state['nodes'] == an
    assert len(af) == 73 and len(an) == len(selected) == 64
    assert ar['status'] == 'finite64root_candidate_complete_pending_independent_primary'
    assert ar['factors'] == 73 and ar['roots'] == 64 and ar['cached_band_products'] == 4672
    assert ar['registration_sha256'] == sha(author/'registration.json')
    assert ar['new_H_R_gradient_log_exp_powerroot_calls'] == ar['old_baseline_recomputations'] == 0
    assert ar['no_powerroot_replay'] and ar['all_sources_unchanged'] and not ar['global_Cn_claim']
    assert reg['roots'] == 64 and reg['factors'] == 73 and reg['cached_band_products'] == 4672
    times = pre['time_endpoints']
    assert len(times) == 73 and times[0] == '1/512'
    h, mu = Q(geos['g002']['profile']['h']), Q(geos['g002']['profile']['mu_infinity'])
    factors = []
    with (output/'factor_checks.jsonl').open('x') as log:
        for j, at in enumerate(times):
            bt = times[j+1] if j+1 < len(times) else None
            a, b = Q(1)/Q(at), Q(0) if bt is None else Q(1)/Q(bt)
            assert Q(0) <= b < a
            # Integrals of z^3 and z^4 over [b,a], not the producer theta formula.
            m0, m1 = (a**4-b**4)/4, (a**5-b**5)/5
            left = (m1-b*m0)/((a-b)*m0)
            right = 1-left
            assert Q(0) < left < Q(1) and left+right == Q(1)
            if bt is None:
                assert left == Q(4, 5)
            ra = Q(catalog[str(h*a)]['upper'])
            ta = 1+mu*a+Q(25, 24)*Q(old['g002'][at]['J_upper'])*a**2
            rb = Q(1) if bt is None else Q(catalog[str(h*b)]['upper'])
            tb = Q(1) if bt is None else 1+mu*b+Q(25, 24)*Q(old['g002'][bt]['J_upper'])*b**2
            rb, tb = min(rb, ra), min(tb, ta)
            rchord, tchord = left*ra+right*rb, left*ta+right*tb
            oldc = Q(old['g002'][at]['factor_upper'])
            old_units = on_grid(oldc, 192)
            exact = min(oldc, rchord, tchord)
            units = ceil_div(exact.n*D, exact.d)
            assert D <= units <= old_units
            row = af[j]
            assert row['a'] == at and row['b'] == bt and Q(row['theta']) == left
            assert Q(row['range_time_upper']) == rchord and Q(row['taylor_time_upper']) == tchord
            assert Q(row['old_upper']) == oldc and Q(row['new_upper']) == Q(units, D)
            factors.append(units)
            log.write(json.dumps({'band': j, 'a': at, 'b': bt, 'left_weight': str(left),
                                 'right_weight': str(right), 'new_factor_units192': units,
                                 'exact_factor_before_ceil': str(exact)}, sort_keys=True)+'\n')
    before = on_grid(pre['saved_complete_UH'], 160)
    gain_sum = 0
    with (output/'root_checks.jsonl').open('x') as log:
        for i, item in enumerate(selected):
            r, c = item['source_node'], item['source112_certificate']
            packet = read(r['packet_path'])
            assert sha(r['packet_path']) == r['packet_sha256']
            assert packet['key'] == item['key'] and packet['point'] == r['point']
            endpoints = packet['endpoints']
            assert [e['t'] for e in endpoints] == times
            vals = [(Q(e['H']['lower']), Q(e['H']['upper'])) for e in endpoints]
            bits = max(dyadic_bits(q) for pair in vals for q in pair)
            ints = [(q.n << (bits-dyadic_bits(q)), u.n << (bits-dyadic_bits(u))) for q, u in vals]
            assert all(0 < lo <= up for lo, up in ints)
            total = 0
            for j, (lo, up) in enumerate(ints):
                nl, nu = ints[j+1] if j+1 < len(ints) else (0, 0)
                assert 0 < lo-nu <= up-nl
                total += (up-nl)*factors[j]
            mean_i = min(c['new_M_i'], ceil_div(total, 1 << (bits+192-160)))
            assert r['old_L_i'] <= mean_i <= c['new_M_i'] <= r['old_U_i']
            gain = (r['weight']*(c['new_M_i']-mean_i)) >> (4*r['level']+4)
            new_upper = item['saved_selfU']-gain
            assert item['saved_selfL'] <= new_upper <= item['saved_selfU']
            gain_sum += gain
            row = an[i]
            assert row['rank'] == i and row['key'] == item['key']
            assert row['old_mean_i'] == c['new_M_i'] and row['new_mean_i'] == mean_i
            assert row['safe_secant_slope_lower'] == '1/16'
            assert row['old_selfU'] == item['saved_selfU'] and row['new_selfU'] == new_upper
            assert row['gain_units'] == gain and row['prefix_complete_UH_units'] == before-gain_sum
            assert row['packet_sha256'] == r['packet_sha256']
            log.write(json.dumps({'key': item['key'], 'common_H_denominator_bits': bits,
                                 'new_mean_units160': mean_i, 'gain_units160': gain,
                                 'prefix_complete_UH_units160': before-gain_sum}, sort_keys=True)+'\n')
    assert state['gain_units'] == gain_sum
    assert Q(ar['old_UH']) == Q(before, N) and Q(ar['new_UH']) == Q(before-gain_sum, N)
    assert Q(ar['exact_UH_gain']) == Q(gain_sum, N)
    assert ar['strict_global_UH_improvement'] == (gain_sum > 0)
    assert ar['mean_improved_roots'] == sum(r['new_mean_i'] < r['old_mean_i'] for r in an)
    for path, h in {**plan['input_locks'], **binding['candidate_locks']}.items():
        assert sha(path) == h, path
    result = {'status': 'finite73factor64root_arithmetic_independently_qualified', 'pid': os.getpid(),
              'factors': 73, 'roots': 64, 'cached_band_products': 4672,
              'old_UH': str(Q(before, N)), 'new_UH': str(Q(before-gain_sum, N)),
              'exact_gain': str(Q(gain_sum, N)), 'strict_UH_improvement': gain_sum > 0,
              'all_sources_before_after_unchanged': True, 'old_power_or_H_R_calls': 0,
              'global_Cn_claim': False, 'seconds': time.monotonic()-start}
    write(output/'result.json', result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reviewed-plan-sha', required=True)
    parser.add_argument('--binding', type=Path, required=True)
    parser.add_argument('--reviewed-binding-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert sha(HERE/'audit_plan.json') == args.reviewed_plan_sha
    assert sha(args.binding) == args.reviewed_binding_sha
    plan, binding = read(HERE/'audit_plan.json'), read(args.binding)
    assert binding['execution_ready'] is True and binding['stage_id'] == 'fixed150_selected64'
    assert binding['audit_plan_sha256'] == args.reviewed_plan_sha
    for name, h in plan['code_locks'].items():
        assert sha(HERE/name) == h
    assert args.output.is_absolute() and args.output.is_relative_to(HERE) and not args.output.exists()
    args.output.mkdir()
    write(args.output/'registration.json', {
        'task': 'audit_integral152_unique_finite_batch', 'pid': os.getpid(),
        'utc': datetime.now(timezone.utc).isoformat(), 'plan_sha256': args.reviewed_plan_sha,
        'binding_sha256': args.reviewed_binding_sha, 'factors': 73, 'roots': 64,
        'cached_band_products': 4672, 'old_power_or_H_R_calls': 0})
    audit(plan, binding, args.output)


if __name__ == '__main__':
    main()
