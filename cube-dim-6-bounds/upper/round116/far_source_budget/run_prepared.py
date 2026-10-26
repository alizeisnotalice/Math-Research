#!/usr/bin/env python3
"""Two new fixed exact116 scalar certificates. Never runs an upstream program."""
from __future__ import annotations
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction as F

HERE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def atomic(name, data):
    p = HERE / name
    q = p.with_name(p.name + '.own-tmp')
    with q.open('w', encoding='utf8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(q, p)


def req(test, label):
    if not test:
        raise ArithmeticError(label)


def dyad(m, e):
    return F(int(m) * (2 ** max(e, 0)), 2 ** max(-e, 0))


def endpoint(p, side):
    return dyad(p[side + '_mantissa'], p[side + '_binary_exponent'])


def power_enclose(x, power, bits=140):
    """Exact dyadic enclosure of x**power; fixed integer bisection, no floats."""
    req(x >= 0 and power > 0, 'positive power input')
    power = F(power)
    y = x ** power.numerator
    den = 1 << bits
    lo, hi = 0, den
    while F(hi, den) ** power.denominator < y:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if F(mid, den) ** power.denominator <= y:
            lo = mid
        else:
            hi = mid
    l, u = F(lo, den), F(hi, den)
    req(l ** power.denominator <= y <= u ** power.denominator,
        'exact root enclosure')
    return l, u


def interval(values):
    return [str(x) for x in values]


def show_interval(values):
    # Displays are never reused as certificate inputs.
    return [format(float(x), '.17g') for x in values]


def cfg_window(inp, pi):
    piL, piU = pi
    a, b, rho, m = map(F, ['1/64', '1/32', '1/64', '7/1760'])
    alpha = F(1, 32)
    req(piL > 3 and piU < F(22, 7), 'saved pi bounds')
    req(alpha * alpha > rho * rho / (60 * 60 * 6), 'alpha > beta')
    B2pi6 = 64 * sum(F(row['ordered_count']) * F(row['R_exact'])
                    for row in inp['R2'])
    req(sum(row['ordered_count'] for row in inp['R2']) == 992,
        'complete ordered pair multiplicity')
    CfL = 32 * F(inp['R1']) / piU**6
    B2L, B3L = B2pi6 / piU**6, F(inp['B3_pi6_exact']) / piU**6
    P6L = F(216, 2**30) / piU**6
    ds = [P6L / b**6,
          CfL / (b**4 * (1 + 1024*b*b)),
          B2L / (b*b * (1 + 1024*b*b) * (rho*rho + 992*b*b)),
          B3L / ((1 + 1024*b*b) * (rho*rho + 992*b*b)
                  * (rho*rho + 960*b*b))]
    for i, d in enumerate(ds):
        req(d > F(1, 5), 'D' + str(i) + '(b)>1/5')
    pmax_lower = F(23, 24)**31 - F(32, 180**30)
    req(pmax_lower > F(15, 64), 'unique-max source event probability')
    source_factor_lower = a / (F(22, 7) * (F(1, 4) + b*b))
    req(source_factor_lower * F(15, 64) > m, 'Mmax window lower')
    req(F(1, 80) / piU > m, 'positive face L window lower')
    volume_lower = 1 - 150 * rho**4
    req(volume_lower > F(1, 2) and F(1, 10) > m, 'true far volume')
    W5 = (b-a)**6 * m**5 / 5
    RL, RU = F(inp['far_allowance']['lower']), F(inp['far_allowance']['upper'])
    req(0 < RL <= RU, '111 complete raw allowance endpoints')
    req(W5 > RU**5, 'strict failure against even upper remaining allocation')
    W = power_enclose(W5, F(1, 5))
    ratio = power_enclose(W5 / RU**5, F(1, 5))
    return {
        'id': 'true_far_source_unique_max_window',
        'status': 'STRICT_OBSTRUCTION_TO_EXPLICIT_NEW_ENVELOPE_FAMILY',
        'fixed_window': [str(a), str(b)], 'rho': str(rho),
        'D_at_b_certified_lower': interval(ds),
        'D_at_b_display': show_interval(ds),
        'B2_pi6_exact_from_saved_R2': str(B2pi6),
        'source_probability_lower': str(pmax_lower),
        'source_probability_minus_15_over64': str(pmax_lower - F(15, 64)),
        'source_mass_window_lower': str(source_factor_lower * F(15, 64)),
        'far_volume_lower': str(volume_lower), 'all_mass_entries_lower': str(m),
        'I_envelope_lower_to_fifth': str(W5),
        'I_envelope_lower_interval': interval(W),
        'I_envelope_lower_display': show_interval(W),
        'obstruction_to_RU_ratio_interval': interval(ratio),
        'obstruction_to_RU_ratio_display': show_interval(ratio),
        'residual_to_RL_fifth': str(W5-RL**5),
        'residual_to_RU_fifth': str(W5-RU**5),
        'actual_I_far_lower_bound': False,
        'full_far_paid': False,
        'stop': 'this family stops; no time refinement or changed source parameters',
    }


def cfg_tail(inp, pi):
    piL, piU = pi
    T = F(1, 8)
    pL, pU = F(216, 2**30) / piU**6, F(216, 2**30) / piL**6
    tailL, tailU = pL/(5*T**5), pU/(5*T**5)
    RL = F(inp['far_allowance']['lower'])
    residual = RL**5 / 20**6 - tailU**6
    req(residual > 0, 'entire T>=1/8 tail uses less than 1/20 far norm')
    remaining = power_enclose(RL, F(5, 6))
    short = (remaining[0] - tailU, remaining[1] - tailL)
    req(short[0] > 0, 'strict positive short-time norm allocation')
    tailI = (power_enclose(tailL, F(6, 5))[0],
             power_enclose(tailU, F(6, 5))[1])
    # Existing local mathematical tool: independent Arb enclosure for the NEW
    # rank6 tail only. No upstream peak, body kernel or grid is recomputed.
    from flint import arb, ctx
    ctx.prec = 384
    z = arb(216) / (arb(2)**30 * arb.pi()**6 * arb(5) * (arb(1)/8)**5)
    tail_outL = F(math.floor(tailL * 10**30), 10**30)
    tail_outU = F(math.ceil(tailU * 10**30), 10**30)
    # The positive rational endpoints alone certify the bound. This check
    # verifies the independent Arb evaluation lies inside a coarser enclosure.
    loA = arb(tail_outL.numerator) / tail_outL.denominator
    hiA = arb(tail_outU.numerator) / tail_outU.denominator
    req(bool(z > loA) and bool(z < hiA), 'new tail Arb crosscheck')
    return {
        'id': 'accepted108_full6D_far_large_time_tail', 'status': 'STRICT_PAID_TAIL',
        'T': str(T), 'P6_interval': interval((pL, pU)),
        'old095_peak': '1/(160pi^6)', 'new_peak': '216/(2^30pi^6)',
        'new_over_old_exact': str(F(160*216, 2**30)),
        'tail_norm_interval': interval((tailL, tailU)),
        'tail_norm_display': show_interval((tailL, tailU)),
        'tail_I_if_alone_interval': interval(tailI),
        'tail_I_if_alone_display': show_interval(tailI),
        'tail_below_one_twentieth_remaining_norm_residual': str(residual),
        'short_time_norm_allocation_interval': interval(short),
        'short_time_norm_allocation_display': show_interval(short),
        'arb384_new_tail_crosscheck': str(z),
        'norm_addition_required_before_q_power': True,
        'full_far_paid': False, 'global_g0_proved': False, 'new_C6_bound': False,
    }


def main():
    start, cpu = time.monotonic(), time.process_time()
    reg = json.loads((HERE/'registration.json').read_text())
    release = json.loads((HERE/'parent_release.json').read_text())
    req(release.get('authorized') is True and str(release.get('round')) == '116',
        'explicit main release')
    req(release['script_sha256'] == sha(__file__), 'released script hash')
    req(release['registration_sha256'] == sha(HERE/'registration.json'),
        'released registration hash')
    req(sha(HERE/'upstream.json') == reg['upstream_sha256'], 'upstream pin')
    req(sha(HERE/'proof.md') == reg['proof_sha256'], 'proof pin')
    req(sha(HERE/'source_manifest.json') == reg['source_manifest_sha256'],
        'source manifest pin')
    sources = json.loads((HERE/'source_manifest.json').read_text())
    for row in sources['sources']:
        req(sha(row['path']) == row['sha256'], 'source pin '+row['role'])
    inp = json.loads((HERE/'upstream.json').read_text())
    flag = HERE/'started.json'
    with flag.open('x') as f:
        json.dump({'pid': os.getpid(), 'utc': utc(), 'script_sha256': sha(__file__),
                   'release_sha256': sha(HERE/'parent_release.json')}, f, indent=2)
        f.write('\n')
    status, done, failure = 'FAILED', [], None
    try:
        pi = (endpoint(inp['pi_saved036'], 'lower'),
              endpoint(inp['pi_saved036'], 'upper'))
        done.append(cfg_window(inp, pi))
        atomic('checkpoint.json', {'completed_configurations': done, 'count': len(done)})
        done.append(cfg_tail(inp, pi))
        atomic('checkpoint.json', {'completed_configurations': done, 'count': len(done)})
        status = 'COMPLETED_2_FIXED_CONFIGURATIONS'
        atomic('results.json', {'status': status, 'configurations': done,
                               'release_sha256': sha(HERE/'parent_release.json'),
                               'new_C6_bound': False, 'full_far_paid': False})
        print(json.dumps({'status': status, 'pid': os.getpid(),
                          'window_I_lower_display': done[0]['I_envelope_lower_display'],
                          'tail_norm_display': done[1]['tail_norm_display']}), flush=True)
    except BaseException as e:
        failure = type(e).__name__ + ': ' + str(e)
        raise
    finally:
        ru = resource.getrusage(resource.RUSAGE_SELF)
        atomic('run_receipt.json', {
            'status': status, 'completed': len(done), 'failure': failure,
            'pid': os.getpid(), 'utc': utc(), 'CPU_seconds': time.process_time()-cpu,
            'wall_seconds_sample_before_final_receipt': time.monotonic()-start,
            'timing_scope': 'start before release checks through completed results; sample excludes this receipt write and interpreter exit',
            'max_rss_bytes_macos': ru.ru_maxrss, 'children': 0,
            'runtime_limit': None, 'CPU_limit': None,
            'source_files_read': len(sources['sources']),
            'old_native_or_peak_or_grid_programs_executed': 0,
        })


if __name__ == '__main__':
    main()
