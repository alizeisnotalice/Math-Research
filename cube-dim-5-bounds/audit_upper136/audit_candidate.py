#!/usr/bin/env python3
"""Prepared136 audit: direct quartic Arb plus exact rational producer-cap check.

One complete candidate only, after main review of plan, inputs, code hashes.
No old scalar, H, R, kernel, gradient, log, exp or Gamma computation.
"""
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import os
import time
from flint import arb, ctx, __version__ as flint_version

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write('\n')

def read_locked(row):
    assert sha(row['path']) == row['sha256'], row['path']
    return json.loads(Path(row['path']).read_bytes())

def selected(data, keys):
    for key in keys:
        data = data[key]
    return data

def ball(q):
    q = Q(q)
    return arb(q.numerator) / arb(q.denominator)

def endpoints(value):
    assert value.is_finite()
    lo, hi = Q(str(value.lower().fmpq())), Q(str(value.upper().fmpq()))
    assert lo <= hi
    return lo, hi

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reviewed-plan-sha', required=True)
    parser.add_argument('--inputs', required=True, type=Path)
    parser.add_argument('--reviewed-inputs-sha', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    planpath = HERE / 'proposed_plan.json'
    assert sha(planpath) == args.reviewed_plan_sha
    assert sha(args.inputs) == args.reviewed_inputs_sha
    plan = json.loads(planpath.read_bytes())
    assert sha(__file__) == plan['audit_code_sha256']
    assert plan['maximum_physical_cases'] == 1
    inputs = json.loads(args.inputs.read_bytes())
    sources = {name: read_locked(row) for name, row in inputs['sources'].items()}
    args.output.mkdir(exist_ok=False)
    write(args.output/'registration.json', {
        'task': 'audit_upper136_one_complete_join_candidate', 'pid': os.getpid(),
        'utc': datetime.now(timezone.utc).isoformat(),
        'code_sha256': sha(__file__), 'plan_sha256': args.reviewed_plan_sha,
        'inputs_sha256': args.reviewed_inputs_sha,
        'sources': inputs['sources'], 'M_pointer': inputs['M_pointer'],
        'precision_bits': plan['precision_bits'], 'maximum_physical_cases': 1,
        'workers': 1, 'wall_deadline': None,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls': 0,
        'old_baseline_recomputations': 0,
    })
    begin = time.monotonic()
    eps, eta, g = map(Q, (inputs['epsilon'], inputs['eta'], inputs['g']))
    assert eps == Q(1,512) and eta == Q(1,10) and g == Q(21,200000)
    assert sources['g105']['qualified_g'] == str(g)
    const = sources['stored076']['derived_constants']
    c3l, c3u = Q(const['C3']['lower']), Q(const['C3']['upper'])
    c1l, c1u = Q(const['C1']['lower']), Q(const['C1']['upper'])
    pil = Q(const['pi_outward_384_dyadic']['lower'])
    piu = Q(const['pi_outward_384_dyadic']['upper'])
    assert 0 < c3l < c3u and 0 < c1l < c1u and 3 < pil < piu < Q(22,7)
    condition_margin = Q(20,83)*c3l - c1u*eps**2/(30*eta**2)
    assert condition_margin > 0, '134 sufficient condition not separated'
    Craw = Q(20,83)*c3u + c1u*eps**2
    assert sources['condition134']['epsilon'] == str(eps)
    assert sources['condition134']['eta'] == str(eta)
    Cbar = Q(sources['condition134']['Cbar_upper_256dyadic'])
    assert Cbar >= Craw, '134 dyadic Cbar is below the reconstructed sufficient value'
    assert inputs['M_qualification']['epsilon'] == str(eps)
    assert inputs['M_qualification']['domain'] == 'complete_F'
    assert inputs['M_qualification']['bound'] == 'H_epsilon<=M'
    M = Q(selected(sources['M132'], inputs['M_pointer']))
    IH = Q(sources['integrals014']['complete_H512_integral_upper'])
    candidate = sources['candidate135']
    assert Q(candidate['epsilon']) == eps
    assert Q(candidate['Cbar_upper']) == Cbar
    assert Q(candidate['M_upper']) == M > 0
    assert Q(candidate['IH_upper']) == IH >= 0
    assert [Q(x) for x in candidate['pi_interval']] == [pil,piu]
    jlo, ju = Q(candidate['cap_lower']), Q(candidate['cap_upper'])
    assert 0 < jlo <= ju
    # Cancellation makes the exact quartic ceiling affine in Cbar and M.
    numerator = 625*eps**6*(2*Cbar+3*M*eps**2)
    j4lo = numerator/(9072*piu**4)
    j4hi = numerator/(9072*pil**4)
    assert jlo**4 <= j4lo, 'producer lower endpoint is too high'
    assert ju**4 >= j4hi, 'producer upper endpoint does not bound quartic ceiling'
    iu = IH+ju
    assert iu == Q(candidate['conditional_complete_I_upper'])
    exact_gate = 32**4*g**5 - 3125*iu**4
    assert exact_gate != 0, 'exact32 gate is not separated'
    assert (exact_gate > 0) == candidate['conditional_cost_upper_below_32']
    # Seven exact rational links, all consuming the same complete physical inputs.
    caps = inputs['short_rational_chain']
    cs, ms, ps, hs, rs, Is, Ks = map(Q, (
        caps['Cbar_upper'], caps['M_upper'], caps['pi_lower'], caps['IH_upper'],
        caps['root_argument_fourth_root_upper'], caps['I_upper'], caps['cost_upper']))
    assert [cs,ms,ps,hs,rs,Is,Ks] == [Q(1,300000),Q(2),Q(157,50),
        Q(44372,10**9),Q(33,32),Q(2271,50000000),Q(799,25)]
    root_argument_short = (2*cs/eps**2+3*ms)/7
    increment_short = 5*eps**2*rs/(6*ps)
    short_margins = {
        'Cbar_below_1over300000': cs-Cbar,
        'M_below2': ms-M,
        'pi_lower_above157over50': pil-ps,
        'IH_below44372over1e9': hs-IH,
        'root_argument_below_33over32_fourth_power': rs**4-root_argument_short,
        'complete_I_below2271over50000000': Is-hs-increment_short,
        'cost_cap799over25_fourth_power_gate': Ks**4*g**5-3125*Is**4,
    }
    assert all(x>0 for x in short_margins.values()), 'short rational chain not separated'
    assert Ks < 32
    ctx.prec = plan['precision_bits']
    ctx.threads = 1
    # Different root path: sqrt twice on the entire quartic, not on 135's mean argument.
    ceiling4 = ball(numerator)/(9072*ball(pil)**4)
    ceiling = ceiling4.sqrt().sqrt()
    cl, cu = endpoints(ceiling)
    assert cl > 0 and cu**4 >= j4hi
    assert cl**4 <= j4hi
    arb_iu = IH+cu
    arb_gate = 32**4*g**5 - 3125*arb_iu**4
    assert arb_gate != 0, 'independent rational32 gate not separated'
    cost = (3125*ball(arb_iu)**4/ball(g)**5).sqrt().sqrt()
    costlo, costhi = endpoints(cost)
    assert costlo > 0 and costlo**4 <= 3125*arb_iu**4/g**5 <= costhi**4
    assert (arb_gate > 0 and costhi < 32) or (arb_gate < 0 and costlo > 32)
    for row in inputs['sources'].values():
        assert sha(row['path']) == row['sha256']
    write(args.output/'results.json', {
        'task': 'audit_upper136', 'status': 'one_join_candidate_audited_pending_primary',
        'pid': os.getpid(), 'seconds': time.monotonic()-begin,
        'flint_version': flint_version, 'precision_bits': ctx.prec,
        'epsilon': str(eps), 'eta': str(eta), 'g': str(g),
        'condition134_positive_margin': str(condition_margin),
        'Cbar_raw_reconstructed': str(Craw), 'Cbar_134_dyadic_consumed': str(Cbar),
        'M_consumed': str(M), 'IH_consumed': str(IH),
        'quartic_ceiling_interval': [str(j4lo),str(j4hi)],
        'producer_cap_endpoints_pass_exact_fourth_powers': True,
        'producer_complete_I_upper': str(iu), 'producer_exact32_gate': str(exact_gate),
        'independent_Jensen_ceiling_interval': [str(cl),str(cu)],
        'independent_complete_I_upper': str(arb_iu),
        'independent_exact32_gate': str(arb_gate),
        'independent_cost_interval': [str(costlo),str(costhi)],
        'cost_display_only': [float(costlo),float(costhi)],
        'cost_strictly_below32': arb_gate>0,
        'short_rational_chain_exact_positive_margins': {k:str(v) for k,v in short_margins.items()},
        'short_rational_chain_complete_I_upper': str(Is),
        'short_rational_chain_cost_upper': str(Ks),
        'actual_J_R_or_I_lower_claim': False, 'new_Cn_claim_without_primary_acceptance': False,
        'shared_sources': ['076 constants','057 normalization proof','132 M endpoint','111/112/121 UH','g105 cover'],
        'all_source_hashes_unchanged': True,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls': 0,
        'old_baseline_recomputations': 0,
        'registration_sha256': sha(args.output/'registration.json'),
    })
    print('PASS: one candidate; exact producer quartic check and independent Arb cost path')

if __name__ == '__main__':
    main()
