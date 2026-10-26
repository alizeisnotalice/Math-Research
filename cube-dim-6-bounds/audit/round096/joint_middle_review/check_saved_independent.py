#!/usr/bin/env python3
"""Verify one saved certificate group, not the author scientific batch.

Only Fraction identities and saved root inequalities. No old F/L values,
time grid, root search, semigroup, kernel or author-program import.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import hashlib
import json
import resource
import signal
import sys
import time

H = Path(__file__).resolve().parent
WALL, CPU = time.monotonic(), time.process_time()
resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('wall limit')))
signal.alarm(165)
checks = 0
counts = {'saved_certificate_groups': 1, 'source_hashes': 0, 'saved_root_certificates': 0,
          'profiles': 0, 'candidate_budget_comparisons': 0,
          'new_parameter_configurations': 0, 'new_root_searches': 0,
          'author_program_executions': 0, 'old_F_L_replays': 0,
          'old_86_cells_replayed': 0, 'kernel_calls': 0}
out = {'status': 'RUNNING'}
def ck(ok, msg):
    global checks
    checks += 1
    if not ok:
        raise AssertionError(msg)
    if time.process_time()-CPU >= 18:
        raise TimeoutError('CPU prelimit')
def load(p):
    return json.loads(Path(p).read_text())
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

try:
    gate = load(H/'source_gate.json')
    for p, expected in gate['source_pins'].items():
        ck(sha(p) == expected, 'source pin '+p)
        counts['source_hashes'] += 1
    author = Path(next(p for p in gate['source_pins'] if p.endswith('/round096/joint_middle_budget/manifest.json'))).parent
    r, reg, manifest, par = [load(author/n) for n in ('results.json', 'registration.json', 'manifest.json', 'parameters.json')]
    ck(sha(author/'manifest.json') == gate['author_manifest_sha256'], 'manifest seal')
    ck(sha(author/'results.json') == gate['author_results_sha256'], 'results seal')
    for name, ent in manifest['files'].items():
        ck((author/name).stat().st_size == ent['bytes'] and sha(author/name) == ent['sha256'], 'manifest artifact '+name)
    source = list(reg['source_sha256'])
    prior = load(next(p for p in source if p.endswith('/round092/middle_envelope_gap/lower_results.json')))
    ledger = load(next(p for p in source if p.endswith('/round062/primary_joint_account/actual_run/results.json')))
    ck(prior['status'] == 'STRICT_FIXED_ENVELOPE_NONCOMPETITIVE', '092 source status')
    q, b, g0 = Q(par['q']), Q(par['b']), Q(par['g0'])
    ck((q, b, g0) == (Q(6,5), Q(1,64), Q(1,65000)), 'fixed model parameters')
    ck(q/4 == Q(r['full_time_integrability_exponent']) == Q(3,10) < 1, 'rank integrability exponent')
    cJ, cW, piU = map(Q, (prior['cJ_lower'], prior['cW_lower'], prior['pi_upper']))
    ck(min(cJ, cW, piU) > 0, 'positive saved constants')
    # Check both whole-mass branches from the inherited rational pi cap.
    ck(cJ <= Q(16384,377)/piU**2, 'uniform marginal coefficient direction')
    ck(cW <= 5*cJ and cW <= Q(256,5)/piU, 'both whole mass branches')
    common = prior['common_interval_checks']
    ck(sorted(x['k'] for x in common) == list(range(1,6)), 'all five projection types')
    for x in common:
        ck(Q(x['z_lower_upper_endpoint']) <= Q(1,4) and Q(x['z_upper_lower_endpoint']) >= Q(3,4), 'common spatial interval containment')
    ck(Q(15,64) <= Q(1,4) < Q(3,4) <= Q(49,64), 'whole projection interval containment')
    coeff = {}
    for k in range(1,6):
        x = prior['f_coefficients'][str(k)]
        f, v, n = Q(x['uniform_F_lower']), Q(x['V_lower']), x['n']
        ck(f > 0 and v > 0 and n == comb(5,k-1), 'source coefficient positivity/multiplicity')
        coeff['k'+str(k)] = (f,v,cJ,n)
    coeff['whole_union'] = (max(x[0] for x in coeff.values()), sum(x[1]*x[3] for x in coeff.values()), cW, 1)
    ck(sum(comb(5,k-1) for k in range(1,6)) == 31, 'tube count')
    supplied = {x['name']:(Q(x['f']),Q(x['v']),Q(x['c']),x['multiplicity']) for x in r['source_coefficients']}
    ck(supplied == coeff and set(r['profiles']) == set(coeff), 'exact source coefficient handoff')
    # Every recorded root is certified by raising rational endpoints to integer powers.
    roots = {}
    for x in r['root_certificates']:
        label, d = x['label'], x['degree']
        rad, lo, hi = map(Q, (x['radicand'],x['lower'],x['upper']))
        ck(label not in roots and isinstance(d,int) and d in (2,5,6), 'root label/degree')
        ck(0 <= lo <= hi and lo**d <= rad <= hi**d, 'root outward inequality '+label)
        roots[label] = (d,rad,lo,hi)
        counts['saved_root_certificates'] += 1
    used = set()
    def rt(label, degree, radicand):
        d, rad, lo, hi = roots[label]
        ck(d == degree and rad == radicand, 'root formula binding '+label)
        used.add(label)
        return lo, hi
    bpow = rt('b-to-q-plus2', 5, b**16)
    witness = {}
    profile_reviews = {}
    for name, (f,v,c,n) in coeff.items():
        p = r['profiles'][name]
        ck((Q(p['f']),Q(p['input_v']),Q(p['c']),p['multiplicity'],Q(p['b'])) == (f,v,c,n,b), 'profile coefficient handoff '+name)
        tau = rt(name+'-space-cutoff', 2, v*f/c)
        m = Q(p['rational_cutoff_m']); z = b-m
        ck(m == min(b,tau[0]) and 0 < m <= b and m*m <= v*f/c, 'safe support cutoff '+name)
        actual_v = c*m*m/f
        ck(actual_v == Q(p['actual_used_volume']) and 0 < actual_v <= v, 'actual support volume '+name)
        # These identities certify the two entire time regimes by monotonicity:
        # 0<t<=m: mass=c*t^2; m<t<=b: mass=c*m^2<=c*t^2.
        # Outside (0,b] mass=0; peak=f<=F and support<=actual_v<=v<=V.
        ck(f*actual_v == c*m*m and c > 0, 'all-time mass-regime identity '+name)
        fp = rt(name+'-f-to-q-minus1', 5, f)
        zp1 = rt(name+'-z-to-q-plus1', 5, z**11)
        zp2 = rt(name+'-z-to-q-plus2', 5, z**16)
        # Negative term takes the upper root for a lower integral bound.
        low_bracket = bpow[0]/((q+1)*(q+2))-b*zp1[1]/(q+1)+zp2[0]/(q+2)
        high_bracket = bpow[1]/((q+1)*(q+2))-b*zp1[0]/(q+1)+zp2[1]/(q+2)
        ck(0 < low_bracket <= high_bracket, 'positive outward primitive '+name)
        lo, hi = 2*c*fp[0]*low_bracket, 2*c*fp[1]*high_bracket
        ck((lo,hi) == (Q(p['q_integral']['lower']),Q(p['q_integral']['upper'])) and 0 < lo <= hi, 'saved integral endpoints '+name)
        witness[name] = (lo,hi)
        profile_reviews[name] = {'actual_volume':str(actual_v),'cutoff_m':str(m),
            'window':'(0,1/64]; identically zero outside',
            'mass_regimes':'c*t^2 for t<=m; c*m^2 for m<t<=b',
            'height_mass_volume':'PASS relative to inherited F/L/V upper constraints',
            'I_lower':str(lo),'I_upper':str(hi)}
        counts['profiles'] += 1
    tube = tuple(sum(comb(5,k-1)*witness['k'+str(k)][j] for k in range(1,6)) for j in (0,1))
    whole = witness['whole_union']
    ck(tube == (Q(r['tube_witness_I_sum']['lower']),Q(r['tube_witness_I_sum']['upper'])), 'weighted tube sum')
    ck(whole == (Q(r['whole_witness_I']['lower']),Q(r['whole_witness_I']['upper'])), 'whole sum')
    new = min(tube[0],whole[0])
    ck(new == Q(r['new_min_relaxation_I_lower']), 'outer min of both lower bounds')
    allowI = Q(ledger['sufficient_remaining_D_F_I_allowance'])
    allowB = Q(ledger['sufficient_remaining_D_F_B_allowance'])
    paid = Q(ledger['paid_C_E_H_I_upper'])
    ck((allowI,allowB,paid) == (Q(r['saved_D_F_I_allowance']),Q(r['saved_D_F_B_allowance']),Q(r['paid_CEH_I_upper'])), '062 exact handoff')
    budgets = {}
    for name, I, rootname, field in [('tube',tube[0],'tube-relaxation-B-lower','tube_relaxation_B_lower'),
                                      ('whole',whole[0],'whole-relaxation-B-lower','whole_relaxation_B_lower'),
                                      ('min',new,'new-min-relaxation-B-lower','new_min_relaxation_B_lower')]:
        lower = rt(rootname,6,I**5)[0]
        ck(lower == Q(r[field]), 'B endpoint '+name)
        ck(I > allowI and lower > allowB, 'both I and B exceed allowance '+name)
        budgets[name] = {'I_lower':str(I),'B_lower':str(lower),'I_gap':str(I-allowI),'B_gap':str(lower-allowB)}
        if name != 'min':
            counts['candidate_budget_comparisons'] += 1
    ck(Q(r['I_allowance_gap']) == new-allowI and Q(r['B_allowance_gap']) == Q(r['new_min_relaxation_B_lower'])-allowB, 'saved gap directions')
    target_unit = rt('conditional-I-target-unit',5,(g0/6)**6)
    target = tuple(64*x for x in target_unit)
    ck(target == (Q(r['conditional_I_target']['lower']),Q(r['conditional_I_target']['upper'])), 'I target outward factor')
    Iplan = paid+new
    ck(Iplan == Q(r['conditional_plan_I_lower']), 'conservative I plan lower')
    factor = rt('conditional-U-factor',5,(6/g0)**6)
    Ulo, Sup = factor[0]*Iplan,target[1]-Iplan
    ck(Ulo == Q(r['conditional_plan_U_lower']) and Sup == Q(r['conditional_plan_S_upper']), 'U/S formulas and outward signs')
    ck(Sup < 0 and Ulo > 64, 'conditional plan strictly fails')
    ck(Q(r['conditional_g0']) == g0 and r['global_g0_proven'] is False, 'conditional g accounting')
    ck(r['actual_density_witness'] is False and r['actual_model_cost_lower_bound'] is False and r['new_C6_bound'] is False, 'exclusion scope flags')
    ck(r['competitive_D_budget_paid'] is False and r['far_budget_paid'] is False and r['full_time_D_requires_no_separate_short_addition'] is True, 'full-time/far accounting')
    ck(len(used) == len(roots) == 30 and set(roots) == used, 'all saved root bindings complete')
    for p, expected in gate['source_pins'].items():
        ck(sha(p) == expected, 'post-run source immutability')
    out.update(status='PASS_SAVED_CERTIFICATE_AND_THEORETICAL_SCOPE',profiles=profile_reviews,
               budgets=budgets,conditional_plan={'g0':str(g0),'globally_proved':False,
               'I_lower':str(Iplan),'U_lower':str(Ulo),'S_upper':str(Sup)},
               inherited_scope='074 F/L/V and 092 uniform lower coefficient proofs reviewed and hash-bound; no old function recalculation',
               verdict='Only both separately solved sharp mass/peak/volume relaxations and their min are noncompetitive under the fixed 062 conservative allocation; no actual model or simultaneous-overlap optimum exclusion.')
except Exception as exc:
    out.update(status='FAIL_FIRST_FAILURE',error=repr(exc))
finally:
    out.update(checks=checks,counts=counts,CPU_seconds_before_save=time.process_time()-CPU,
               wall_seconds_before_save=time.monotonic()-WALL,
               self_sha256=sha(Path(__file__)),source_gate_sha256=sha(H/'source_gate.json'))
    with (H/'review_results.json').open('x') as f:
        json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:out[k] for k in ('status','checks','counts','CPU_seconds_before_save','wall_seconds_before_save')}))
    signal.alarm(0)
if out['status'] != 'PASS_SAVED_CERTIFICATE_AND_THEORETICAL_SCOPE':
    sys.exit(1)
