#!/usr/bin/env python3
"""Exact rational finite guard for a fixed-delay continuous-winner split.

No finite receiver test certifies a spatial integral or a general endpoint bound.
"""
from pathlib import Path
from fractions import Fraction as F
import datetime, hashlib, importlib.util, json, math, time

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'fixed_delay_winner_guard_20261007_results.json'
OLD = HERE / 'greedy_witness_cover_guard_20261007.py'
spec = importlib.util.spec_from_file_location('read_only_greedy_fixture', OLD)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)  # __main__ is NOT executed; prior outputs untouched.
ETA0 = F(49, 65536)


def digest(value):
    return hashlib.sha256(json.dumps(value, default=str, sort_keys=True).encode()).hexdigest()


def layered(n, theta):
    """A complete positive source with strict inside/boundary/outer layers.

    All atom masses, positions and lambda are fixed before evaluating receivers.
    This synthetic guard source may depend on the registered theta; receivers do not.
    """
    a = F(n)
    R = 2*a
    inside = theta*R*(1-F(1, 4*n))
    c = min((a/R)**n, (inside/R)**n)/16
    q = (inside/R)**n/3
    d = theta**n/6
    atoms = [engine.vector(n), engine.vector(n, inside/2),
             engine.vector(n, theta*R/2), engine.vector(n, R/2)]
    weights = [c, q, d, 1-c-q-d]
    assert all(w > 0 for w in weights) and sum(weights) == 1
    receivers = [dict(name='layer_origin', x=engine.vector(n)),
                 dict(name='layer_other_face_plus', x=engine.vector(n, 0, a/4)),
                 dict(name='layer_other_face_minus', x=engine.vector(n, 0, -a/4)),
                 dict(name='layer_shift_small', x=engine.vector(n, a/F(32*n))),
                 dict(name='layer_shift_negative', x=engine.vector(n, -a/F(32*n))),
                 dict(name='layer_near_outer_atom', x=engine.vector(n, a))]
    return a, atoms, weights, receivers, F(1)/(3*R**n)


def evaluate(n, theta, family, a, atoms, weights, receivers, lam):
    W = sum(weights)
    assert W == 1 and all(w > 0 for w in weights)
    rows = []
    for receiver in receivers:
        x = receiver['x']
        winner, candidates, arrivals = engine.true_window_winner(x, atoms, weights, n, a)
        R, m, u = winner['R'], winner['mass'], winner['u']
        inner = theta*R
        eligible = inner >= a  # beta >= -n log(theta), with equality allowed.
        band = 2*lam < u <= 4*lam
        opened = [j for j, arrival in enumerate(arrivals) if arrival < inner]
        closed = [j for j, arrival in enumerate(arrivals) if arrival <= inner]
        boundary = [j for j, arrival in enumerate(arrivals) if arrival == inner]
        open_mass = sum((weights[j] for j in opened), F(0))
        closed_mass = sum((weights[j] for j in closed), F(0))
        ratio = open_mass/m if m else None
        closed_ratio = closed_mass/m if m else None
        row = dict(name=receiver['name'], x=x, winner=winner,
                   all_continuous_candidates=candidates, atom_arrival_sides=arrivals,
                   hardband=band, inner_side=inner, inner_scale_allowed=eligible,
                   branch='long' if eligible else 'low',
                   inner_open_members=opened, inner_closed_members=closed,
                   inner_boundary_members=boundary, inner_open_mass=open_mass,
                   inner_closed_mass=closed_mass, open_mass_ratio=ratio,
                   closed_mass_ratio=closed_ratio, theta_power=theta**n,
                   raw_open_inequality_holds=(not m or open_mass <= theta**n*m),
                   raw_closed_inequality_holds=(not m or closed_mass <= theta**n*m))
        if eligible:
            # Independent direct inner response comparison. It is an allowed
            # continuum scale even when it is not one of the candidate arrivals.
            assert a <= inner <= 2*a
            inner_response = closed_mass/inner**n
            assert inner_response <= u
            assert open_mass <= closed_mass <= theta**n*m
            row.update(inner_response=inner_response,
                       inner_response_le_winner=True, long_tail_verified=True)
            if band:
                # This is the g=1 majorant of any actual traffic with 0<=g<=1.
                # Open inner is precisely delay > effective_v. Boundary equals
                # the cutoff and is excluded, without tolerance or tie splitting.
                traffic_majorant = open_mass/R**n
                assert traffic_majorant <= theta**n*u <= 4*theta**n*lam
                row.update(g_le_one_traffic_majorant=traffic_majorant,
                           hardband_tail_upper=4*theta**n*lam,
                           pointwise_hardband_tail_verified=True)
        else:
            assert R < a/theta
            row.update(low_branch_radius_upper=a/theta,
                       low_branch_radius_verified=True,
                       long_tail_verified=None,
                       pointwise_hardband_tail_verified=None)
        rows.append(row)
    counts = dict(receivers=len(rows), hardband=sum(r['hardband'] for r in rows),
                  long=sum(r['branch']=='long' for r in rows),
                  hardband_long=sum(r['branch']=='long' and r['hardband'] for r in rows),
                  low=sum(r['branch']=='low' for r in rows),
                  hardband_low=sum(r['branch']=='low' and r['hardband'] for r in rows),
                  low_raw_tail_failures=sum(r['branch']=='low' and not r['raw_open_inequality_holds'] for r in rows),
                  inner_boundary_rows=sum(bool(r['inner_boundary_members']) for r in rows),
                  full_winner_boundary_rows=sum(bool(r['winner']['boundary_members']) for r in rows),
                  inner_exactly_a_rows=sum(r['inner_side']==a for r in rows))
    core = dict(n=n, family=family, a=a, b=2*a, theta=theta, theta_power=theta**n,
                W=W, common_lambda=lam, atoms=atoms, weights=weights,
                receiver_records=rows, counts=counts,
                rational_short_source_once_budget_upper=(1+n*(1-theta)/theta)*W,
                rational_log_upper='-log(theta) <= (1-theta)/theta',
                exact_tail_coefficient=4*theta**n,
                target_tail_coefficient=4*ETA0,
                tail_coefficient_at_most_fixed_target=theta**n <= ETA0,
                effective_v_diagnostic_float=-n*math.log(float(theta)),
                effective_v_and_budget_scope='v=-n log(theta). The general source-once low-branch fee is (1+v)W: r<a contributes <=W and a<=r<=a/theta contributes <=v W. This general integral argument is recorded, not certified by finitely many receiver rows. The displayed rational fee is only an exact upper bound for (1+v)W.')
    core['input_sha256'] = digest(dict(n=n, a=a, theta=theta, atoms=atoms,
                                     weights=weights, lambda_=lam,
                                     receivers=receivers))
    return core


def finite_j_counterexample(n):
    a = F(n)
    b = 2*a
    theta = F(3, 4)
    x = engine.vector(n)
    atoms, weights = [x], [F(1)]
    lam = F(1)/(3*b**n)
    finite_R, finite_m = b, F(1)
    finite_u = finite_m/finite_R**n
    inner = theta*finite_R
    assert inner >= a and inner < b
    ratio = F(1)
    assert 2*lam < finite_u <= 4*lam
    assert ratio > theta**n
    assert finite_u > 4*theta**n*lam
    continuous, candidates, arrivals = engine.true_window_winner(x, atoms, weights, n, a)
    assert continuous['R'] == a
    return dict(n=n, a=a, b=b, theta=theta, full_source_mass=F(1), atoms=atoms,
                weights=weights, receiver=x, common_lambda=lam, finite_J=[b],
                finite_winner=dict(R=b, mass=finite_m, u=finite_u),
                inner_side=inner, inner_in_continuous_window=True, inner_in_finite_J=False,
                open_inner_captured_ratio=ratio, theta_power=theta**n,
                failed_tail_fraction_bound=True, g_one_traffic=finite_u,
                claimed_hardband_tail_upper=4*theta**n*lam,
                failed_hardband_traffic_bound=True,
                continuous_winner=continuous, continuous_candidates=candidates,
                explanation='J={b} forces R=b. All source mass is strictly inside theta*b, so ratio=1>theta^n. The inner scale belongs to [a,b] but is absent from J. Continuous maximization instead chooses R=a, a low-branch row; no contradiction to the continuous-winner guard.')


def main():
    start = time.perf_counter()
    records = []
    for round_, (n, target_theta) in enumerate([(4,F(1,8)), (16,F(5,8)), (64,F(7,8))], 1):
        # target-theta tails <=49/65536 exactly; theta=3/4 is a separate
        # registered diagnostic permitting nonempty long rows at n=4.
        for theta_label, theta in [('target_coefficient_guard',target_theta),
                                   ('moderate_shrink_diagnostic',F(3,4)),
                                   ('inner_endpoint_guard',F(1,2))]:
            assert theta_label != 'target_coefficient_guard' or theta**n <= ETA0
            a, h, atoms, weights, boxes, receivers, lam = engine.construct(n)
            for family, fixture in [('read_only_overlap_fixture',(a,atoms,weights,receivers,lam)),
                                    ('registered_layered_fixture',layered(n,theta))]:
                record = evaluate(n,theta,family,*fixture)
                record.update(round=round_,theta_label=theta_label)
                records.append(record)
                print('ROUND',round_,'n',n,theta_label,family,record['counts'],flush=True)
    counterexamples = [finite_j_counterexample(n) for n in [4,16,64]]
    counts = {key:sum(r['counts'][key] for r in records) for key in records[0]['counts']}
    assert counts['hardband_long'] > 0 and counts['hardband_low'] > 0
    assert counts['low_raw_tail_failures'] > 0 and counts['inner_boundary_rows'] > 0
    assert counts['inner_exactly_a_rows'] > 0
    data = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                records=records,finite_J_singleton_counterexamples=counterexamples,counts=counts,
                arithmetic='Exact Fraction coordinates, weights, arrival sides, full/inner open and closed membership, candidate responses, winners, common lambdas, classification and power comparisons. Floats occur only in optional effective_v diagnostics, timestamps and runtime.',
                fixed_target=dict(exp_minus_vstar=ETA0,exp_minus_vstar_times_four=4*ETA0,
                                  vstar='log(65536/49)',
                                  n4_target_long_branch_empty_reason='2^4*49/65536<1, so vstar>4 log2 and beta<=4 log2 in [a,2a]. This is checked by exact rational comparison.',
                                  n4_target_long_impossible_exact=(2**4*ETA0 < 1)),
                winner_scope='Exact REAL continuous-window maximal responses for finite rational positive atoms: endpoints a,b and every atom arrival side within [a,b]; complete equal-distance groups; response ties choose smaller R. Direct membership agrees with independently grouped cumulative mass. This is not a finite radius grid and is not an original finite-J certification.',
                branch_scope='Long iff theta*R>=a; inner<=b automatically. At long, closed inner response<=winner response implies closed mass/m<=theta^n; open inner excludes every equal-delay boundary atom. At low, no tail bound is inferred even when its raw inequality happens to hold. The radius guard R<a/theta is recorded.',
                traffic_scope='The recorded open_inner_mass/R^n is a source-weighted pointwise g=1 majorant for 0<=actual g<=1. It is not numerical actual FIRST, not a fitted kernel, and not Gamma energy. Hardband 2lambda<u<=4lambda gives traffic majorant<=4theta^n lambda. Integrating this pointwise inequality on E is the separately audited general traffic argument.',
                short_branch_scope='The symbolic source-once budget (1+v)W, v=-nlog(theta), uses dx=r^n dt dnu for normalized cone measure on [-1/2,1/2]^n. r<a pays W via (r/a)^n, a<=r<=a/theta pays v W. Finitely many receiver tests do not certify the spatial integral or all-source geometry.',
                source_scope='All positive masses retained. Each fixture has one lambda for all its receivers and it is never retuned per receiver. The registered layered fixture depends on n and theta solely to exercise exact edge cases. The old fixture is imported read-only without executing its main, and its source and lambda are unchanged.',
                confidence_scope='Deterministic exact finite-input certificates; no Monte Carlo, confidence interval or asymptotic regression. No implication of a general geom endpoint bound.',
                imported_script_sha256=hashlib.sha256(OLD.read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                runtime_seconds=time.perf_counter()-start)
    with OUTPUT.open('x') as f:
        json.dump(data,f,default=str,ensure_ascii=False,indent=2)
    print('OUTPUT',OUTPUT,'COUNTS',counts,'RUNTIME',data['runtime_seconds'],flush=True)


if __name__ == '__main__':
    main()
