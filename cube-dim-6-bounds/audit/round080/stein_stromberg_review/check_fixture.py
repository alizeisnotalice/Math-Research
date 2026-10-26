"""Four fixed independent diagnostic checks, not a rerun of 080's five facts."""
import sys, signal, time, resource
sys.dont_write_bytecode = True
START = time.monotonic()
CPU = time.process_time()
resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
def halt(sig, frame):
    raise TimeoutError('A080 work limit; save once and stop')
signal.signal(signal.SIGALRM, halt)
signal.setitimer(signal.ITIMER_REAL, 169)
from pathlib import Path
from fractions import Fraction as Q
import json, os, hashlib, traceback
import sympy as sp

H = Path(__file__).resolve().parent
read = lambda p: json.loads(Path(p).read_text())
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name, data):
    (H / name).write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
if (H / 'run_receipt.json').exists():
    raise RuntimeError('One invocation only; saved receipt exists')
reg = read(H / 'registration.json')
A = Path(reg['author_directory'])
out = {'status': 'INCOMPLETE', 'configurations': [],
       'fixture_dimension': 1, 'fixture_count': 1, 'fixture_ball_count': 5,
       'author_five_fact_reevaluations': 0, 'parameter_grid_points': 0,
       'body32_kernel_calls': 0, 'old_baseline_evaluations': 0}
receipt = {'status': 'RUNNING', 'pid': os.getpid(),
           'attempted_configurations': 0, 'completed_configurations': 0,
           'scientific_program_invocations': 1, 'automatic_retries': 0,
           'CPU_limit_seconds': 30, 'wall_limit_including_cleanup_seconds': 180}
def checkpoint():
    receipt['wall_seconds'] = time.monotonic() - START
    receipt['cpu_seconds'] = time.process_time() - CPU
    save('results.json', out)
    save('run_receipt.json', receipt)
def begin(i):
    assert i == receipt['attempted_configurations'] + 1 <= 4
    receipt['attempted_configurations'] = i
    checkpoint()
def done(i, **kw):
    out['configurations'].append({'id': i, 'pass': True, **kw})
    receipt['completed_configurations'] = i
    checkpoint()
def pins():
    for p, digest in reg['source_sha256'].items():
        assert sha(p) == digest, p
    assert sha(reg['source_pdf']) == reg['source_pdf_sha256']
    for f, record in read(A / 'manifest.json')['files'].items():
        assert sha(A / f) == record['sha256'], f
        assert (A / f).stat().st_size == record['bytes'], f
    bind = read(H / 'execution_binding.json')
    assert sha(H / 'registration.json') == bind['registration_sha256']
    assert sha(H / 'proof_review.md') == bind['proof_review_sha256']
    for f, digest in bind['code_sha256'].items():
        assert sha(H / f) == digest, f

def dilation(ball, factor):
    return ball['center'] - factor * ball['radius'], ball['center'] + factor * ball['radius']
def inside(x, interval):
    return interval[0] < x < interval[1]
def measure(segments):
    return sum((right - left for left, right in segments), Q(0))
def subtract(interval, previous):
    pieces = [interval]
    for cover in previous:
        next_pieces = []
        for left, right in pieces:
            cl, cr = cover
            if cr <= left or right <= cl:
                next_pieces.append((left, right))
            else:
                if left < cl:
                    next_pieces.append((left, min(right, cl)))
                if cr < right:
                    next_pieces.append((max(left, cr), right))
        pieces = next_pieces
    return pieces
def union_measure(intervals):
    covered = []
    total = Q(0)
    for interval in intervals:
        total += measure(subtract(interval, covered))
        covered.append(interval)
    return total

try:
    begin(1)
    pins()
    author = read(A / 'results.json')
    assert author['status'] == 'PASS_FIVE_FIXED_EXACT_FACTS'
    assert author['configurations'] == 5
    assert author['code_sha256'] == sha(A / 'check_exact.py')
    assert author['theorem_sha256'] == sha(A / 'theorem.md')
    assert read(A / 'run_receipt.json')['exit'] == 0
    # These are dependencies read from saved data, not recomputed facts.
    out['author_saved_fact_strings_read_only'] = {k: author[k] for k in (
        'exp_prefix', 'log_two_prefix', 'log_seven_prefix',
        'middle_product_lower', 'uniform_allocated_cost_strict_lower')}
    beta = Q(reg['fixture']['beta'])
    family = [{'id': b['id'], 'center': Q(b['center']), 'radius': Q(b['radius'])}
              for b in reg['fixture']['family_in_tie_order']]
    selected, remaining, trace = [], list(family), []
    def score(center, collection):
        return sum((b['weight'] for b in collection if inside(center, dilation(b, beta + 1))), Q(0))
    while remaining:
        eligible = [b for b in remaining if score(b['center'], selected) <= 1]
        if not eligible:
            break
        # Python max keeps the first in the frozen tie order.
        chosen = max(eligible, key=lambda b: b['radius'])
        admission = score(chosen['center'], selected)
        increments = subtract(dilation(chosen, Q(1)), [dilation(b, Q(1)) for b in selected])
        increment_volume = measure(increments)
        chosen.update(admission=admission, increments=increments,
                      increment_volume=increment_volume,
                      weight=increment_volume / (2 * beta * chosen['radius']))
        trace.append({'chosen': chosen['id'], 'radius': str(chosen['radius']),
                      'admission_score': str(admission),
                      'ineligible_before_choice': [b['id'] for b in remaining if score(b['center'], selected) > 1],
                      'increment_measure': str(increment_volume), 'weight': str(chosen['weight'])})
        selected.append(chosen)
        remaining.remove(chosen)
    assert [b['id'] for b in selected] == ['A', 'B', 'D', 'E']
    assert [b['id'] for b in remaining] == ['C']
    assert [b['radius'] for b in selected] == sorted([b['radius'] for b in selected], reverse=True)
    assert [b['increment_volume'] for b in selected] == [Q(8), Q(8), Q(1), Q(2)]
    assert selected[-1]['admission'] == 1
    assert 'C' in trace[2]['ineligible_before_choice']
    assert all(score(b['center'], selected) > 1 for b in remaining)
    out['greedy_trace'] = trace
    out['selected_increment_intervals'] = {
        b['id']: [[str(left), str(right)] for left, right in b['increments']] for b in selected}
    done(1, name='finite_greedy_exact_selection_and_persistence',
         selected=['A', 'B', 'D', 'E'], unselected=['C'],
         blocked_radius='3', subsequently_selected_radius='1',
         exact_boundary_admission='1', disjoint_increment_total='19')

    begin(2)
    cover_certificates = []
    for b in remaining:
        large = [j for j in selected if j['radius'] >= b['radius']
                 and inside(b['center'], dilation(j, beta + 1))]
        assert sum((j['weight'] for j in large), Q(0)) > 1
        bl, br = dilation(b, Q(1))
        for j in large:
            jl, jr = dilation(j, beta + 2)
            assert jl < bl < br < jr
        cover_certificates.append({'unselected': b['id'], 'large_controllers': [j['id'] for j in large],
                                   'large_weight_sum': str(sum((j['weight'] for j in large), Q(0)))})
    x = Q(reg['fixture']['multiplicity_point'])
    contributing = [j for j in selected if inside(x, dilation(j, beta))]
    minimum_radius = min(j['radius'] for j in contributing)
    assert minimum_radius == 1
    k = next(j for j in reversed(selected) if j in contributing and j['radius'] == minimum_radius)
    outer = [j for j in contributing if j['radius'] >= beta]
    assert all(selected.index(j) < selected.index(k) for j in outer)
    assert all(inside(k['center'], dilation(j, beta + 1)) for j in outer)
    outer_weight = sum((j['weight'] for j in outer), Q(0))
    assert outer_weight == Q(2, 3) <= k['admission'] <= 1
    shell = [j for j in contributing if 1 <= j['radius'] < beta]
    assert len(outer) + len(shell) == len(contributing)
    shell_weight = sum((j['weight'] for j in shell), Q(0))
    enclosing_radius = (beta + 1) * beta
    for j in shell:
        assert all(x - enclosing_radius < left < right < x + enclosing_radius for left, right in j['increments'])
    assert shell_weight == 1
    shell_volume_bound = (beta + 1) / beta * beta
    assert shell_weight <= shell_volume_bound
    b, r = sp.symbols('beta r', positive=True)
    assert sp.expand((b + 1) * r - (b * r + b)) == r - b
    assert sp.expand(((b + 1) * r - (b * r + b)).subs(r, b)) == 0
    out['unselected_cover_certificates'] = cover_certificates
    done(2, name='outer_radius_and_disjoint_shell_certificate',
         x=str(x), minimum_selected='E', outer_ids=[j['id'] for j in outer],
         outer_weight=str(outer_weight), shell_ids=[j['id'] for j in shell],
         shell_weight=str(shell_weight), shell_volume_upper=str(shell_volume_bound),
         threshold_identity='(beta+1)r-(beta*r+beta)=r-beta',
         beta_boundary_in_outer=True, spatial_sampling_grid=False)

    begin(3)
    b, m, a, u, v, c = sp.symbols('beta m a u v c', positive=True)
    derivative = sp.diff(sp.exp(u) / u, u)
    assert sp.simplify(derivative - sp.exp(u) * (u - 1) / u ** 2) == 0
    assert sp.simplify(m * sp.exp(a / m) - a * sp.exp(a / m) / (a / m)) == 0
    first = 1 + ((b + 2) / b) ** 6
    assert sp.simplify(sp.diff(first, b) + 12 * (b + 2) ** 5 / b ** 7) == 0
    volume_ratio = c * ((b + 1) * v) ** 6 / (c * (b * u) ** 6)
    assert sp.simplify(volume_ratio - ((b + 1) / b) ** 6 * (v / u) ** 6) == 0
    shell_ratio = b ** (1 / m)
    assembled = (1 + (b + 2) ** 6 / b ** 6) * (1 + m * ((b + 1) / b) ** 6 * shell_ratio ** 6)
    target = (1 + ((b + 2) / b) ** 6) * (1 + m * ((b + 1) / b) ** 6 * b ** (6 / m))
    assert sp.simplify(assembled - target) == 0
    done(3, name='independent_symbolic_all_parameter_structure',
         derivative='exp(u)*(u-1)/u^2',
         minimization='a*exp(u)/u >= a*e for u=a/m>0',
         first_factor_derivative='-12*(beta+2)^5/beta^7',
         shell_volume_ratio_residual='0', allocated_six_dimensional_formula_residual='0',
         author_five_rational_facts_reevaluated=0,
         parameter_grid_points=0, beta_equals_one_claim=False)

    begin(4)
    support = (Q(-10), Q(12))
    lam = Q(reg['fixture']['lambda'])
    def step_integral(interval):
        return max(Q(0), min(support[1], interval[1]) - max(support[0], interval[0]))
    witnesses = []
    for j in family:
        star = dilation(j, beta)
        average = step_integral(star) / (star[1] - star[0])
        assert average > lam
        witnesses.append({'id': j['id'], 'small_radius': str(j['radius']),
                          'witness_radius': str(beta * j['radius']), 'average': str(average)})
    weighted_integral = sum((j['weight'] * step_integral(dilation(j, beta)) for j in selected), Q(0))
    increment_total = sum((j['increment_volume'] for j in selected), Q(0))
    assert weighted_integral == increment_total == Q(19)
    assert lam * increment_total < weighted_integral
    full_mass = support[1] - support[0]
    cover_factor = 1 + (beta + 2) / beta
    multiplicity_factor = 1 + (beta + 1) / beta * beta
    assert weighted_integral <= multiplicity_factor * full_mass
    original_union_volume = union_measure([dilation(j, Q(1)) for j in family])
    assert original_union_volume <= cover_factor * increment_total
    assert lam * original_union_volume <= cover_factor * multiplicity_factor * full_mass
    out['common_scale_witnesses'] = witnesses
    done(4, name='exact_weighted_L1_witness_transfer_on_same_fixture',
         function='indicator[-10,12]', function_L1=str(full_mass), lambda_value=str(lam),
         weighted_integral=str(weighted_integral), increment_total=str(increment_total),
         original_union_measure=str(original_union_volume),
         restricted_to_finite_fixture=True,
         general_L1_and_all_scales='proved in proof_review.md, not inferred from fixture')
    pins()
    out.update(status='PASS_SPECIFIED_ALLOCATION_DIAGNOSTIC_ONLY',
               author_manifest_sha256=sha(A / 'manifest.json'),
               author_theorem_sha256=sha(A / 'theorem.md'),
               analytical_beta_domain='beta>1', analytical_m_domain='integer m>=1',
               specified_formula_strict_lower='1602/25', new_C6_bound=False,
               true_covering_optimum_lower_bound=False, general_covering_method_excluded=False,
               source_five_facts_read_only=True,
               fixture_scope='one 1D five-interval family; not numerical evidence for all 6D bodies')
    receipt.update(status='PASS_STOPPED', exit_code=0)
    checkpoint()
    receipt['results_sha256'] = sha(H / 'results.json')
    save('run_receipt.json', receipt)
    print(json.dumps({'status': out['status'], 'configurations': 4,
                      'cpu_seconds': receipt['cpu_seconds'], 'wall_seconds': receipt['wall_seconds'],
                      'fixture_count': 1, 'author_five_fact_reevaluations': 0}))
except BaseException as exc:
    receipt.update(status='FAILED_STOP_NO_RETRY', exit_code=1,
                   error=repr(exc), traceback=traceback.format_exc())
    checkpoint()
    raise
finally:
    signal.setitimer(signal.ITIMER_REAL, 0)
