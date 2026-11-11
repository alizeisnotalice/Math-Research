"""Independent acceptance for r090 scalar Taylor remainder certificates."""
import argparse
import hashlib
import json
import multiprocessing as mp
import re
import sys
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
from pathlib import Path

REL_TOL = D('1e-45')
K_VALUES = (16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024)
COMPARE_FIELDS = ('mu', 'sigma', 'L', 'x_max', 'log_remainder', 'log_cubic_bound',
                  'remainder_to_cubic_bound', 'z', 'h', 'delta0', 'delta0_over_h',
                  'normalized_identity_value')


def number(x):
    if x is None:
        raise AssertionError('unexpected null numeric field')
    return D(re.sub(r'`[^*]*', '', str(x)).replace('*^', 'E'))


def finite(x):
    z = number(x)
    assert z.is_finite(), x
    return z


def relative_difference(a, b):
    a, b = finite(a), finite(b)
    if a == 0 and b == 0:
        return D(0)
    denominator = max(abs(a), abs(b))
    assert denominator > 0
    return abs(a-b)/denominator


def qdecimal(x):
    q = Q(x)
    return D(q.numerator)/D(q.denominator)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_replay(stored, replay):
    assert len(stored) == len(replay)
    for old, new in zip(stored, replay):
        assert {k: v for k, v in old.items() if k != 'pid'} == {
            k: v for k, v in new.items() if k != 'pid'
        }


def resource_check(z):
    assert z['exit_code'] == 0 and z['stop_reason'] is None
    assert D(str(z['reserve_GiB'])) >= D(1)
    reserve = D(str(z['reserve_GiB']))*D(1024)
    budget = D(str(z['owned_RSS_budget_GiB']))*D(1024)
    assert D(str(z['min_system_available_MiB'])) >= reserve
    assert D(str(z['peak_owned_RSS_MiB'])) <= budget
    assert 'estimate' in z['scope'].lower()
    return {
        'exit_code': z['exit_code'], 'stop_reason': z['stop_reason'],
        'reserve_GiB': z['reserve_GiB'], 'owned_RSS_budget_GiB': z['owned_RSS_budget_GiB'],
        'peak_owned_RSS_MiB': z['peak_owned_RSS_MiB'],
        'min_system_available_MiB': z['min_system_available_MiB'],
        'samples': z['samples'], 'seconds': z['seconds'],
        'scope': z['scope'], 'available_memory_is_estimate': True,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory')
    ap.add_argument('output')
    ap.add_argument('--project-root', required=True)
    args = ap.parse_args()
    directory, project_root = Path(args.directory), Path(args.project_root)
    code_dir = Path(__file__).resolve().parent
    names = ('python/experiments90.json', 'python/resources90.json',
             'wolfram/verification90.json', 'wolfram/resources90.json')
    for name in names:
        assert (directory/name).is_file(), directory/name
    py = json.loads((directory/names[0]).read_text())
    pyres = json.loads((directory/names[1]).read_text())
    wf = json.loads((directory/names[2]).read_text())
    wfres = json.loads((directory/names[3]).read_text())

    project_code = project_root/'code'
    sys.path.insert(0, str(project_code))
    sys.path.insert(0, str(code_dir))
    from experiment90 import job
    from verify85 import constants85, sqrt_interval

    assert py['pass_all'] is True and wf['pass'] is True
    assert py['working_precision'] == [80, 160]
    assert len(wf['symbolic']) == 12 and all(v is True for v in wf['symbolic'].values())
    assert wf['kernels_closed'] is True
    assert 1 <= wf['configured_workers'] <= 12
    assert py['configured_workers'] == 12
    assert py['counts'] == {'dimensions': 13, 'numeric_precision_rows': 26, 'grid_points_enumerated': 0}
    assert py['actual_workers'] == len(set(py['worker_pids'])) and 0 < py['actual_workers'] <= 12
    assert 0 < wf['per_case_memory_limit_MiB'] <= 1024
    pyresource = resource_check(pyres)
    wfresource = resource_check(wfres)

    expected_mu = constants85()['mu_interval']
    assert py['constants']['mu_interval'] == expected_mu
    sigma_iv = sqrt_interval(tuple(Q(x) for x in expected_mu))
    assert py['constants']['sigma_interval'] == [str(x) for x in sigma_iv]
    pyrows, wfrows = py['rows'], wf['rows']
    assert len(pyrows) == len(wfrows) == len(K_VALUES)
    assert [r['k'] for r in pyrows] == list(K_VALUES)
    assert [r['k'] for r in wfrows] == list(K_VALUES)
    assert [r['N'] for r in pyrows] == [str(1 << k) for k in K_VALUES]

    # Replay the full exact/numeric Python job set using a fresh spawn pool.
    jobs = [(k, expected_mu) for k in K_VALUES]
    with mp.get_context('spawn').Pool(12) as pool:
        replay = pool.map(job, jobs)
    compare_replay(pyrows, replay)
    replay_workers = len({r['pid'] for r in replay})

    kernel_ids = {r['kernel_id'] for r in wfrows}
    assert kernel_ids and all(isinstance(k, int) and k > 0 for k in kernel_ids)
    assert len(kernel_ids) <= wf['configured_workers']

    same_precision_diffs = {'python': [], 'wolfram': []}
    cross_diffs = []
    field_max = {field: D(0) for field in COMPARE_FIELDS}
    max_remainder_interval_width = D(0)
    max_cubic_interval_width = D(0)
    grid_counts = []
    max_precision_residual = D(0)
    for pr, wr in zip(pyrows, wfrows):
        k = pr['k']
        assert wr['k'] == k
        assert pr['uniform_remainder_is_error_bound'] is True
        assert pr['grid_enumerated'] is False
        assert pr['grid_count_upper'] < 1 << k
        grid_counts.append(pr['grid_count_upper'])
        py_prec, wf_prec = pr['precision_rows'], wr['precision_rows']
        assert [r['working_precision'] for r in py_prec] == [80, 160]
        assert [r['working_precision'] for r in wf_prec] == [80, 160]
        rem_lo, rem_hi = map(qdecimal, pr['remainder_log_interval'])
        cubic_lo, cubic_hi = map(qdecimal, pr['cubic_bound_log_interval'])
        assert rem_lo <= rem_hi and cubic_lo <= cubic_hi
        max_remainder_interval_width = max(max_remainder_interval_width, rem_hi-rem_lo)
        max_cubic_interval_width = max(max_cubic_interval_width, cubic_hi-cubic_lo)
        assert pr['uniform_remainder_log_upper'] == pr['remainder_log_interval'][1]
        for sw_index, rows in enumerate((py_prec, wf_prec)):
            for wp_row in rows:
                lr = finite(wp_row['log_remainder'])
                lc = finite(wp_row['log_cubic_bound'])
                ratio = finite(wp_row['remainder_to_cubic_bound'])
                assert rem_lo <= lr <= rem_hi, ('remainder interval', k, wp_row['working_precision'], lr)
                assert ratio > 0
                assert ratio <= D(1)+D('1e-40')
                assert cubic_lo <= lc <= cubic_hi
                assert abs(finite(wp_row['normalized_identity_value'])-1) < REL_TOL
                assert D(0) < finite(wp_row['x_max']) < D('.25')
                assert finite(wp_row['h']) > 0 and finite(wp_row['delta0']) > 0
                for field in COMPARE_FIELDS:
                    finite(wp_row[field])
            for field in COMPARE_FIELDS:
                diff = relative_difference(rows[0][field], rows[1][field])
                assert diff < REL_TOL, (k, 'within software precision', field, diff)
                same_precision_diffs['python' if sw_index == 0 else 'wolfram'].append(diff)
        for pi in range(2):
            p_row, w_row = py_prec[pi], wf_prec[pi]
            assert w_row['pass'] is True
            assert finite(w_row['log_tail_residual']) < D(10)**(-w_row['working_precision'])
            for field in COMPARE_FIELDS:
                diff = relative_difference(p_row[field], w_row[field])
                assert diff < REL_TOL, (k, 'cross software', p_row['working_precision'], field, diff)
                cross_diffs.append(diff)
                field_max[field] = max(field_max[field], diff)
            max_precision_residual = max(max_precision_residual,
                                         finite(w_row['log_tail_residual']))

    # Numerical statistics and provenance hashes.
    local_sources = ('experiment90.py', 'verify90.wl', 'check90.py')
    public_dependencies = ('math_guard.py', 'verify85.py', 'verify86.py', 'verify84.py',
                           'sliver79.py', 'experiment88.py')
    source_hashes = {f: sha((project_code/f)
                            if f in public_dependencies else (code_dir/f))
                     for f in local_sources+public_dependencies}
    input_hashes = {f: sha(directory/f) for f in names}
    document_hashes = {f: sha(directory/f) for f in ('core90.md', 'audit90.md')}
    prior = ['rounds/r011/report.md', 'rounds/r021/report.md', 'rounds/r022/report.md', 'rounds/r079/rational79.md', 'rounds/r089/core89.md', 'rounds/r089/audit89.md']
    prior_hashes = {f: sha(project_root/f) for f in prior}
    out = {
        'pass_all': True,
        'scope': 'Exact scalar Taylor-remainder bounds and scalar diagnostics only; no numerical proof of the physical witness theorem and no finite-N physical onset.',
        'case_counts': {'dimensions': len(pyrows), 'numeric_precision_rows': 26,
                        'symbolic_gates': len(wf['symbolic']), 'grid_points_enumerated': 0},
        'dimensions_k': list(K_VALUES),
        'Python_configured_workers': py['configured_workers'],
        'Python_actual_workers_stored': py['actual_workers'],
        'Python_replay_actual_workers': replay_workers,
        'Python_replay_worker_pids': sorted({r['pid'] for r in replay}),
        'Wolfram_configured_workers': wf['configured_workers'],
        'Wolfram_actual_kernel_ids': sorted(kernel_ids),
        'Wolfram_kernels_closed': wf['kernels_closed'],
        'Wolfram_symbolic_gates': wf['symbolic'],
        'working_precision': [80, 160],
        'relative_tolerance': str(REL_TOL),
        'max_within_python_precision_difference': str(max(same_precision_diffs['python'], default=D(0))),
        'max_within_wolfram_precision_difference': str(max(same_precision_diffs['wolfram'], default=D(0))),
        'max_cross_software_relative_difference': str(max(cross_diffs, default=D(0))),
        'max_cross_software_difference_by_field': {k: str(v) for k, v in field_max.items()},
        'max_remainder_log_interval_width': str(max_remainder_interval_width),
        'max_cubic_log_interval_width': str(max_cubic_interval_width),
        'max_wolfram_log_tail_residual': str(max_precision_residual),
        'grid_count_upper_range': [min(grid_counts), max(grid_counts)],
        'grid_enumerated': False,
        'Python_resource_monitor': pyresource,
        'Wolfram_resource_monitor': wfresource,
        'available_memory_is_estimated_not_hard_reserved': True,
        'source_sha256': source_hashes,
        'input_sha256': input_hashes,
        'document_sha256': document_hashes,
        'prior_dependency_sha256': prior_hashes,
        'goal_status': 'active',
        'pause_request_revoked': True,
    }
    Path(args.output).write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: v for k, v in out.items() if not k.endswith('sha256')}, indent=2))


if __name__ == '__main__':
    with localcontext() as ctx:
        ctx.prec = 180
        main()
