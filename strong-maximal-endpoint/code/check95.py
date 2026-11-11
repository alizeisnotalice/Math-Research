"""Portable exact replay and two-language validation for r095 checks."""
import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import sys
from fractions import Fraction as Q
from pathlib import Path

EXPECTED_FILES = (
    'python/certificates95.json', 'python/resources95.json',
    'wolfram/verification95.json', 'wolfram/resources95.json',
)
OWN_SOURCES = ('experiment95.py', 'verify95.wl', 'check95.py')
P_VALUES = (Q(1, 31), Q(1, 127))


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve_file(name, data_dir, project_root, extra=()):
    candidates = [data_dir/name, data_dir.parent/name,
                  project_root/'work/r095'/name,
                  project_root/'rounds/r095'/name,
                  project_root/'code'/name,
                  project_root/'public-code'/name,
                  project_root/'public_code'/name,
                  project_root/'code/public'/name,
                  project_root/'public'/name]
    candidates.extend(extra)
    for path in candidates:
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError(f'could not resolve {name}; checked portable source/data roots')


def resource(row):
    assert row['exit_code'] == 0 and row['stop_reason'] is None
    assert float(row['reserve_GiB']) >= 1
    assert float(row['min_system_available_MiB']) >= 1024*float(row['reserve_GiB'])
    assert float(row['peak_owned_RSS_MiB']) <= 1024*float(row['owned_RSS_budget_GiB'])
    assert row['samples'] > 0
    return row


def import_source(path):
    spec = importlib.util.spec_from_file_location('experiment95_replay', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def q(value):
    return Q(str(value))


def check_overlap(py, wf):
    assert py['patterns_enumerated'] == wf['patterns_enumerated'] == 64
    assert py['pattern_rows_sha256'] == wf['pattern_rows_sha256']
    for field in ('raw_mass_expectation', 'surviving_mass_expectation',
                  'lost_mass_expectation', 'expected_ordered_pair_bound'):
        assert q(py[field]) == q(wf[field])
    assert q(py['raw_mass_expectation']) == Q(6, 31)
    assert q(py['surviving_mass_expectation']) == Q(148230000, 887503681)
    assert q(py['lost_mass_expectation']) == Q(23544906, 887503681)
    assert q(py['expected_ordered_pair_bound']) == Q(30, 961)
    assert q(py['lost_mass_expectation']) <= q(py['expected_ordered_pair_bound'])
    assert py['exact_loss_le_pair_bound'] is wf['exact_loss_le_pair_bound'] is True
    assert py['pointwise_loss_le_l_times_l_minus_1'] is wf['pointwise_loss_le_l_times_l_minus_1'] is True


def check_constants(py, wf):
    fields = ('c', 'eta', 'delta', 'lambda_upper',
              'minimum_expected_success_fraction',
              'packing_simultaneous_success_probability_lower',
              'guard_simultaneous_success_probability_lower',
              'long_edge_deleted_fraction_lower')
    for field in fields:
        assert q(py[field]) == q(wf[field]), (field, py[field], wf[field])
    assert q(py['minimum_expected_success_fraction']) == Q(20455, 32768)
    assert q(py['packing_simultaneous_success_probability_lower']) == Q(485691, 1048576)
    assert q(py['guard_simultaneous_success_probability_lower']) == Q(45, 64)
    assert q(py['long_edge_deleted_fraction_lower']) == Q(3, 256)
    assert q(py['long_edge_deleted_fraction_lower']) > q(py['eta'])
    assert py['strictly_exceeds_eta'] is wf['strictly_exceeds_eta'] is True
    assert py['derived_with_exact_fraction_arithmetic'] is wf['derived_with_exact_fraction_arithmetic'] is True


def check_generalization(py, wf):
    assert py['case_count'] == wf['case_count'] == 256
    assert py['canonical_case_rows_sha256'] == wf['canonical_case_rows_sha256']
    assert py['loss_le_pair_bound_all_cases'] is wf['loss_le_pair_bound_all_cases'] is True
    for field in ('minimum_surviving_mass', 'maximum_surviving_mass',
                  'minimum_pair_bound_minus_loss', 'maximum_loss_over_pair_bound'):
        a, b = py[field], wf[field]
        assert (q(a['value']), a['m'], a['k'], q(a['p'])) == \
               (q(b['surviving'] if 'surviving' in b else b['value']),
                b['m'], b['k'], q(b['p']))
    assert q(py['minimum_pair_bound_minus_loss']['value']) >= 0
    assert 0 <= q(py['maximum_loss_over_pair_bound']['value']) <= 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory', help='r095 directory containing python/, wolfram/, core95.md and audit95.md')
    ap.add_argument('output')
    ap.add_argument('--project-root', required=True,
                    help='repository root; source resolver also accepts code/, public-code/ and rounds/r095/')
    args = ap.parse_args()
    data_dir = Path(args.directory).resolve()
    project_root = Path(args.project_root).resolve()

    exp_path = resolve_file('experiment95.py', data_dir, project_root)
    module = import_source(exp_path)
    py, pyres, wf, wfres = [json.loads((data_dir/name).read_text()) for name in EXPECTED_FILES]
    pyres, wfres = resource(pyres), resource(wfres)
    assert py['pass_all'] is True and wf['pass'] is True
    assert py['overlap']['patterns_enumerated'] == 64
    assert py['generalization']['case_count'] == 256
    assert py['configured_workers'] in (1, 2)
    assert py['actual_workers'] <= py['configured_workers']
    assert len(py['worker_pids']) == py['actual_workers']
    assert wf['configured_workers'] in (1, 2) and wf['kernels_closed'] is True
    assert len(wf['actual_kernel_ids']) <= wf['configured_workers']
    assert all(int(kid) > 0 for kid in wf['actual_kernel_ids'])

    check_overlap(py['overlap'], wf['overlap'])
    check_constants(py['constants'], wf['constants'])
    check_generalization(py['generalization'], wf['generalization'])
    contract = py['log_contract']
    wf_contract = wf['log_contract']
    assert contract['statement'] == wf_contract['statement']
    assert contract['domain'] == wf_contract['domain']
    assert contract['numeric_log_evaluation_used'] is False
    assert wf_contract['numeric_log_evaluation_used'] is False
    facts = wf_contract['exact_ordering_facts']
    assert len(facts) == 4 and all(facts.values())

    # Replay the actual exact Python entry points, then compare every
    # deterministic result to the stored sections. Worker PIDs are irrelevant.
    replay = {name: (value, pid) for name, value, pid in
              map(module.run_job, ('overlap', 'constants', 'generalization'))}
    for key in ('overlap', 'constants', 'generalization'):
        stored = py[key]
        replayed = replay[key][0]
        assert stored == replayed, key

    own_hashes = {}
    for name in OWN_SOURCES:
        path = resolve_file(name, data_dir, project_root,
                            extra=(Path(__file__).resolve().parent/name,))
        own_hashes[name] = file_hash(path)
    core_path = resolve_file('core95.md', data_dir, project_root)
    audit_path = resolve_file('audit95.md', data_dir, project_root)
    document_hashes = {'core95.md': file_hash(core_path), 'audit95.md': file_hash(audit_path)}
    input_hashes = {name: file_hash(data_dir/name) for name in EXPECTED_FILES}
    for name in ('report.md','certification95.md','sources95.md'):
        if (data_dir/name).is_file(): document_hashes[name]=file_hash(data_dir/name)
    prior_names=('rounds/r021/report.md','rounds/r091/core91.md','rounds/r092/core92.md','rounds/r093/core93.md','rounds/r094/core94.md','rounds/r094/audit94.md')
    prior_hashes={name:file_hash(project_root/name) for name in prior_names if (project_root/name).is_file()}

    result = {
        'pass_all': True,
        'goal_status': 'active',
        'general_upper_closed': False,
        'prior_dependency_sha256': prior_hashes,
        'scope': ('Exact finite point-model checks and the symbolic logarithmic boundary contract only; '
                  'no high-dimensional physical source construction, no general H_N coefficient, '
                  'and no finite-dimensional onset.'),
        'overlap_pattern_count': 64,
        'generalization_case_count': 256,
        'Python_configured_workers': py['configured_workers'],
        'Python_stored_workers': py['actual_workers'],
        'Python_replay_workers': len({row[1] for row in replay.values()}),
        'Wolfram_configured_workers': wf['configured_workers'],
        'Wolfram_kernel_ids': wf['actual_kernel_ids'],
        'Wolfram_kernels_closed': wf['kernels_closed'],
        'Python_resources': pyres,
        'Wolfram_resources': wfres,
        'resources_are_sampled_estimates_not_OS_reservations': True,
        'overlap_pattern_rows_sha256': py['overlap']['pattern_rows_sha256'],
        'generalization_case_rows_sha256': py['generalization']['canonical_case_rows_sha256'],
        'constants': py['constants'],
        'generalization_extrema': {key: py['generalization'][key] for key in
                                   ('minimum_surviving_mass', 'maximum_surviving_mass',
                                    'minimum_pair_bound_minus_loss', 'maximum_loss_over_pair_bound')},
        'log_contract_exact_not_numerically_evaluated': contract,
        'source_sha256': own_hashes,
        'document_sha256': document_hashes,
        'input_sha256': input_hashes,
        'portable_invocation': 'python check95.py <r095-directory> <validation.json> --project-root <repository-root>',
        'portable_source_resolution': 'r095 directory, work/r095, rounds/r095, code, public-code, public_code',
        'high_dimensional_source_construction_claimed': False,
        'finite_dimension_onset_claimed': False,
        'general_H_N_coefficient_claimed': False,
    }
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('Python_resources', 'Wolfram_resources')}, indent=2))


if __name__ == '__main__':
    main()
