"""One finite pure-field equivalence and byte archive; no numerical helpers.

This imports only the field/JSON portion of portable180. It never imports the
frozen mathematical reader or invokes any Fraction/exp/q/corner/forest logic.
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import time
from portable180 import aggregate, digest, encoded, load, original_leaf_view, original_source_view, repository_parameters, sha

HERE = Path(__file__).resolve().parent
ORIGINAL = HERE.parent


def write_new(path, value):
    with Path(path).open('xb') as handle:
        handle.write(encoded(value)+b'\n')


def main(plan_sha):
    start = time.monotonic()
    planpath = HERE/'archive_plan.json'
    assert sha(planpath) == plan_sha
    plan = load(planpath)
    assert not (HERE/'archive_registration.json').exists()
    assert not (HERE/'evidence').exists()
    for name, expected in plan['code_locks'].items():
        assert sha(HERE/name) == expected
    for relative, expected in plan['selected_original_files'].items():
        p = ORIGINAL/relative
        assert sha(p) == expected['sha256'] and p.stat().st_size == expected['bytes']
    original_plan = load(ORIGINAL/'plan.json')
    for spec in original_plan['specs']:
        expected = original_plan['sources'][spec['path']]
        assert sha(spec['path']) == expected['sha256'] and Path(spec['path']).stat().st_size == expected['bytes']
    for relative, expected in plan['repository_locks'].items():
        p = Path(plan['repo_base'])/relative
        assert sha(p) == expected['sha256'] and p.stat().st_size == expected['bytes']
    write_new(HERE/'archive_registration.json',{
        'pid':os.getpid(),'unix':time.time(),'archive_plan_sha256':plan_sha,
        'source_contract':107530,'support_contract':3812,'leaf_contract':111201,
        'operation':'pure JSON field equality, SHA, byte copy/readback only',
        'new_mathematical_cases':0,'raw_state_replay':False,'wallcap':None,
        'automatic_retry':False})
    counts = {'sources':0,'supports':0,'leaves':0}
    try:
        sources,leaves,supports,common = repository_parameters(plan['repo_base'])
        seen_sources,seen_leaves,seen_supports = set(),set(),set()
        for spec in original_plan['specs']:
            with gzip.open(spec['path'],'rt') as handle:
                original = json.load(handle)
            for task in original['sources']:
                cid = task['source155']['canonical_id']
                assert cid not in seen_sources
                assert original_source_view(task,common['upper']) == sources[cid], ('source_field_mismatch',cid)
                seen_sources.add(cid)
                counts['sources'] += 1
            for task in original['leaves']:
                idx = task['leaf151']['leaf_index']
                assert idx not in seen_leaves
                assert original_leaf_view(task) == leaves[idx], ('leaf_field_mismatch',idx)
                seen_leaves.add(idx)
                counts['leaves'] += 1
            for task in original['supports']:
                node = tuple(task['support']['node'])
                assert node not in seen_supports
                assert task == supports[node], ('support_field_mismatch',node)
                seen_supports.add(node)
                counts['supports'] += 1
            write_new(HERE/f'equivalence-shard{spec["worker"]:02}.json',{
                'original_input':spec['path'],'sha256':sha(spec['path']),
                'source_cases':len(original['sources']),'leaf_cases':len(original['leaves']),
                'support_cases':len(original['supports']),'status':'exact_field_equality',
                'new_mathematical_cases':0})
        assert seen_sources == set(sources) and seen_leaves == set(leaves) and seen_supports == set(supports)
        assert counts == {'sources':107530,'supports':3812,'leaves':111201}
        result = {
            'status':'all_parameter_views_exactly_equal_no_new_mathematics',
            'pid':os.getpid(),'counts':counts,
            'parameter_view_sha256':{'sources':aggregate(sources),'leaves':aggregate(leaves),'supports':aggregate(supports)},
            'all_original_canonical_leaf_support_sets_exact':True,
            'Delta115_policy':'input view is original Delta110 + same common_log_U; original subtraction identity inherits actual31006, not re-evaluated',
            'unrecoverable_mathematical_fields':[], 'numerical_supplement_files':[],
            'new_EXP_q_division_corners_forest_or_H_cases':0,
            'raw_state_binding':'inherited frozen actual31006, 109446 raw bindings; no 9-state replay',
            'seconds':time.monotonic()-start}
        write_new(HERE/'equivalence_result.json',result)
        write_new(HERE/'repository_locks.json',plan['repository_locks'])
        for relative, expected in plan['selected_original_files'].items():
            destination = HERE/'evidence'/relative
            destination.parent.mkdir(parents=True,exist_ok=True)
            assert not destination.exists()
            shutil.copyfile(ORIGINAL/relative,destination)
            assert sha(destination) == expected['sha256'] and destination.stat().st_size == expected['bytes']
        oldmanifest = load(ORIGINAL/'manifest.json')
        omitted = {}
        for relative, record in oldmanifest['retained_files'].items():
            if relative in plan['selected_original_files']:
                continue
            assert relative.endswith('/chunks.jsonl')
            sid = str(int(Path(relative).parent.name.removeprefix('worker')))
            omitted[relative] = {**record,
                'last_verified_chain':oldmanifest['worker_journal_integrity'][sid]['chain'],
                'chunks':oldmanifest['worker_journal_integrity'][sid]['chunks'],
                'policy':'retained at original local path; whole-byte SHA/length and completed chain inherited; not copied'}
        assert len(omitted) == 15
        write_new(HERE/'omitted_local_evidence.json',{
            'original_manifest_sha256':sha(ORIGINAL/'manifest.json'),
            'original_input_packages':{s['path']:original_plan['sources'][s['path']] for s in original_plan['specs']},
            'journals':omitted,'old_raw_state_locks':'preserved verbatim in evidence/plan.json',
            'mathematical_parameters':'all recoverable from repo151/155; exact views proved by equivalence_result.json'})
        # Read back every selected source again after copying; no math evaluation.
        for relative, expected in plan['selected_original_files'].items():
            assert sha(ORIGINAL/relative) == expected['sha256']
            assert sha(HERE/'evidence'/relative) == expected['sha256']
        for spec in original_plan['specs']:
            assert sha(spec['path']) == original_plan['sources'][spec['path']]['sha256']
        names = sorted(p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file())
        manifest = {'status':'portable_archive_complete_without_new_mathematics',
            'actual_archive_pid':os.getpid(),'original_actual_math_pid':31006,
            'archive_seconds':time.monotonic()-start,'source_files_byte_unchanged':True,
            'selected_original_files':len(plan['selected_original_files']),
            'all_selected_bytes_readback_equal':True,
            'files':{name:{'sha256':sha(HERE/name),'bytes':(HERE/name).stat().st_size} for name in names},
            'copied_H_gradient_or_numerical_datasets':False,
            'new_mathematical_cases':0,'verify_mathematics_entry_executed':False,
            'shared_original_H_enclosures':True,'commit_or_push':False}
        write_new(HERE/'portable_manifest.json',manifest)
        print(json.dumps({'status':manifest['status'],'pid':os.getpid(),
            'seconds':manifest['archive_seconds'],'equivalence_result_sha256':sha(HERE/'equivalence_result.json'),
            'portable_manifest_sha256':sha(HERE/'portable_manifest.json'),
            'new_mathematical_cases':0}),flush=True)
    except BaseException as error:
        import traceback
        write_new(HERE/'archive_failure.json',{'pid':os.getpid(),'exception':repr(error),
            'traceback':traceback.format_exc(),'completed_counts':counts,
            'automatic_retry':False,'new_mathematical_cases':0})
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reviewed-plan-sha',required=True)
    main(parser.parse_args().reviewed_plan_sha)
