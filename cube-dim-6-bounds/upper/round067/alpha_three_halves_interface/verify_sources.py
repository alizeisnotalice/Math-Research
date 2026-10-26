#!/usr/bin/env python3
"""Metadata only: hash old originals and copy exact JSON slots; never run old code."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

HERE = Path(__file__).resolve().parent
BASE = Path('/Users/zhengzhihao/Desktop/T')
ROOTS = {
    'C006': BASE/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round006/alternative',
    'F006_current': BASE/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/audit/round006/alpha15_review',
    'F006_historical': BASE/'Math-Research-cube6-literature/cube-dim-6-bounds/audit/round006/alpha15_review',
    'F006_original_path': BASE/'Math-Research-cube-dim-6/cube-dim-6-bounds/audit/round006/alpha15_review',
    'C007': BASE/'Math-Research-cube6-audit/cube-dim-6-bounds/upper/round007/alpha15_deterministic',
    'C008': BASE/'Math-Research-cube6-audit/cube-dim-6-bounds/upper/round008/alpha15_coercivity',
    'F008': BASE/'Math-Research-cube6-literature/cube-dim-6-bounds/audit/round008/alpha15_product',
    'A008': BASE/'Math-Research-cube6-theory/cube-dim-6-bounds/upper/round008/alpha15_global_g',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(name, obj):
    (HERE/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')


inventory = []
checks = []
for tag, root in ROOTS.items():
    for path in sorted(root.iterdir()):
        if path.is_file():
            inventory.append({'id': tag+'/'+path.name, 'path': str(path),
                              'sha256': sha(path), 'bytes': path.stat().st_size,
                              'execution_in_067': False})


def pin(owner, pointer, target, expected, alias=None):
    target = Path(target)
    actual = sha(target) if target.is_file() else None
    used = target
    status = 'PASS_EXACT_PATH'
    if actual != expected:
        assert alias is not None, (owner, pointer, str(target), expected, actual)
        used = Path(alias)
        assert sha(used) == expected, (owner, pointer, str(used))
        status = 'PASS_EXACT_CONTENT_AT_EXPLICIT_HISTORICAL_ALIAS'
    checks.append({'owner': str(owner), 'json_pointer': pointer,
                   'declared_target': str(target), 'expected_sha256': expected,
                   'current_declared_path_sha256': actual, 'resolved_target': str(used),
                   'resolved_sha256': sha(used), 'status': status})


# Producer/auditor file seals, including their full original JSON outputs.
for tag in ('C006', 'C007', 'C008', 'F008'):
    manifest = ROOTS[tag]/'manifest.json'
    obj = read(manifest)
    for name, entry in obj['files'].items():
        pin(manifest, '/files/'+name+'/sha256', ROOTS[tag]/name, entry['sha256'])
        assert (ROOTS[tag]/name).stat().st_size == entry['bytes']

# F006 has two byte-distinct scripts: only SRC changed from an absolute
# producer directory to a clone-relative directory; report bytes are identical.
for tag in ('F006_current', 'F006_historical', 'F006_original_path'):
    result = ROOTS[tag]/'check_results.json'
    obj = read(result)
    pin(result, '/script_sha256', ROOTS[tag]/'check.py', obj['script_sha256'])
    for name, expected in obj['input_hashes'].items():
        pin(result, '/input_hashes/'+name, ROOTS['C006']/name, expected)
old = (ROOTS['F006_historical']/'check.py').read_text()
new = (ROOTS['F006_current']/'check.py').read_text()
assert old.replace("SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round006/alternative')",
                   "SRC=HERE.parents[2]/'upper/round006/alternative'") == new
assert sha(ROOTS['F006_current']/'report.md') == sha(ROOTS['F006_historical']/'report.md')

for tag, field in (('C007','readonly_inputs'), ('C008','readonly_sources')):
    manifest = ROOTS[tag]/'manifest.json'
    for target, expected in read(manifest)[field].items():
        pin(manifest, '/'+field+'/'+target, target, expected,
            ROOTS['F006_historical']/Path(target).name if tag == 'C007' else None)

for tag, filename in (('F008','check_results.json'), ('A008','input_manifest.json'), ('A008','results.json')):
    source = ROOTS[tag]/filename
    obj = read(source)
    for target, expected in obj['input_sha256'].items():
        pin(source, '/input_sha256/'+target, target, expected)
    if tag == 'F008':
        pin(source, '/script_sha256', ROOTS[tag]/'check.py', obj['script_sha256'])
    else:
        pin(source, '/source_sha256', ROOTS[tag]/'check_bounds.py', obj['source_sha256'])

source = ROOTS['A008']/'checked_results.json'
obj = read(source)
pin(source, '/results_sha256', ROOTS['A008']/'results.json', obj['results_sha256'])
pin(source, '/checker_sha256', ROOTS['A008']/'verify_normalization.py', obj['checker_sha256'])
pdf = HERE/'sources/Bogdan-Grzywny-Ryznar_0905.2626v4.pdf'
pin('067 primary-source acquisition record', '/download_sha256', pdf,
    'fe711b64b9de875b38aba1c95f25b0ce0c435eab981f1d671c6c72db1ae3e15b')

for name in ('alternative_alpha/transfer.md', 'audit/014_alpha32_transfer_and_isotropic_audit.md'):
    path = BASE/'Math-Research-cube5-lower-audit-cont/cube-dim-5-bounds'/name
    inventory.append({'id': '5D_abstract_background_only/'+name, 'path': str(path),
                      'sha256': sha(path), 'bytes': path.stat().st_size,
                      'execution_in_067': False, 'numerical_results_imported': False})

slots = []


def copy_slots(tag, filename, pointers, meaning):
    path = ROOTS[tag]/filename
    obj = read(path)
    for pointer in pointers:
        value = obj
        for key in pointer.split('/')[1:]:
            value = value[int(key)] if isinstance(value, list) else value[key]
        slots.append({'source': str(path), 'source_sha256': sha(path),
                      'json_pointer': pointer, 'value_verbatim': value,
                      'meaning': meaning, 'recomputed_in_067': False})


copy_slots('F006_current', 'check_results.json', ['/status','/no_new_C6_bound'],
           'Previously checked analytic interface, conditional on quantitative g/I; not a new 067 proof.')
copy_slots('F008', 'check_results.json', ['/status','/strict_rational_bounds/G_norm','/strict_rational_bounds/I'],
           'Existing strict full-positive-face upper bounds for the same alpha=3/2 model.')
copy_slots('F008', 'check_results.json', ['/strict_rational_bounds/C007_deleted_parameter_tail',
           '/C007/ell','/C007/R','/C007/total_tail','/C007/interior_integral_certified'],
           'One fixed-point deleted-parameter-tail certificate only; interior expectation remains unpaid.')
copy_slots('A008','results.json',['/status','/best_strict_rational_global_g_lower','/global_g_proof_status',
           '/generic_C6_upper_retained','/budget_obstruction_scope'],
           'Existing full-cube g lower certificate, not an upper bound on true g.')
copy_slots('A008','checked_results.json',['/status','/best_global_g_strict_lower',
           '/actual_g_upper_bound_claimed','/actual_cost_lower_bound_claimed','/new_C6_bound',
           '/checked_schemes'],
           'Old four optimized certificate allocations paired with the fixed F008 norm envelope; not true optimal model cost.')

save('source_inventory.json', {
    'created_utc': datetime.now(timezone.utc).isoformat(),
    'status': 'PASS_SOURCE_CHAIN_WITH_EXPLICIT_C007_HISTORICAL_ALIAS',
    'kind': 'hash_and_JSON_read_only_no_scientific_computation',
    'actual_scientific_configurations': 0,
    'files': inventory, 'pin_checks': checks,
    'historical_metadata_not_revalidated_as_current': {
        'C006/manifest.json/readonly_main_inputs':
            'Historical mutable workflow snapshots; not used as current scientific premises.',
        'A008/executed_proof_input_sha256':
            'Historical pre-amendment proof hash retained verbatim in original files. Final theorem pin verified; no claim that historical bytes equal final bytes.'},
    'version_note': 'C007 old script/results pins resolve to exact literature-clone copies; current F006 path uses clone-relative SRC. Both script/result chains verified, and report is byte-identical.',
})
save('existing_results.json', {
    'status': 'VERBATIM_SAVED_RESULT_SLOTS_ONLY',
    'actual_scientific_configurations': 0,
    'display_decimals_parsed_or_reconstructed': False,
    'slots': slots,
})
print(json.dumps({'status':'PASS_METADATA_ONLY', 'source_files':len(inventory),
                  'source_pin_checks':len(checks), 'copied_exact_slots':len(slots),
                  'actual_scientific_configurations':0, 'old_programs_executed':0}))
