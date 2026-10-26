"""Metadata-only original/archive SHA checks and verbatim saved-slot extraction."""
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import time

resource.setrlimit(resource.RLIMIT_CPU, (30, 30))
signal.alarm(170)
start = time.monotonic()
out = Path(__file__).resolve().parent
top = Path('/Users/zhengzhihao/Desktop/T')
archive = top / 'Math-Research-cube-dim-6/cube-dim-6-bounds'
specs = [
    ('C012', 'Math-Research-cube6-curvature-astra', 'upper/round012/augmented_frame', 'delivery_receipt.json'),
    ('A012', 'Math-Research-cube6-theory', 'upper/round012/closest_subspaces', None),
    ('F012', 'Math-Research-cube6-literature', 'audit/round012/augmented_frame_theory', 'manifest.json'),
    ('A015', 'Math-Research-cube6-theory', 'upper/round015/band_taylor_cover', None),
    ('F015', 'Math-Research-cube6-literature', 'audit/round015/band_taylor_review', 'manifest.json'),
    ('B018', 'Math-Research-cube6-upper-astra', 'upper/round018/augmented_obstruction_probe', 'manifest.json'),
]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(name, data):
    (out / name).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')

records, checks, gaps, roots = [], [], [], {}
for tag, clone, rel, seal in specs:
    original = top / clone / 'cube-dim-6-bounds' / rel
    roots[tag] = original
    if not original.is_dir():
        gaps.append({'original_directory_missing': str(original)})
        continue
    for p in sorted(original.iterdir()):
        if not p.is_file() or p.name.startswith('.'):
            continue
        archived = archive / rel / p.name
        rec = {'tag': tag, 'path': str(p), 'bytes': p.stat().st_size,
               'sha256': sha(p), 'archive_path': str(archived),
               'archive_exists': archived.is_file()}
        rec['archive_sha256'] = sha(archived) if archived.is_file() else None
        rec['archive_identical'] = rec['sha256'] == rec['archive_sha256']
        records.append(rec)
        if not rec['archive_identical']:
            gaps.append({'archive_difference': rec})
    if seal:
        doc = json.loads((original / seal).read_text())
        files = doc['files']
        items = files if isinstance(files, list) else [dict(path=k, **v) for k, v in files.items()]
        for item in items:
            p = original / item['path']
            actual = sha(p) if p.is_file() else None
            rec = {'seal_path': str(original / seal), 'file_path': str(p),
                   'expected_sha256': item['sha256'], 'actual_sha256': actual,
                   'match': actual == item['sha256']}
            checks.append(rec)
            if not rec['match']:
                gaps.append({'sealed_file_mismatch': rec})

slots = []
def copy_slots(tag, name, keys):
    p = roots[tag] / name
    data = json.loads(p.read_text())
    for k in keys:
        slots.append({'source_path': str(p), 'source_sha256': sha(p),
                      'json_pointer': '/' + k, 'verbatim_value': data[k]})

if not gaps:
    copy_slots('A012', 'results.json', [
        'status', 'body_directions', 'augmented_directions', 'beta_test',
        'pairs_per_point', 'augmented_pairs_per_point', 'geometrically_distinct_planes',
        'augmented_distinct_planes_total', 'body_axis_planes_add_to_S2',
        'all_beta_positive_span_proof', 'tubes_disjoint_claimed'])
    copy_slots('C012', 'fan_audit.json', [
        'status', 'beta', 'T', 'Gram', 'base_cones', 'full_signed_permutation_group',
        'complete_group_terms', 'ray_occurrences_checked', 'unique_integer_rays',
        'all_coordinate_forms_nonnegative_on_each_base_cone',
        'all_38_linear_forms_constant_weak_sign_per_cone',
        'rays_determinants_unchanged', 'new_fan_sha256'])
    copy_slots('C012', 'checked_results.json', [
        'status', 'rows', 'cumulative_augmented_full_passes', 'global_m_certified'])
    copy_slots('B018', 'checked_results.json', [
        'status', 'beta', 'T', 'frame_sha256', 'numerators', 'denominator',
        'K_interval_exact', 'coarse_rational_K_enclosure',
        'm64_minus_K_upper_exact', 'rejects_strict64', 'rejects625',
        'true_augmented_S2', 'native_confluent_terms', 'native_zero_terms',
        'scope', 'new_C6_bound', 'global_minimum_certified'])
    copy_slots('B018', 'manifest.json', ['native_passes', 'postprocess_native_calls', 'active_handles'])
write('source_binding.json', {'status': 'PASS' if not gaps else 'GAP',
      'scope': 'metadata-only; no scientific replay; original and main archive compared',
      'source_files': records, 'internal_seal_checks': checks, 'gaps': gaps})
write('saved_slots.json', {'mode': 'verbatim_saved_json_values_not_recomputed', 'slots': slots})
receipt = {'status': 'PASS_METADATA_ONLY' if not gaps else 'PROVENANCE_GAP',
           'pid': os.getpid(), 'elapsed_seconds': time.monotonic() - start,
           'cpu_seconds': time.process_time(), 'source_files': len(records),
           'internal_seal_checks': len(checks), 'saved_slots': len(slots),
           'new_scientific_configurations': 0, 'new_kernel_evaluations': 0,
           'new_beta_candidates': 0, 'old_author_programs_executed': 0,
           'active_child_processes': [], 'gaps': len(gaps)}
write('metadata_receipt.json', receipt)
print(json.dumps(receipt, ensure_ascii=False))
if gaps:
    raise SystemExit(1)
