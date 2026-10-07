#!/usr/bin/env python3
"""Regenerate round 42 in a temporary directory; print a lean verification receipt."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from fractions import Fraction as F
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]

def main():
    with tempfile.TemporaryDirectory(prefix='euclidean42_') as folder:
        tmp = Path(folder)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1')
        def run(name, *extra):
            output = tmp / (name + '.json')
            subprocess.run([sys.executable, str(HERE/(name+'.py')), '--output', str(output), *extra],
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
            result = json.loads(output.read_text())
            assert result['status'] == 'passed'
            return result
        row = run('row_probe')
        assert row['new_case_count'] == 44 and len(row['cases']) == 46
        assert row['independent_maximum_case_verification']['status'] == 'passed'
        assert F(row['maximum_row']['Rj_over_Xj']) == F(49365547633297546998997619,308583131727672796791177216)
        flux = run('flux_probe', '--row-input', str(tmp/'row_probe.json'))
        assert len(flux['cases']) == 4
        micro = [dict(label=c['label'], **next(r for r in c['rows'] if r['j']==4)) for c in flux['cases'][:3]]
        assert F(micro[1]['positive_flow'])-F(micro[1]['negative_flow']) == F(81,268435456)
        geometry = run('geometry_probe')
        assert [F(v['mass_over_lens_length']) for v in geometry['lens']] == [F(63,2),F(255,2),F(1023,2)]
        high = run('highdim_probe')
        assert high['configuration_count'] == 25 and high['replicate_count'] == 50
        candidates = [dict(family=c['family'], **r) for c in high['cases'] for r in c['replicates'] if r['best_row_ratio_estimate'] is not None]
        best = max(candidates, key=lambda r:r['best_row_ratio_estimate'])
        dependencies = sorted((PROJECT/'runtime').rglob('*.py')) + [HERE.parent/'041'/'square_probe.py'] + sorted(HERE.glob('*.py'))
        receipt = dict(status='passed', round=42, date='2026-10-07',
            python=sys.version.split()[0], numpy=np.__version__,
            exact=dict(new_cases=44, reference_cases=2, batches=3,
                       maximum_row_label=row['maximum_row_label'], maximum_row=row['maximum_row'],
                       independent_verification=row['independent_maximum_case_verification'],
                       flux_cases=4, flux_rows=sum(len(c['rows']) for c in flux['cases']), micro_j4=micro,
                       shell_dimensions=[v['n'] for v in geometry['shell']], lens=geometry['lens']),
            exploratory=dict(configurations=25, replicates=50, samples_per_replicate=32768,
                seeds='420711+100*case_index+replicate_index', max_dimension=2048,
                analytic_baselines=5, best_observed=dict(family=best['family'], n=best['n'], seed=best['seed'],
                    j=best['best_observed_row'], ratio_estimate=best['best_row_ratio_estimate']),
                limitations=high['warning']),
            main_bound_proved=False, row_counterexample_found=False,
            sha256={str(f.relative_to(PROJECT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in dependencies})
        print(json.dumps(receipt, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
