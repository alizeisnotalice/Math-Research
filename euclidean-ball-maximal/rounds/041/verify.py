"""Run the final round-41 probes without retaining regenerable ledgers."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    commands = []
    with tempfile.TemporaryDirectory(prefix='euclidean-round41-') as folder:
        folder = Path(folder)
        records = {}
        for script in ('square_probe.py', 'gap_probe.py'):
            output = folder / (script + '.json')
            run = subprocess.run([sys.executable, str(here/script), '--output', str(output)],
                                 env=env, text=True, capture_output=True, check=True)
            commands.append({'script': script, 'exit_code': run.returncode})
            records[script] = json.loads(output.read_text())
        square = records['square_probe.py']
        gap = records['gap_probe.py']
        assert square['status'] == gap['status'] == 'passed'
        assert square['new_case_count'] == 19 and len(square['cases']) == 21
        assert square['independent_maximum_verification']['status'] == 'passed'
        assert [v['k'] for v in gap['cases']] == [10, 12, 14]
        assert [v['I4k'] for v in gap['cases']] == ['495/262144', '1971/1048576', '15759/8388608']
        best = square['maximum_Q_over_X']
        files = [here/'square_probe.py', here/'gap_probe.py', here/'verify.py']
        files += sorted((root/'runtime/work/general_input_20261003').glob('*.py'))
        receipt = {
            'status': 'passed', 'commands': commands,
            'arithmetic': 'Exact rational arithmetic; decimal fields are displays only.',
            'square': {
                'new_batches': 3, 'new_cases': 19, 'old_benchmarks': 2,
                'maximum_label': best['label'],
                'maximum': {k: best[k] for k in ('M', 'Q', 'X', 'Q_over_X', 'maxA')},
                'independent_maximum_verification': square['independent_maximum_verification'],
            },
            'gap': [{k:v for k,v in row.items() if k != 'source'} for row in gap['cases']],
            'source_hashes': {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'scope': 'One-dimensional finite input checks; arbitrary-depth no-decay result is analytic in report.md.',
            'main_dimension_independent_weak_bound_proved': False,
            'all_regenerable_full_ledgers_discarded_after_verification': True,
        }
        print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
