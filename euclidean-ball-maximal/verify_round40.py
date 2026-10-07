"""Reproduce round 40 in disposable storage, without altering archived inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regenerate', action='store_true', help='Compatibility flag; generation is always performed.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    archive = root / 'rounds/040'
    manifest = json.loads((archive / 'provenance.json').read_text())
    with tempfile.TemporaryDirectory(prefix='euclidean-round40-') as scratch:
        scratch = Path(scratch)
        work = scratch / 'work/general_input_20261003'
        output = scratch / 'output/general_input_20261003'
        work.mkdir(parents=True)
        output.mkdir(parents=True)
        for source in sorted((root / 'runtime/work/general_input_20261003').glob('*.py')):
            assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest[source.name]['archived_sha256']
            shutil.copyfile(source, work / source.name)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        commands = []
        def run(name, *options):
            result = subprocess.run([sys.executable, str(work/name), *options], env=env,
                                    text=True, capture_output=True, check=True)
            commands.append({'script': name, 'arguments': list(options), 'exit_code': result.returncode})
        run('stage40_band_search.py')
        run('stage40_band_search.py', '--verify-saved')
        run('stage40_independent_band_counterexample.py')
        data = json.loads((output/'stage40_band_search.json').read_text())
        independent = json.loads((output/'stage40_independent_band_counterexample.json').read_text())
        assert data['verification']['all_saved_summary_count'] == 28
        assert data['verification']['complete_scope_count'] == 7
        assert independent['E1'] == '14879934083/1125899906842624'
        assert independent['bridge_support_pieces'] == 110
        stats = data['statistics']
        assert stats['certified_E1_positive_count'] == 11
        assert stats['certified_E1_zero_count'] == 17
        assert stats['threshold_reached_count'] == 0
        independent.pop('independent_band_cells', None)
        receipt = {'status': 'passed', 'mode': 'regenerate',
                   'commands': commands, 'verification': data['verification'],
                   'counts': {k:v for k,v in stats.items() if k.endswith('count')},
                   'independent': independent, 'main_theorem_proved': False}
        print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
