"""Replay frozen code in a new directory; never overwrites archived evidence."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def prepare(target):
    target = Path(target).resolve()
    if target.exists():
        raise RuntimeError('Output directory must not exist')
    target.mkdir(parents=True)
    for name in ['round038_highdim', 'round038_conditional', 'round038_pb_bound']:
        shutil.copytree(HERE/name, target/name)
    # Preserve historical source bytes and their registered hashes; map paths
    # in an isolated replay tree rather than editing any frozen program.
    for source, folder in [(HERE.parent/'037/cell/probe.py', 'round037_quantile_cell'),
                           (HERE.parent/'034/block/block_probe.py', 'round034_block_mixture')]:
        (target/folder).mkdir()
        shutil.copy2(source, target/folder/source.name)
    return target


def execute(root, relative):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    subprocess.run([sys.executable, str(root/relative)], check=True, env=env,
                   cwd=root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['highdim', 'conditional', 'pb', 'source'])
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    root = prepare(args.output)
    if args.mode == 'highdim':
        for f in ['HIGH_DIM_PLAN.json', 'HIGH_DIM_REGISTRATION.json', 'HIGH_DIM_RESULTS.json']:
            (root/'round038_highdim'/f).unlink()
        execute(root, 'round038_highdim/probe_highdim.py')
    elif args.mode == 'conditional':
        for f in ['RESULTS.json', 'LAYERED_RESULTS.json']:
            (root/'round038_conditional'/f).unlink()
        execute(root, 'round038_conditional/probe.py')
        execute(root, 'round038_conditional/layered_postprocess.py')
    elif args.mode == 'pb':
        for f in ['registration.json', 'results.json', 'checks.json']:
            (root/'round038_pb_bound'/f).unlink()
        execute(root, 'round038_pb_bound/probe.py')
        execute(root, 'round038_pb_bound/check_bounds.py')
    else:
        (root/'round038_highdim/SOURCE_VALIDATION_normalized_reference.json').unlink()
        execute(root, 'round038_highdim/validate_source_tables.py')


if __name__ == '__main__':
    main()
