"""Portable audit entry: validate/import by default; compute only with --run.

Usage: python -B portable.py --project-root /path/to/repository
       python -B portable.py --project-root /path/to/cube-dim-5-bounds --run
The four frozen audit artifacts are never modified. --run writes to a new
temporary directory and prints its path; the historical author path is remapped.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile

HERE = Path(__file__).resolve().parent
FROZEN = {
    'independent_check.py': 'd6d395c74fefc87174effc2bfb4f50ac7a9ffe96a44e2027d9b2f23a7a7b3f5e',
    'independent_check.json': 'a57dc3bea1a1d54f1482acca338ec763862976fee0f48eae4aa56d31b901cdff',
    'registration.md': 'ca81eb48b4dad0109736f5cdf4ab42308cfb424ee693c567517af46a19be1ccd',
    'report.md': 'b3fca5107ccf51c5cead9f235841dadbf4dd3c8d441b2978e8de3d24731f1b1c',
}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', required=True, type=Path,
                        help='repository root or cube-dim-5-bounds directory containing astra_upper001/002')
    parser.add_argument('--run', action='store_true', help='explicitly run the finite audit after validation')
    args = parser.parse_args()
    project = args.project_root.resolve()
    if (project/'cube-dim-5-bounds').is_dir():
        project = project/'cube-dim-5-bounds'
    for name, expected in FROZEN.items():
        if digest(HERE/name) != expected:
            raise ValueError('Frozen audit artifact hash mismatch: '+name)
    receipt = json.loads((HERE/'independent_check.json').read_text())
    path_map = {}
    for legacy, expected in receipt['author_input_sha256'].items():
        parts = Path(legacy).parts[-2:]
        if parts[0] not in ('astra_upper001','astra_upper002') or parts[1] not in (
                'registration.md','report.md','verify.py','verification.json'):
            raise ValueError('Unexpected locked source name: '+legacy)
        local = project.joinpath(*parts)
        if digest(local) != expected:
            raise ValueError('Locked author input hash mismatch: '+str(local))
        path_map[legacy] = local
    spec = importlib.util.spec_from_file_location('independent_g_audit', HERE/'independent_check.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # Definitions only; its __main__ guard remains false.
    module.AUTHOR = project
    # The frozen main also checks historical JSON hash keys through Path(path).
    # Map only those eight locked names, leaving all other Path calls intact.
    module.Path = lambda value: path_map.get(str(value), Path(value))
    if args.run:
        output = Path(tempfile.mkdtemp(prefix='audit-g011-run-'))
        module.ROOT = output
        module.HERE = output
        shutil.copy2(HERE/'registration.md', output/'registration.md')
        module.main()
        print(json.dumps({'mode':'run','output_directory':str(output)}))
    else:
        with tempfile.TemporaryDirectory(prefix='audit-g011-import-') as output:
            module.ROOT = Path(output)
            module.HERE = Path(output)
            print(json.dumps({'mode':'validate_and_import_only','calculation_executed':False,
                              'frozen_artifacts_verified':len(FROZEN),
                              'locked_author_inputs_verified':len(path_map),
                              'project_root':str(project)}))

if __name__ == '__main__':
    main()
