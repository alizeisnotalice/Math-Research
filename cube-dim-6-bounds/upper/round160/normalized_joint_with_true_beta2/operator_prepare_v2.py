#!/usr/bin/env python3
"""Fail-fast engineering wrapper: fixed name schema, fixture, then metadata recovery."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
sys.dont_write_bytecode=True
SOURCE_V2=Path(__file__).resolve().parent
def validate_manifest_v2(root,manifest):
    # Package manifests use name. External source manifests use path separately.
    if not isinstance(manifest.get('files'),list)or not manifest['files']:raise ValueError('PACKAGE_FILES_ARRAY_REQUIRED')
    for row in manifest['files']:
        if set(row)!={'name','bytes','sha256'}or not isinstance(row['name'],str)or Path(row['name']).name!=row['name']:raise ValueError('PACKAGE_ENTRY_REQUIRES_LOCAL_NAME_NOT_PATH')
        path=Path(root)/row['name']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']or path.stat().st_size!=row['bytes']:raise ValueError('PACKAGE_HASH_OR_SIZE:'+row['name'])
def main_v2():
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);args=ap.parse_args();root=Path(args.root).resolve()
    # Any validation exception exits before any child prepare is invoked.
    for name in ['PACKAGE_MANIFEST.json','V2_PACKAGE_MANIFEST.json']:validate_manifest_v2(SOURCE_V2,json.loads((SOURCE_V2/name).read_text()))
    obs=json.loads((SOURCE_V2/'PARTIAL_FAILURE_OBSERVATION.json').read_text())
    if str(root)!=obs['target_root']:raise ValueError('ONLY_FROZEN_PARTIAL_ROOT')
    if(root/'prepare_v2').exists()or(root/'operator_v2_scope_fixture_receipt.json').exists():raise ValueError('NO_AUTOMATIC_OPERATOR_RETRY')
    # check=True is the failure barrier; no shell and no command chaining.
    subprocess.run([sys.executable,'-B',str(SOURCE_V2/'scope_regression_fixture_v2.py'),'--receipt',str(root/'operator_v2_scope_fixture_receipt.json')],check=True)
    subprocess.run([sys.executable,'-B',str(SOURCE_V2/'repair_prepare_v2.py'),'--root',str(root)],check=True)
if __name__=='__main__':main_v2()
