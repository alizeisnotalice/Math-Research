#!/usr/bin/env python3
"""V2 gate only; delegate to unchanged original phaseA runner and original output paths."""
from pathlib import Path
from importlib.util import spec_from_file_location as spec_v2
from importlib.util import module_from_spec as module_v2
import sys
sys.dont_write_bytecode=True
from repair_prepare_v2 import read_v2,sha_v2,require_v2,assert_snapshot_v2,assert_v2_files_v2,contract_v2
ROOT_V2=Path(__file__).resolve().parent
def run_v2():
    require_v2(not sys.argv[1:],'usage: run_phaseA_v2.py (after NEW exact-hash dual gate only)')
    assert_v2_files_v2(ROOT_V2);assert_snapshot_v2(ROOT_V2,ROOT_V2)
    receipt=read_v2(ROOT_V2/'prepare_v2/exit_receipt.json')
    require_v2(receipt['status']=='V2_METADATA_RECOVERY_COMPLETE_SCIENCE_UNRUN','COMPLETED_V2_METADATA_ONLY_RECEIPT')
    for key,path in [('registration_sha256',ROOT_V2/'phaseA/registration.json'),('release_requirements_sha256',ROOT_V2/'phaseA/release_requirements.json'),('scope_fixture_receipt_sha256',ROOT_V2/'operator_v2_scope_fixture_receipt.json')]:require_v2(receipt[key]==sha_v2(path),'V2_RECEIPT_BINDS_OUTPUT:'+key)
    requirements=contract_v2(ROOT_V2);require_v2(read_v2(ROOT_V2/'phaseA/release_requirements.json')==requirements,'SAME_V2_RELEASE_CONTRACT')
    release=read_v2(ROOT_V2/'phaseA/parent_release.json')
    for k,v in requirements.items():require_v2(release.get(k)==v,'NEW_V2_RELEASE:'+k)
    require_v2(release.get('main_approved')is True and release.get('independent_review_approved')is True,'NEW_EXACT_HASH_DUAL_GATE_REQUIRED')
    # The original runner performs its unchanged source156 saved acceptance,
    # exact-model checks, single launch marker and <=150-call numerical selector.
    spec=spec_v2('c160_preserved_phaseA_v2_entry',ROOT_V2/'run_phaseA.py');mod=module_v2(spec);spec.loader.exec_module(mod);mod.run()
if __name__=='__main__':run_v2()
