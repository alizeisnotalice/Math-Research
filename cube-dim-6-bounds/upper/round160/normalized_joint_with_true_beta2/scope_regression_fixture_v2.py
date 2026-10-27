#!/usr/bin/env python3
"""Engineering only: symbol binding and stub dispatch, no Fraction/model evaluation."""
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import argparse,hashlib,json,symtable,sys
sys.dont_write_bytecode=True
SOURCE_V2=Path(__file__).resolve().parent
def require_fixture_v2(ok,msg):
    if not ok:raise AssertionError(msg)
def load_fixture_module_v2(name,path):
    spec=spec_from_file_location(name,path);mod=module_from_spec(spec);spec.loader.exec_module(mod);return mod
def run_fixture_v2(receipt_path):
    oldtable=symtable.symtable((SOURCE_V2/'prepare_for_terra.py').read_text(),'v1_prepare','exec')
    oldphase=next(t for t in oldtable.get_children()if t.get_name()=='phaseA')
    require_fixture_v2(oldphase.lookup('importlib').is_local(),'ORIGINAL_IMPORTLIB_DEFECT_MUST_BE_REPRODUCIBLE_BY_BINDING')
    table=symtable.symtable((SOURCE_V2/'repair_prepare_v2.py').read_text(),'v2_prepare','exec')
    functions={t.get_name():t for t in table.get_children()}
    require_fixture_v2(functions['version_metadata_v2'].lookup('metadata_version_v2').is_global(),'METADATA_ALIAS_GLOBAL_BINDING')
    require_fixture_v2(functions['load_common_v2'].lookup('spec_from_file_location_v2').is_global()and functions['load_common_v2'].lookup('module_from_spec_v2').is_global(),'UTIL_ALIAS_GLOBAL_BINDING')
    repair=load_fixture_module_v2('c160_engineering_repair_fixture',SOURCE_V2/'repair_prepare_v2.py')
    called=[];repair.metadata_version_v2=lambda name:called.append(name)or('stub-'+name)
    require_fixture_v2(repair.version_metadata_v2()=={'numpy_version':'stub-numpy','scipy_version':'stub-scipy'}and called==['numpy','scipy'],'METADATA_ALIAS_ACTUAL_STUB_DISPATCH')
    loader=SimpleNamespace(exec_module=lambda m:setattr(m,'engineering_loaded',True))
    repair.spec_from_file_location_v2=lambda name,path:SimpleNamespace(loader=loader)
    repair.module_from_spec_v2=lambda spec:SimpleNamespace()
    require_fixture_v2(repair.load_common_v2(Path('/ENGINEERING_STUB_ONLY')).engineering_loaded is True,'UTIL_ALIAS_ACTUAL_STUB_DISPATCH')
    operator=load_fixture_module_v2('c160_engineering_operator_fixture',SOURCE_V2/'operator_prepare_v2.py')
    rejected=False
    try:operator.validate_manifest_v2(SOURCE_V2,{'files':[{'path':'common.py','bytes':1,'sha256':'invalid'}]})
    except ValueError as err:rejected=str(err)=='PACKAGE_ENTRY_REQUIRES_LOCAL_NAME_NOT_PATH'
    require_fixture_v2(rejected,'WRONG_PATH_SCHEMA_FAILS_BEFORE_FILE_OR_CHILD_ACCESS')
    subprocess_calls=[];operator.subprocess=SimpleNamespace(run=lambda *a,**k:subprocess_calls.append((a,k)))
    with TemporaryDirectory(prefix='c160_scope_engineering_')as temp:
        temp=Path(temp);(temp/'PACKAGE_MANIFEST.json').write_text(json.dumps({'files':[{'path':'common.py','bytes':1,'sha256':'invalid'}]}))
        operator.SOURCE_V2=temp;oldargv=sys.argv;sys.argv=['operator_fixture','--root',str(temp/'unused_target')]
        stopped=False
        try:operator.main_v2()
        except ValueError as err:stopped=str(err)=='PACKAGE_ENTRY_REQUIRES_LOCAL_NAME_NOT_PATH'
        finally:sys.argv=oldargv
    require_fixture_v2(stopped and subprocess_calls==[],'ACTUAL_OPERATOR_MAIN_FAILS_BEFORE_FIXTURE_OR_PREPARE_CHILD')
    result={'status':'ENGINEERING_SCOPE_FIXTURE_PASS','original_importlib_local_shadow_observed':True,'new_metadata_and_util_alias_bindings_global':True,'actual_stub_alias_dispatch_passed':True,'wrong_package_path_schema_rejected':True,'actual_operator_main_fail_fast_before_children':True,'mathematical_functions_called':0,'prepare_or_science_children':0,'repair_prepare_v2_sha256':hashlib.sha256((SOURCE_V2/'repair_prepare_v2.py').read_bytes()).hexdigest(),'operator_wrapper_v2_sha256':hashlib.sha256((SOURCE_V2/'operator_prepare_v2.py').read_bytes()).hexdigest()}
    with(Path(receipt_path)).open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--receipt',required=True);args=ap.parse_args();run_fixture_v2(args.receipt)
