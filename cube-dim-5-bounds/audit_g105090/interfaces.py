"""Read-only binding to the accepted081 integer/tree primitives."""
import sys,hashlib,importlib.util,json
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
REG=json.loads((HERE/'registration.json').read_text())


def load_verified(label,path,expected):
    p=Path(path)
    assert hashlib.sha256(p.read_bytes()).hexdigest()==expected,(label,'accepted source changed')
    spec=importlib.util.spec_from_file_location(label,p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


common=load_verified('accepted_common081',REG['accepted081_common'],REG['accepted081_common_sha256'])
arith=load_verified('accepted_arith081',REG['accepted081_arith'],REG['accepted081_arith_sha256'])
