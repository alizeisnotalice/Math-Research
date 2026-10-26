"""Hash-verified, read-only accepted integer arithmetic and geometry."""
import sys, json, hashlib, importlib.util
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
REG=json.loads((HERE/'registration.json').read_text())


def verified(label,path,expected):
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


common=verified('mixed100_readonly_common081',REG['accepted081_common'],REG['accepted081_common_sha256'])
arith=verified('mixed100_readonly_arith081',REG['accepted081_arith'],REG['accepted081_arith_sha256'])
