"""Only frozen JSON/integer geometry is loaded; no H/log arithmetic module."""
import sys,json,hashlib,importlib.util
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
REG=json.loads((HERE/'registration.json').read_text())
path=REG['accepted081_common']
assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==REG['accepted081_common_sha256']
spec=importlib.util.spec_from_file_location('reuse122_common',path)
common=importlib.util.module_from_spec(spec);spec.loader.exec_module(common)
