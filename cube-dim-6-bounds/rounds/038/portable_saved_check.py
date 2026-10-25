"""Check identical copied bytes then execute original path-anchored pure verifier.

No writer, native call, kernel, additional candidate, or changes to frozen code.
"""
import sys
sys.dont_write_bytecode=True
import json,hashlib,importlib.util
from pathlib import Path
P=Path(__file__).resolve().parents[2]
D=P/'upper/round038/green_floor_diagnostic'
O=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round038/green_floor_diagnostic')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(O/'manifest.json')=='4f8e79acbe2a7a36c5bc6ea9eed59c3b140d856d8f7244fd1801071137675e54'
m=json.loads((O/'manifest.json').read_text())
for rel,digest in m['whitelist_sha256'].items():assert sha(O/rel)==sha(D/rel)==digest,rel
assert sha(O/'manifest.json')==sha(D/'manifest.json')
s=importlib.util.spec_from_file_location('immutable_A038_saved_checker',D/'check_saved.py')
z=importlib.util.module_from_spec(s);s.loader.exec_module(z)
# Candidate rows intentionally store absolute previous-run paths. Original
# HERE preserves their declared meaning after identical-copy hash checks.
z.HERE=O
result=z.check()
print(json.dumps({'portable_copy_byte_equivalence':True,'immutable_original_path_semantics_preserved':True,'original_verifier_sha256':sha(O/'check_saved.py'),'actual_result':result,'new_kernel_calls':0},ensure_ascii=False),flush=True)
