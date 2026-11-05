"""Portable adaptation of the R076 final version/eroder guard.

Original reviewer script SHA256:
cc59a57aae645eb024e736c6332a62d2b8d8e4b903661d124f2adf9beac0de43
Scope of adaptation: read the final public proof and final-review-manifest only;
reconstruct old and intermediate guard bytes in memory, then preserve the
original version, exact patch, clip inequalities and counterexample assertions.
No source text is written or edited.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import difflib, hashlib, json
here=Path(__file__).resolve().parent
root=here.parents[1]
final=(root/"proofs/076-finite-depth-supply.md").read_bytes()
manifest=json.loads((here/"final-review-manifest.json").read_text())
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(final)==manifest["final_sha256"]
ops=manifest["operations"]
# Recreate the pre-guard source by reversing all declared hunks, in memory.
old_lines=final.decode().splitlines(keepends=True)
for z in reversed(ops):
    a=z["new_start_line0"]; b=z["new_end_line0"]
    assert "".join(old_lines[a:b])==z["new_text"]
    old_lines[a:b]=z["old_text"].splitlines(keepends=True)
old="".join(old_lines).encode()
assert sha(old)==manifest["old_sha256"]
# The intermediate eroder-guard file is final with only the later mathcal typo fix reversed.
fix_ops=[z for z in ops if "mathcal" in z["old_text"] and "mathcal" in z["new_text"]]
assert len(fix_ops)==1
z=fix_ops[0]
guard_lines=final.decode().splitlines(keepends=True)
a=z["new_start_line0"];b=z["new_end_line0"]
assert "".join(guard_lines[a:b])==z["new_text"]
guard_lines[a:b]=z["old_text"].splitlines(keepends=True)
guard="".join(guard_lines).encode()
assert sha(guard)==manifest["intermediate_guard_sha256"]
assert len(ops)==2
olds, finals=old.decode().splitlines(keepends=True),final.decode().splitlines(keepends=True)
changes=[op for op in difflib.SequenceMatcher(None,olds,finals,autojunk=False).get_opcodes() if op[0]!="equal"]
assert changes==[("replace",171,172,171,172),("insert",225,225,225,239)]
assert finals[171].replace("\\mathcal","mathcal",1)==olds[171]
assert "\\tag{13a}" in "".join(finals[225:239])
assert "\\Gamma_i(\\mathcal N_h)\\subset\\Gamma_i(A)" in "".join(finals[225:239])
clip=lambda H,h:min(H,max(F(0),h))
scalar_cases=0
for H,hN,hB,e in product((F(0),F(1,2),F(1),F(2)),repeat=4):
    nloss=clip(H,hN)-clip(H,hN-e)
    aloss=clip(H,hN+hB)-clip(H,hN+hB-e)
    assert F(0)<=nloss<=e and F(0)<=aloss<=e
    assert clip(H,hN+hB-e)-clip(H,hN-e)<=hB
    scalar_cases+=1
H,hN,hB,e=F(1),F(1),F(1),F(1,2)
counterexample={"normal_loss":str(clip(H,hN)-clip(H,hN-e)),"full_loss":str(clip(H,hN+hB)-clip(H,hN+hB-e))}
assert counterexample=={"normal_loss":"1/2","full_loss":"0"}
result={"portable_rewrite":True,"original_script_sha256":"cc59a57aae645eb024e736c6332a62d2b8d8e4b903661d124f2adf9beac0de43","source_inputs":["final public proof","final-review-manifest.json"],"reconstructed_old_sha256":sha(old),"expected_old_sha256":manifest["old_sha256"],"reconstructed_guard_sha256":sha(guard),"expected_guard_sha256":manifest["intermediate_guard_sha256"],"final_sha256":sha(final),"exact_diff_operations":changes,"scalar_clip_guard_cases":scalar_cases,"same_eroder_counterexample":counterexample,"all_assertions_passed":True,"geometry_certified":False,"source_files_modified":False}
print(json.dumps(result,ensure_ascii=False,indent=2))
