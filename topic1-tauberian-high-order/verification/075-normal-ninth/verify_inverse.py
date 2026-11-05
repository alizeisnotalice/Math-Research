"""Portable byte-level inverse verifier for the R075 shared-eroder patch.

Run from the staged/public repository tree. No source file is modified.
"""
from pathlib import Path
import hashlib, json, re
root=Path(__file__).resolve().parents[2]
proof=(root/"proofs/075-normal-ninth.md").read_bytes()
diff=(Path(__file__).resolve().parent/"shared_eroder.diff").read_bytes()
manifest=json.loads((Path(__file__).resolve().parent/"shared_eroder_check.json").read_text())
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(proof)==manifest["revised_sha256"]
text=diff.decode().splitlines(keepends=True)
hunks=[]; i=0
while i<len(text):
    if not text[i].startswith("@@ "):
        i+=1; continue
    m=re.match(r"@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@",text[i])
    assert m, text[i]
    old_start=int(m.group(1)); old_count=int(m.group(2) or 1)
    new_start=int(m.group(3)); new_count=int(m.group(4) or 1)
    i+=1; expected=[]; replacement=[]
    while i<len(text) and not text[i].startswith("@@ "):
        line=text[i]
        if line.startswith(" "):
            expected.append(line[1:]); replacement.append(line[1:])
        elif line.startswith("+"):
            expected.append(line[1:])
        elif line.startswith("-"):
            replacement.append(line[1:])
        elif line.startswith("\\ No newline"):
            pass
        elif line.startswith("--- ") or line.startswith("+++ "):
            break
        else:
            raise AssertionError(f"unrecognized unified-diff line: {line!r}")
        i+=1
    hunks.append((new_start,new_count,expected,replacement))
lines=proof.decode().splitlines(keepends=True)
for new_start,new_count,expected,replacement in reversed(hunks):
    at=max(0,new_start-1)
    assert lines[at:at+new_count]==expected, (new_start, lines[at:at+new_count], expected)
    lines[at:at+new_count]=replacement
restored="".join(lines).encode()
result={"portable":True,"hunks_reversed":len(hunks),"final_sha256":sha(proof),"restored_old_sha256":sha(restored),"expected_old_sha256":manifest["original_sha256"],"restored_old_bytes":len(restored),"expected_old_bytes":manifest["original_bytes"],"inverse_byte_hash_match":sha(restored)==manifest["original_sha256"],"source_modified":False}
assert result["inverse_byte_hash_match"] and len(restored)==manifest["original_bytes"]
print(json.dumps(result,ensure_ascii=False,indent=2))
