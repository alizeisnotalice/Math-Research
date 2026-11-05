"""Read-only, portable inverse check for the three R078 notation edits."""
from pathlib import Path
import hashlib
import json
import re

here = Path(__file__).resolve().parent
project = here.parents[1]
proof = (project / "proofs/078-end-face-reserve.md").read_bytes()
diff_bytes = (here / "notation.diff").read_bytes()
meta = json.loads((here / "freeze.v2.json").read_text())
old_meta = json.loads((here / "freeze.json").read_text())
sha = lambda data: hashlib.sha256(data).hexdigest()
assert sha(proof) == meta["source"]["revised_sha256"]
assert sha(diff_bytes) == meta["source"]["diff_sha256"]

new = proof.decode("utf-8").splitlines(keepends=True)
lines = diff_bytes.decode("utf-8").splitlines(keepends=True)
assert lines[0].startswith("--- ") and lines[1].startswith("+++ ")
old_out = []
new_pos = 0
i = 2
hunks = 0
while i < len(lines):
    if not lines[i].startswith("@@ "):
        i += 1
        continue
    m = re.match(r"@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@", lines[i])
    assert m, lines[i]
    old_start = int(m.group(1))
    old_count = int(m.group(2) or "1")
    new_start = int(m.group(3))
    new_count = int(m.group(4) or "1")
    target = new_start - 1
    assert target >= new_pos
    old_out.extend(new[new_pos:target])
    new_pos = target
    i += 1
    old_seen = new_seen = 0
    while i < len(lines) and not lines[i].startswith("@@ "):
        line = lines[i]
        if line.startswith("\\ No newline at end of file"):
            i += 1
            continue
        assert line and line[0] in " +-", repr(line)
        body = line[1:]
        if line[0] == " ":
            assert new_pos < len(new) and new[new_pos] == body
            old_out.append(body)
            new_pos += 1
            old_seen += 1
            new_seen += 1
        elif line[0] == "-":
            old_out.append(body)
            old_seen += 1
        else:
            assert new_pos < len(new) and new[new_pos] == body
            new_pos += 1
            new_seen += 1
        i += 1
        if old_seen == old_count and new_seen == new_count:
            break
    assert old_seen == old_count and new_seen == new_count
    assert old_start > 0
    hunks += 1
old_out.extend(new[new_pos:])
old = "".join(old_out).encode("utf-8")
assert hunks == meta["source"]["replacements"] == 3
assert len(proof) == meta["source"]["revised_bytes"]
assert len(old) == old_meta["source"]["bytes"]
assert sha(old) == meta["source"]["original_sha256"]
assert sha(old) == old_meta["source"]["sha256"]
result = {
    "read_only": True,
    "executed_command": "python3 verify_inverse.py > verify_inverse.json",
    "checker_sha256": sha(Path(__file__).read_bytes()),
    "source_inputs": {
        "final_proof_sha256": sha(proof),
        "diff_sha256": sha(diff_bytes),
        "metadata_sha256": sha((here / "freeze.v2.json").read_bytes()),
        "original_freeze_sha256": sha((here / "freeze.json").read_bytes()),
    },
    "hunks_reversed_in_memory": hunks,
    "restored_original_sha256": sha(old),
    "expected_original_sha256": meta["source"]["original_sha256"],
    "restored_original_bytes": len(old),
    "expected_original_bytes": old_meta["source"]["bytes"],
    "inverse_byte_hash_match": sha(old) == meta["source"]["original_sha256"],
    "source_files_modified": False,
}
print(json.dumps(result, ensure_ascii=False, indent=2))
