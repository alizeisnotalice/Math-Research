#!/usr/bin/env python3
"""Read-only replay of the main candidate verifier into this audit directory.

The audited program is executed unchanged in its arithmetic and geometric
logic.  Only its root for read-only ledgers and its generated-result path are
redirected, so this replay does not write the main clone.
"""
from pathlib import Path

AUDIT = Path(__file__).resolve().parent
MAIN = Path(__file__).resolve().parents[3]
SOURCE_PATH = MAIN / 'audit/lower003/rational_inward_candidate.py'
SOURCE_BYTES = SOURCE_PATH.read_bytes()
source = SOURCE_BYTES.decode()
needle = "ROOT=Path(__file__).resolve().parents[2]"
assert source.count(needle) == 1
source = source.replace(needle, f"ROOT=Path({str(MAIN)!r})")
source = source.replace(
    "Path(__file__).with_name('rational_inward_candidate_results.json')",
    "Path(__file__).with_name('audit_replay_main_candidate_results.json')",
)
# Preserve the hash of the reviewed main source even though __file__ points
# here to keep the program-generated JSON inside this independent clone.
needle = "hashlib.sha256(Path(__file__).read_bytes()).hexdigest()"
assert source.count(needle) == 1
source = source.replace(needle, "hashlib.sha256(SOURCE_BYTES).hexdigest()")
globals_for_replay = {
    '__name__': '__main__',
    '__file__': str(AUDIT / 'rational_inward_candidate.py'),
    'SOURCE_BYTES': SOURCE_BYTES,
}
exec(compile(source, '<read-only-main-candidate-replay>', 'exec'), globals_for_replay)
