#!/usr/bin/env python3
"""Read-only replay of the main k16 inward-rational checker in this clone."""
from pathlib import Path

AUDIT = Path(__file__).resolve().parent
MAIN = Path(__file__).resolve().parents[3]
SOURCE_PATH = MAIN / 'audit/lower003/rational_inward_k16.py'
SOURCE_BYTES = SOURCE_PATH.read_bytes()
source = SOURCE_BYTES.decode()
needle = "ROOT=Path(__file__).resolve().parents[2]"
assert source.count(needle) == 1
source = source.replace(needle, f"ROOT=Path({str(MAIN)!r})")
source = source.replace(
    "Path(__file__).with_name('rational_inward_k16_results.json')",
    "Path(__file__).with_name('audit_replay_main_k16_results.json')",
)
needle = "hashlib.sha256(Path(__file__).read_bytes()).hexdigest()"
assert source.count(needle) == 1
source = source.replace(needle, "hashlib.sha256(SOURCE_BYTES).hexdigest()")
exec(compile(source, '<read-only-main-k16-replay>', 'exec'), {
    '__name__': '__main__',
    '__file__': str(AUDIT / 'rational_inward_k16.py'),
    'SOURCE_BYTES': SOURCE_BYTES,
})
