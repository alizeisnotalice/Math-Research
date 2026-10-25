#!/usr/bin/env python3
"""Screen only the 100 length-five indices omitted by the prior capped run."""
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import sys
import time

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import length5_screen as S


def main():
    started = time.monotonic()
    all_new = [ell for ell in product(range(1, 6), repeat=6) if 5 in ell]
    all_new.sort(key=lambda e: (-prod(e), e))
    omitted = all_new[S.NEW_ELL_CAP:]
    assert len(all_new) == 11529 and len(omitted) == 100
    a = S.spacing()
    # Reference data are constructed only to test the newly omitted boxes.
    old = [(ell, S.witness_box(a, ell)) for ell in product(range(1, 5), repeat=6)]
    old_keep = S.maximal([(ell, box) for ell, box in old if box is not None])
    records = [(ell, S.witness_box(a, ell)) for ell in omitted]
    nonempty = [(ell, box) for ell, box in records if box is not None]
    direct = [(ell, box) for ell, box in nonempty
              if any(S.contained(box, oldbox) for _, oldbox in old_keep)]
    non_direct = [(ell, box) for ell, box in nonempty
                  if not any(S.contained(box, oldbox) for _, oldbox in old_keep)]
    out = {"status": "completed_numeric_omitted_only", "full_new_index_count": len(all_new),
           "previously_screened": S.NEW_ELL_CAP, "screened_here": len(omitted),
           "omitted_list_sha256": hashlib.sha256(repr(omitted).encode()).hexdigest(),
           "new_raw_nonempty": len(nonempty), "new_numerically_directly_contained": len(direct),
           "new_not_numerically_single_box_contained": [list(ell) for ell, _ in non_direct],
           "elapsed_seconds": time.monotonic() - started,
           "code_sha256": hashlib.sha256(HERE.read_bytes()).hexdigest(),
           "reference_algorithm_sha256": hashlib.sha256((HERE.parent / 'length5_screen.py').read_bytes()).hexdigest()}
    HERE.with_name('length5_remaining100_results.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
