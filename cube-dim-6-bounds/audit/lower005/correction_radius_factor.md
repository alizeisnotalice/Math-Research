# Correction: invalid initial containment conclusion

The initial `containment_check.py` incorrectly used the interval for
`(prod ell)^(1/6)` as the cube radius.  The actual common radius is half
that value.  Thus its original result file
`containment_check_results.json` (SHA-256
`871f563663982d579c2b189b84438b0f2b14ffe99f0d5ed7cc8a74650c31aff8`)
and the statement that every new cap-5/cap-6 box was contained are invalid
and withdrawn.

The span-pruning result is unaffected: its condition uses the cube diameter
`2r=(prod ell)^(1/6)` and is therefore correctly normalized.  The repaired
containment program explicitly divides the sixth-root interval by two and
writes a separate `containment_check_corrected_results.json`, preserving the
erroneous artifact for auditability.  Only the corrected result may be used.
