# R003 finite receiver mechanism: portable reproduction

Place these scripts under `round003/mechanism/` in a project root containing the published `checks/lower/` and `round002/verification/` files. Run from that root:

```sh
python3 -u round003/mechanism/decompose.py --output-dir /tmp/r003-reproduction
OPENBLAS_NUM_THREADS=1 python3 -u round003/mechanism/independent_verify.py --output-dir /tmp/r003-reproduction
wolframscript -file round003/mechanism/scalar.wl /tmp/r003-reproduction
```

Omitting `--output-dir` places the generated files in `round003/mechanism/generated/`. The second script independently resums the support mass of every saved nonempty receiver and uses the coordinate-compression union algorithm. It requires NumPy. The first script writes regenerable `decomposition.json` with all 6811 records and compact `summary.json`; only the latter is included in the portable package. Each bbox is closed and captures all old and new periodic atoms in its stated coordinate interval; receiver radii are strictly below the admissible roots and then clipped to the fixed period cell.

This reproduces the finite selected bbox family at reach=1396/125 and radius discretization Dr=2^40, not the entire true halo. Let A0 be the old seed union, A1 the union after updating all old-coordinate bboxes with positive old receiver area, A12 the additional union of those whose old receiver was empty, and Afull the additional union of bboxes with at least one new-coordinate boundary. Exact fractions, counts and record hash are in `summary.json`. A0 and Afull agree exactly with the previously certified seeds. The ordered area differences reflect overlaps and are not intrinsic or unique attribution. The count zero for old-empty activation is supported by the exhaustive enumeration in `decompose.py`; the independent verifier checks all saved nonempty receivers but does not independently enumerate the unrecorded empty bboxes.

`portable-decompose.log`, `portable-independent.log` and `portable-scalar.log` record an actual run from a separate root layout. The original and replay runs have identical mathematical fields and identical ordered records, with SHA256 `eafcc4722314fb28c3bc6c2eba1cd47992198ce4b5a6814b04a3b542d57ca7fc`; `portable-receipt.json` gives the comparison. The original and replay full JSON file hashes differ due to actual worker PID and elapsed-time metadata; neither full generated JSON file is required for publication. The raw replay is retained locally under `replay/` and is excluded from the publication whitelist. `SHA256SUMS` records the published sources, outputs, and identity of the public inputs and checker sources.

The optional `two_atom_first_order.wl` verifies only the algebraic implications of a separately established two-atom receiver geometry. It is not used for the 51-atom union decomposition.
