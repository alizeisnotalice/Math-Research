# Independent exact calibrated derivative reproduction

Place `independent_verify.py` beside the published `round004/jet/jet.py` and
`jet.json`. It reads the published `round004/weights/candidate_05.json` directly.
Python with NumPy is required.

From the project root:

```sh
python3 -u round003/mechanism/decompose.py --output-dir /tmp/r003-generated
OPENBLAS_NUM_THREADS=1 python3 -u round004/jet/independent_verify.py --records-dir /tmp/r003-generated --output-dir /tmp/r004-independent
```

Defaults are `round003/mechanism/generated/` for records and
`round004/jet/generated/` for output. The verifier writes
`independent_verify.json` only to the selected output directory and prints the
same result. It does not regenerate the R003 records or rerun either finite
candidate area checker.

This retains the previously audited independent algorithm: direct counts from
the infinite periodic lattice and a two-dimensional compressed coverage grid,
instead of the primary jet program's segment tree. Only paths, command-line
arguments and output location changed; see `independent-portable.diff`.

The actual relocated run used regenerated R003 records and returned the same
601-character rational derivative, baseline area, strict gap and radicand bound.
Its algorithm time was 0.919 seconds; process time including Python startup was
1.173 seconds. See `independent-portable-receipt.json` and `independent_verify.json`.
This is the derivative of the fixed strict calibrated receiver family, with the
scope of the original R004 audit.

For integration, copy the script to `round004/jet/independent_verify.py`. The
source, compact actual output, log, receipt and diff may be retained as evidence.
The local `replay/` directory and regenerated 2.85 MB records are not publication
inputs.
