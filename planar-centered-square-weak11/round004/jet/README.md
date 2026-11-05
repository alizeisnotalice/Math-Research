# Calibrated signed first jet

This is a portable interface for the previously audited finite jet calculation. The 6811 R003 records are regenerated rather than packaged. The project root needs `round003/mechanism/`, `round002/verification/candidate002.json`, `round004/weights/candidate_05.check2.json`, and the existing `checks/lower/lower_check.py`.

From the project root, generate and audit in a separate writable directory:

```sh
python3 -u round003/mechanism/decompose.py --output-dir /tmp/r003-generated
python3 -u round004/jet/jet.py --records-dir /tmp/r003-generated --output-dir /tmp/r004-generated
python3 -u round004/jet/audit.py --records-dir /tmp/r003-generated --output-dir /tmp/r004-generated
OPENBLAS_NUM_THREADS=1 python3 -u round004/jet/check_calibrated.py --records-dir /tmp/r003-generated --output-dir /tmp/r004-generated
wolframscript -file round004/jet/scalar.wl /tmp/r004-generated
```

Without options, the default records directory is `round003/mechanism/generated/`, and output is `round004/jet/generated/`; an explicit output directory avoids touching previously recorded results. NumPy is required by the independent coordinate-compression core.

The strict calibrated finite family is based on candidate A and uses the signed mass direction 34→27. For each retained bbox, rhohat²=rhoD²+t*w/4 preserves the old strictly positive gap to its actual captured mass. The script lexicographically clips/sweeps endpoint jets, including seams and zero-width jets. The exact score is in `jet.json`, numerically dA/dt≈0.1865589932; total mass stays fixed. This is not the derivative of the full ideal halo and does not assert a global optimum. The `stats.declipped` entry is an unused prototype counter fixed at zero; it provides no information about clipping.

The portable run's 601-character exact score and every mathematical field match the audited original; see `portable-receipt.json`. It also reran direct copy-count checks, the independent strict t=10^-9 evaluation, and Mathematica assertions. Their actual portable outputs are `portable-audit.log`, `portable-check.log`, and `portable-scalar.log`; `portable-jet.log` is retained as compact replay stdout; the full generated geometry records remain local and can be regenerated. Source and dependency hashes are in `SHA256SUMS`.

The second exact derivative implementation is `independent_verify.py`; see [its replay guide](independent-README.md). Its arithmetic source is byte-identical to the independently relocated script. The independent receipt and diff use distinct names to avoid overwriting the primary replay evidence.
