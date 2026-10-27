# One Terra execution, two separate gates

Sol supplies only five source files plus the code whitelist/manifest. Do not run
this directory in Sol's clone. Copy the exact approved source bytes into the fresh
Terra root:
`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/audit/upper160_saved`.
Do not change producer inputs, roots, coefficients or receipts.

After main approves prepare-only, Terra runs the registered Python with `-B prepare.py`
from that root. This freezes metadata/hashes only, and does not import an evaluator,
run Fraction/mathematical checks, compile libraries, create roots or verify inequalities.
Send the actual inputs/source_manifest/registration/prepare_receipt hashes for final
independent and main review. Preparation is not mathematical PASS.

Only after a concrete dual gate may main write parent_release.json with every field
returned by check_saved.release_requirements(), plus main_approved=true and
independent_review_approved=true. Do not call the helper to generate the gate in Sol;
read its definition and the prepared registration instead. The one authorized
execution is the registered Python `-B observe.py`, nice10, all BLAS/OMP thread
variables1. There is no timeout. The observer starts one checker, and that checker
starts exactly15 fixed local workers. These are local processes, not recursive agents.
Never run `check_saved.py --worker` manually. Never delete markers or automatically
retry a failed or interrupted batch.

Preserve observer_started/original_process/observer_terminal_receipt; checker start,
exit/certificate/audit_model; each issued worker's intent/process/start/checkpoint,
result/exit/reap and stdout/stderr. A missing record is UNKNOWN, not a completed job.
The output is a saved certification of one new M candidate only. It does not accept
a new currentbest by itself: main must read the actual terminal certificate and make
that scientific acceptance separately. No C6/global-g assertion follows.
