# C160 handoff: prepare only, phases remain UNRUN

Sol wrote this packet without executing it, a mathematical preflight,
an optimizer, a root or a saved numerical checker. Administrative source
reads and hashes are the only executed operations. Source/code review
and Terra preparation must precede any primary science decision.

## Unique phase A prepare command

Use exactly this explicit independent target root. This command only
copies the source whitelist, derives the exact upward dyadic96 cap and
freezes phaseA inputs/registration/release requirements. It does not
launch optimization or import Joint. Run only after the primary's
prepare authorization. No mathematical validation is claimed from Sol.

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round160/normalized_joint_with_true_beta2/prepare_for_terra.py --root /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round160/normalized_joint_with_true_beta2 --phase A
```

Source whitelist: PROOF.md, common.py, prepare_for_terra.py,
run_phaseA.py, run_phaseB.py, HANDOFF_TO_TERRA.md,
STATIC_SOURCE_PINS.json, plus PACKAGE_MANIFEST.json (the manifest of
the preceding seven files, without a self hash). Preparation copies no
other code. External read-only mathematical sources are exactly the
role-addressed list in STATIC_SOURCE_PINS.json. No live129 or reserved
density is read. The source cap is the complete C156 result field,
not its display file. The original full128 candidate is the warm start.

## Required phase A actual gate

The generated `phaseA/release_requirements.json` is a contract, not
approval. The primary creates `phaseA/parent_release.json` containing
all exact fields of that contract plus:

* main_approved=true and independent_review_approved=true;
* source156_saved_audit.result_sha256 equal to the frozen C156 hash;
* source156_saved_audit.independent_saved_review_passed=true;
* source156_saved_audit.main_accepted=true;
* source156_saved_audit.evidence_files: a nonempty array of full path,
  SHA-256 and bytes rows identifying the independent terminal saved
  review/seal. Every cited file hash is checked before science.

The gate also explicitly accepts the complete base128 caps hash. Missing
C156 saved acceptance blocks science. A preflight review alone is
insufficient. The primary/operator must verify that the cited saved
audit actually establishes validity of the full C156 upper expression;
the program does not infer a PASS by searching arbitrary review prose.

After that separate gate, phaseA's *only* actual command is:

```sh
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python /Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round160/normalized_joint_with_true_beta2/run_phaseA.py
```

The runner enforces nice>=10, one process, one SLSQP, at most150 total
objective/gradient calls and no wall/CPU limit. It makes no strict M
claim, creates its own exclusive marker and never prepares or launches
phaseB. First error stops without automatic retry.

## Conditional phase B, separately selected and gated

Do not run this section as a continuation of preparation or phaseA.
The primary first examines the one actual new candidate. PhaseB
preparation requires phaseA's saved numerical comparison to improve
old128 and a primary `phaseB_selection.json` with these fields:

```json
{
  "round": "160",
  "main_selected": true,
  "phaseB_prepare_authorized": true,
  "phaseA_result_sha256": "SHA256_OF_ACTUAL_PHASE_A_RESULT",
  "candidate_canonical_sha256": "FIELD_FROM_ACTUAL_PHASE_A_RESULT",
  "evidence": "Concrete numerical margin and reason to certify this one candidate."
}
```

Then use the same prepare script with `--phase B --selection` and the
absolute path to that selected JSON, retaining the exact same `--root`.
This only derives exact a=x/cap, verifies cost accounting, and registers
15 disjoint prefix jobs. It makes no H/root call. A *new* independent
preflight and primary release must bind the generated
`phaseB/release_requirements.json`, selection/candidate/job hashes and
the same accepted C156 saved review. There is no inherited actual gate.

Only after that second gate may Terra run `run_phaseB.py` with the same
single-thread environment. Fifteen worker processes each inherit nice10;
there is no time cap. The batch is one new rational-a candidate with
735471 prefix plus16 tail logical root requests. No old root cache from
a different a is consumed. Existing CDF endpoint laws and symbolic
feature tables are reused; all probability mass is retained.

The outputs are exact comparisons of H_U^5 with old128, Theta(g0) and
Theta(gstar), with cost<=1. The complete M-slot certificate still needs
its independent saved audit. Passing it proves neither global g nor
any new C6 bound. Do not change the caps/precision/cutoff or launch a
second optimization on failure without a new registered task.
