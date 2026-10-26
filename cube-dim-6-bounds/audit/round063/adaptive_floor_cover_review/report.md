# F063 release review: implementation rejected; geometry retained

**REJECT_IMPLEMENTATION_FOR_RELEASE** for the frozen C063 preparation manifest `45f2defd70681f5d4cc826c58d2062f688bded2aca1d16dba7b487f22f4d287c`.

Seven independent configurations ran once under a hard 180-second alarm, completing in 0.03727674999390729 seconds. All author files stayed unchanged. Three forbidden transitions were actually accepted through public `append_snapshot` and `recover` calls with correct predecessor hashes:

1. An unresolved root changed to [0,1/2]^5 is accepted and still reports root coverage, although y=(1,1,1,1,1) is omitted.
2. A stale initial state appended after an inflight reservation is accepted. The reservation count falls from 1 to 0 and recovery reserves the identical `r.a0.n000001` token and input again.
3. A TEST_ONLY snapshot changed to MAIN_RELEASED with `release=None` is accepted and becomes reservable, bypassing the `authorize` helper's gate.

No snapshot bytes were rewritten or hashes forged. The defects concern missing state-transition validation, not the security of hashing or a claim that trusted main's actual initial root was wrong. There is no launcher, and no unauthorized native process or false global result was produced.

The mathematical tightening/emptiness/splitting proof remains correct under its stated root and legal-transition hypotheses. Ordinary API transitions enforce the eight-node/eight-attempt cap, reject mismatched cutoff-window jets, and retain passed leaves without recomputation. External native and exponential certificate authentication is clearly assigned to main; it must remain separate from metadata/hash syntax checks. Parent's final F061 acceptance closes that source gate only when bound by a genuine future release; the two finite reference cases do not become numerical reruns of future inputs.

Before release, a new core version must enforce the fixed root and immutable metadata, validate allowable state transitions, preserve topology and attempt/reservation histories, prevent orphan reissue, and prohibit test-to-production promotion through every public persistence path. Preserve this frozen version and issue a new seal for repair. The independent reasoning, exact code locations and repair contracts are in `proof_review.md`.

Production/native/kernel calls: **0**. Author twelve-case program runs: **0**. MATLAB and delegation: **0**. No global g, actual cover, or new C6 certificate. No stage, commit or push.

`results.json` records all seven scenarios, accepted snapshot hashes, the omitted root witness and identical orphan token/input. `verify_core.py` is the independent reproducer. Its synthetic constant-function fixtures are explicitly TEST_ONLY; their bounds use elementary rational inequalities and no exponential/native evaluation. Temporary state directories were created only here and removed after recording their outcomes.

Replay in a clean copy without results.json:

```sh
python3 -B verify_core.py
```

The script imports only the pinned pure core, never the author's check_fixed.py. It refuses to overwrite an existing completed review and checks all thirteen preparation pins before and after execution. `run_receipt.json` records the original successful audit process, and `manifest.json` seals the final review artifacts. Exit 0 means the seven registered audit scenarios executed as expected; it does not mean the reviewed implementation was accepted.
