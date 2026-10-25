# Lower round 003: bounded spacing-perturbation screen

**Frozen comparator.** The candidate's reliable lower endpoint must strictly
exceed the full `R6-001` exact upper endpoint in the current main ledger,
not its `2.06513594276568` display.  The read-only ledger SHA-256 and the
reference verifier SHA-256 are recorded by the screen output.

**Family.** Start with the certified one-atom rectangular lattice and its
96 retained common-radius witness blocks.  Vary only its six spacings by
rational perturbations.  This keeps the Lebesgue setting, one common cube
radius per witness, and the same finite-block/smoothing transfer interface.
The screen recomputes boxes and containment; it does not assume the old
96-box union remains valid.

**Finite candidates (64 maximum).** Evaluate 16 predeclared rational
perturbations: `p0,p1 in {-1,1}/500` for the first two spacings and the
four tail spacings receive either all `+1/500` or all `-1/500`, or the
two-pair pattern `(s,s,-s,-s)` with `s in {-1,1}/500`; the four choices
are combined with the four `(p0,p1)` choices.  Candidate zero and any
previous round-001 point are not evaluated.  The screen has a 20-minute
wall-time stop; it stops earlier if all 64 finish.  It records a numeric
value only as a discovery diagnostic.

**Promotion.** Only a candidate whose floating diagnostic clears the
frozen full exact upper endpoint by at least `10^-7` will be promoted to a
fresh `Q(u,v)`/rational-perturbation exact union certificate, followed by
the existing finite-L1 transfer proof with its own strict margins.  No
candidate is a bound before that certification.

**Acceptance.** Report every tested candidate and the maximum diagnostic,
or the explicit screened range if none promotes.  No optimizer, random
restart, satellite sweep, or enlarged witness grid is authorized by this
batch.
