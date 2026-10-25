# Lower-bound batch registration: round 001

**Objective.** Evaluate one prescribed six-dimensional instance of the
locked radical rectangular-lattice family and either obtain an exact
Lebesgue `L1` lower certificate or stop this family for the round.

**Finite batch.** Exactly one parameter tuple is permitted:
`(1+u^4/2,(2+u^5v)/4,(1+u)/2,(1+u)/2,(1+u)/2,(1+u)/2)` with
`u^6=2,v^6=3`; all count vectors lie in `{1,2,3,4}^6`.  No optimizer,
random search, satellite perturbation, or enlarged grid is in this batch.

**Acceptance.** Keep only nonempty witness boxes passing all 64 common-
radius corner containments; compute the full six-dimensional union in the
exact number field with root isolation; prove finite-block smoothing and
compare the resulting exact rational lower endpoint with the inherited
four-dimensional full certificate endpoint.

**Stop condition.** A successful exact run is immediately frozen for an
independent audit.  A failed containment, ambiguous sign, or an endpoint
not exceeding the inherited seed is reported as a failed candidate; no
parameter sweep follows automatically.

**Outcome.** The run exited zero.  It retained 96 of 127 candidate boxes,
gave the strict lower endpoint `206513594276568/10^14`, and passed the
stated acceptance criterion.  See `round001_certificate_result.json` and
`round001_transfer_proof.md`.
