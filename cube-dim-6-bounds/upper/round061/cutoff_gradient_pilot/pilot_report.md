# C061: two cutoff-gradient pilot centers and two whole boxes pass

The released two-call primary branch completed successfully. Both the fixed A038
low point and the prescribed genuine S3 point have a certified positive-band
center interval strictly above `1/65000`; both yield a passing whole face box at
the second predeclared halfwidth `h=1/1000`. This is a finite pilot with two local
boxes, not a global g-floor certificate, a new C6 bound, or closure of the main
partial-I argument.

All data use the original unit-weight body32 Cauchy model and the same band
`F=B_(1/250,1/2)`. Each center's value, all six physical derivatives and endpoint
densities came from one unchanged C016 Arb512 native call. Its five face
derivatives were selected from that same result. The whole-box certificates use
the F057-reviewed A057 fixed quadratic `Q=I5/32`, with the exponential evaluated
at Arb384 and rational interval preprocessing. No full-G gradient, other
cutoff, other model or new Hessian calculation is used.

The following decimals are display only; exact rational/dyadic endpoints are
retained in `pilot_summary.json`, the native/jet files and the box certificates.

| Case | Center band value, display only | Lower bound at h=1/500 | Lower bound at h=1/1000 | Final box |
| --- | ---: | ---: | ---: | --- |
| A038 low point 3 | 1.5775912155370232e-5 | 1.5278877211203848e-5 | 1.5677763798311330e-5 | PASS |
| S3 `(1,1,1,1,3/5,0)` | 1.5728260114371485e-5 | 1.5236195399366246e-5 | 1.5632187853588424e-5 | PASS |

The exact target is `1/65000` (approximately `1.5384615384615384e-5`). At the
larger halfwidth `1/500`, the chosen tangent certificate is inconclusive; this
does not prove that the true band drops below target. At halfwidth `1/1000`,
every outward lower endpoint at all 32 corners is at least the exact target.
Each passing box consists of the five intervals
`[max(0,c_i-1/1000),min(1,c_i+1/1000)]`, for physical coordinates `i=2,...,6`.
For S3 these intervals are exactly
`[999/1000,1]` three times, `[599/1000,601/1000]`, and `[0,1/1000]`.
The theorem covers the whole box, including ties, coordinate zero and projection
collisions; it does not claim that the two boxes cover the remaining domain.

The low center is the exact A038 `points[2]`, read from its pinned original
result, with denominator `100000000000000`. Its source identity and full
numerators are retained in the original registration, preparation result and
same-window jet. The prescribed S3 center was already checked in preparation
to have a rank-three body-column representation and to lie outside the paid
S2 tube.

The low-primary native reported 967680 logical group terms, 483840 sign-pair
evaluations, 16931 distinct exact nodes and no confluent/zero-node terms. The
S3-primary native reported the same complete group counts, 347 distinct exact
nodes, **405600 confluent terms and 191808 zero-node terms**. Thus the actual S3
call exercised the existing Hermite and zero-projection path at full-group
scale; these terms were not dropped. Both natives returned finite dyadic
value/gradient/endpoint intervals and an Euler interval containing zero. A
separate exact-rational interval reconstruction of
`x.grad(F)+5F+T*p_T-a*p_a` also contained zero, and both endpoint-density lower
bounds were positive.

The inherited C016 limitation remains explicit: the full six-component gradient
has no independent full-group cross-certificate. The accepted implementation
evidence consists of the saved independent scalar DD5/6 checks, the physical
chain formula, exact arithmetic and consistency checks. Euler by itself does
not establish every gradient component. These pilot results do not claim to
close that independent-gradient validation gap.

The actual main release is `production_release.json`, SHA256
`70dd963fcf2bbfed135ff309832af2c78f7fc3eedabc6cc2fbf5307e12bca2bd`.
It binds the seven ready files, the accepted F057 PASS evidence and the
180-second policy. The original preparation manifest and its files remain
unchanged. The first native's strict center pass selected S3-primary directly;
neither fallback case ran. Each center reused its single jet for two whole-box
checks, giving four boxes and 128 corner evaluations in total, with zero
additional kernel calls.

| Actual case | Original tool session | Native PID | Native exit | Native seconds | Wrapper wall seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| low-primary | 28851 | 63619 | 0 | 4.56374 | 4.653933749999851 |
| S3-primary | 94531 | 63705 | 0 | 7.40414 | 7.453802917007124 |

Both original tool sessions completed normally and were observed through their
same handles. The native PIDs were absent in the final read-only process check.
Each complete runner execution was below 180 seconds. The original tool outputs
are in `tool_execution_receipts.json`; native commands, exact input metadata,
PIDs, exit codes, source/binary/fan/release hashes and timestamps are in each
case's `receipt.json`. Native stdout, stderr and each attempted box certificate
are retained, including the inconclusive first halfwidth. No process was killed
or restarted, and no other agent's process was touched.

Totals: 2 native kernels, 0 fallback kernels, 0 retries, 0 timeouts, 0 new
compiles, 0 old-seed recalculations, 0 old failed-window reruns and 0 C055 reruns.
Native reported times sum to about 11.96788 seconds; wrapper wall times sum to
12.107736667006975 seconds. There were no MATLAB calls, delegation, staging,
commits or pushes. The finite branch is complete and no further kernel is
authorized by this completed pilot. Main retains the proof, coverage and
ledger integration decisions.
