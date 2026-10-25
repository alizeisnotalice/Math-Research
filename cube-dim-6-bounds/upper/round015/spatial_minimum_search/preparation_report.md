# B015: bounded spatial diagnostic prepared, not released

**Current result: preparation PASS, zero native evaluations.** No value of K,
gradient, or p has been computed at a new point in B015. No global minimum,
lower bound on K over the cube, or improvement of C6 is claimed. C6<=64 remains
the existing upper bound. The main agent's explicit release, after its C014
final review, is still required; this report does not grant it.

The prepared batch fixes T=6/125 and uses at most two seeds, at most eight
new native attempts per seed, at most sixteen total, and twenty minutes after
release. It is a SciPy SLSQP candidate search in five free coordinates, with
exactly x1=1, ordered gaps at least 1e-8 and x6 at least 1e-8. Each seed has a
local coordinate box of radius 1/20. K and its gradient are scaled by 1e5 only
inside the optimizer. These local boxes make no assertion about their
complement. The full pre-registered construction, limits, and stops are in
`registration.md`; `search.py` implements them.

The first exact seed is

```
[100000000000000, 99999999000000, 96784774597848,
 85268521812303, 48714088834840, 5157107002384] / 100000000000000.
```

It is the exact nearest-denominator-1e14 rounding of B010 point2, with
ties-to-even rounding. Each rounding error is recorded as an exact Fraction
and bounded by 1/(2e14). At least one error is nonzero, so no old point
certificate is inherited. The second seed retains the exact old A007 point
and denominator 1000000007. Its gradient at T=.048 is new data, not supplied
by B011's earlier H-only pass. The native full output currently read from
C014 concerns T=.032, not this proposed cutoff.

An exact preflight generated the complete signed/permuted zero loci from the
fan's 18 distinct rays. Modulo integer gcd and global sign there are 9308
forms. Both seeds avoid every form exactly. Their minimum nonzero primitive
form numerators are 1000000 and 1163160. For the first seed and all proposed
denominator-1e14 cube points, the conservative raw-node numerator bound is
1.2e15 and the largest raw denominator is 1.6e16, within signed int64. Each
future candidate is checked again; no near-zero tolerance is used. The
preflight detects 1321 zero forms at the deliberately degenerate all-ones
test input without evaluating any kernel there.

The driver records exact proposal construction and permits only six gap-space
zero-node repairs. With gaps `(1-x2,x2-x3,...,x5-x6,x6)`, choose the largest gap
j (smallest index on ties), add k*11^i/1e14 to each other gap i, and subtract
their sum from j, for k=1,...,6 independently from the quantized base. Every
repair is a distinct exact input, checked for cube/order/local-box/range and
genericity, with every rejection recorded. An old c repair explicitly changes
its exact input and inherits no certificate. Six failures stop as unsupported;
there is no pseudo-penalty objective. C014's confluent branch remains available
for repeated nonzero nodes. Neither zero nor repeated terms are skipped.

The future point evidence must be the full native512 dyadic output, with
967680 group terms, 483840 sign pairs, zero_node_terms=0, finite ordered narrow
intervals for K, all six gradient coordinates and p_T, positive K and p_T, and
Euler containment. The driver additionally reconstructs the exact rational
endpoint interval of x.grad K+5K+T p_T. Every stdout, stderr, input, source and
binary hash, PID, return code and duration will be preserved. Repeated exact
points are memoized; nearby points are not. Exceptions stop the batch and
preserve its checkpoint. There is no automatic rerun or precision escalation.

The strict-64 cutoff test is K upper <=3/256000, followed by immediate stop.
For the cost-62.5 target m>=3/250000, strict rejection requires K upper
<3/250000; exact equality is only a threshold contact. The driver records
both facts separately and may continue within the same budget after only a
62.5 rejection. Budget exhaustion or optimizer convergence without a
counterexample leaves the global question unresolved.

## Actual preparation verification

`prepare.py` ran successfully (tool chunk fddef0, exit0) and wrote only this
directory. It read and hashed the old point inputs, 21-cone fan, current C014
source/binary/build/validation output, and runtime libraries. The frame/fan
convention has SHA256
`9ac81d7ab65d06da48eb5b8746ac2a9803a6eec0d98cd7801db60847134d10b5`.

The initial `check_preparation.py` ran successfully (PID31880, tool chunk dfa98c,
exit0, about0.067 seconds). After the parent requested gap-space repairs, the
final checker was rerun (PID32209, tool chunk039ad7, exit0, about0.086 seconds;
see `checked_preparation.json` and `manifest.json`).
It checked syntax, all source hashes, seed rounding and
range/genericity, exact-key distinction, and the decoder's Euler computation
on the existing C014 result. It actually called the Python driver's locked
`--run` entry: exit3 and `NOT_RELEASED: no release.json; zero native calls`.
Neither checkpoint nor native output directory was created. No engine binary
was executed and no background process remains. No pycache was created in
old or new directories.

Dependencies actually inspected: Python3.15.0, SciPy1.18.1, NumPy2.5.4,
python-flint0.9.0; the registered native engine binds FLINT3.6.0 and its
runtime hashes. Skill L03 was used for separating exact input/output evidence
from heuristic proposals; that skill supplies no independent certification
of an interval API or of C014 itself. Engine mathematical correctness remains
subject to the parent/F review already in progress.

`release.json` does not exist. Release must bind the final C014 report SHA
and the current preparation/driver/registration hashes, with the parent's
explicit authorization. All prior B008–B011 directories remain unchanged.
