# B015 spatial minimum diagnostic — preparation only

Status: **PREPARED; NOT RELEASED; ZERO NATIVE EVALUATIONS**.
Ownership is this directory only. The main agent must first review the final
C014 engine report and explicitly release this batch. Neither an engine report
appearing on disk nor elapsed time grants release. `search.py --run` refuses to
launch any process when `release.json` is missing or its pinned hashes fail.

## Fixed scope

Base body32 Cauchy frame, dimension 6, Lebesgue measure; single cutoff
T=6/125. Face x1=1 with x1>x2>...>x6>0. The cube maximal operator has common
coordinate radius, all positive scales, and general L1 data. This is a bounded
point diagnostic, not a global minimum computation or cube-cover certificate.
No new cutoff/frame, separate G pass, Monte Carlo, grid, MATLAB, subprocess
agent, commit, or push is allowed. Other agents own the engine and cell proofs.

Thresholds are K_T>T/4096=3/256000 for the proposed strict C6<64 route, and
m0=T/4000=3/250000 for the sufficient cost 62.5 route. A validated point with
K upper <=3/256000 immediately terminates the whole batch and rules out this
fixed cutoff's strict-64 route. K upper <3/250000 rules out the 62.5 target
only, and search for the 64 obstruction may continue within the original cap.
An exact equality at m0 does not disprove the non-strict m>=m0 claim; the driver
records the <= threshold flag separately from strict 62.5 rejection.

## Exact seeds and optimizer

At most two seeds, in this order:

1. B010 point2, componentwise nearest rational with denominator 10^14,
   using exact Fraction arithmetic and ties-to-even rounding:
   `[100000000000000,99999999000000,96784774597848,85268521812303,48714088834840,5157107002384]/100000000000000`.
   This is a new exact input; it inherits no earlier point certificate.
2. The A007 reference point, unchanged exact denominator 1000000007:
   `[1000000007,998836847,892738545,807464540,467291939,96968636]/1000000007`.
   Its gradient at this cutoff is a new field requiring a complete native pass.

SciPy SLSQP uses the five free coordinates, a supplied gradient, ordered gaps
at least 1e-8, x6>=1e-8, and a coordinate box of radius 1/20 about each seed.
The local box is a search restriction only, with no implication about points
outside it. Function and gradient are both multiplied by 1e5 for scaling.
Options are `maxiter=16, ftol=1e-10`; no finite-difference gradient calls.

Each seed has at most 8 NEW native attempts, including failed launches. There
are at most 16 total native attempts, 80 optimizer requests, and 1200 seconds
total wall time from release execution. An individual native child gets at
most min(120 seconds, remaining time). A seed cap may transfer control to the
second registered seed. Unused seed allowance is not added to another seed.
Optimizer failure, unsupported point, hash mismatch, interval failure, or any
other exception stops the whole batch with a checkpoint; there is no retry or
automatic third seed. Ordinary optimizer convergence is heuristic only.

The first callback uses its exact registered seed. Later binary64 proposals
are converted through their exact binary fractions to nearest denominator
10^14 integers; exact chain/box feasibility rounding handles boundary ulps.
If this changes a proposed free coordinate by more than 1e-10, stop. This is a
documented candidate construction, not a tolerance identification of inputs.
Every evaluated rational input is retained verbatim. Float midpoints are fed
only to SLSQP, never used for a mathematical inequality.

## Exact generic-node preflight

Read the audited 21-cone integer fan, extract its 18 distinct (q,ray) rows,
and generate every signed/permuted integer ray form. Divide each form by its
integer gcd and identify its global negative. The resulting 9308 forms have
exactly the same zero union as every native Fourier node. The use of all 64
signs here is harmless: the native 32-pair reduction already includes the
global negative. Zero tests are exact integer dot products with no tolerance.
Repeated NONZERO nodes remain supported by C014's confluent branch; they are
not excluded or skipped.

The first seed's smallest absolute primitive form numerator is 1000000;
the second's is 1163160. Both have zero zero-forms. At denominator 10^14 the
conservative raw numerator bound is 1.2e15 and den*q<=1.6e16, strictly inside
signed int64; the native cross-product comparator uses int128. Every future
candidate is rechecked before calling native code.

If a proposal has an exact zero node, quantize its exact input to denominator
10^14 and form the six exact gaps `(1-x2,x2-x3,...,x5-x6,x6)`. Choose the largest
gap j, with the smallest index on ties. For each k=1,...,6 independently from
this same base, add k*11^i/10^14 to every gap i other than j, and subtract their
sum from gap j. Reconstruct the five coordinates from these gaps. Each repair
is a NEW exact input. If the old c requires this step, its base rounding errors
are recorded and no earlier certificate is inherited. Check every gap>=1e-8,
face, local box, int64 and all zero forms before accepting any candidate. Each
failed candidate is wholly rejected and its exact failure reasons recorded;
there is no pseudo-penalty objective and no kernel evaluation at that input.
If none is valid and generic, stop as unsupported; do not return a mathematical
counterexample or skip terms. The six repairs are not a line search.

Memoization uses `(frameSHA, T, reduced-common nums, reduced-common den)`.
Exact duplicates reuse the same complete output; nearby points never merge.
The frame SHA includes the 32 unit-weight body directions, Fourier convention,
Lebesgue measure, and fan SHA. Prior G/H-only values do not supply a new K+grad
native result. The two planned native tuples are not the C014 reference tuple,
whose cutoff was 4/125.

## Native evidence and release

The provisional engine source, binary, fan, runtime libraries and existing
validation file hashes are recorded in `preparation.json`. Release additionally
must bind the final C014 report and the four prepared driver/registration files.
No engine is copied, rebuilt, modified, or executed during preparation.

Before each eventual pass, persist exact input, command, source/binary/fan SHA,
then the child PID and start time. Preserve complete stdout/stderr, exit status,
duration and hashes. Accept only precision512, full967680/signpairs483840,
zero-node count0, valid confluent count, all ordered finite dyadic intervals,
positive K and p, sufficiently narrow endpoints, and native Euler containment.
Reconstruct independently with Fraction the interval for
`x.grad K + 5 K + T p_T`, and require it to contain zero. All process handles
must be waited, including killed timeouts/interruption. An existing checkpoint
prevents automatic rerun/resume; parent must separately authorize any new batch.

No native evaluations are authorized by this file. Preparation may validate
the decoder on the EXISTING C014 output and prove that the release gate refuses
to run. Neither counts as a new point evaluation. Completion of 16 points with
no obstruction remains globally unresolved and cannot improve C6.
