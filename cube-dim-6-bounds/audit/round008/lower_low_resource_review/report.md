# Round 008: independent low-resource lower-candidate acceptance

**PASS.** For the single supplied candidate, the complete reviewed chain
proves

\[
C_6>\frac{2067397333267}{10^{12}}=2.067397333267.
\]

Here C6 is the weak-(1,1) constant for Lebesgue measure, centered
axis-parallel cubes with a common radius in all coordinates, the ordinary
general L1 class, and the supremum over all positive radii. This is a new
instance of the existing ordinary-lattice construction, not a new mechanism.
The reviewer did not participate in the round008 construction. This audit
made one primary integer replay, no extra union scan, no search, no MATLAB
call, no delegation, and no stage/commit/push. Lower work stops after this
acceptance, preserving the user's upper-bound priority.

## 1. Candidate and source geometry

Let u,v be the positive roots u^6=2 and v^6=3, and set

\[
a_1=1+u^4/2,\qquad
a_2=u^2v^2/3+u^2v^4/6,\qquad
a_3=\cdots=a_6=(1+u)/2.
\]

The source is one **unit** atom at each point a*k, k in Z^6. Its
fundamental cell has source mass one. For each ell in {1,2,3,4}^6, use the
central consecutive index block with coordinate start
q=-floor(ell/2) for odd ell and q=1-ell/2 for even ell. Its mass is
m=product(ell), and its one common radius is r=m^(1/6)/2. This same r is
used in all six coordinates; no independent rectangle radii occur.

In the positive half cell H=product[0,aj/2], intersection of the source
containment slabs gives odd-length interval
[0,min(aj/2,r-(ell-1)aj/2)] and even-length interval
[max(0,ell*aj/2-r),aj/2]. A center in the resulting box is within r of
every selected source atom in every coordinate, hence its cube contains
the selected m masses. Its volume is (2r)^6=m, giving average at least one.
Additional atoms only increase that average.

The rectangular unit lattice is invariant under all six coordinate sign
reflections. The 64 reflected half cells have disjoint interiors and fill
a fundamental cell. Therefore the reflected **witness union** has volume
64 times its half-cell union; period mass remains one. This normalization
does not assert that the witness family exhausts a maximal level set.

The independent checker imported no D/E field or geometry routines. It
read the exported affine endpoints, implemented its own Fraction vectors
and root enclosures, and verified all **6,096** source endpoint conditions
at t=0 and t=1, plus **4,572** cell-face/positive-width conditions. Every
source-center difference and half-cell face is affine in t, so validity
at these endpoints implies validity throughout the supplied closed chamber.
It also checked 120 adjacent order gaps, with the single expected boundary
tie. Thus no source geometry or ordering assumption is left pending.

## 2. Exact union and independent polynomial evaluation

The primary program uses rational spacing bounds and lower rational
sixth-root radii. Every retained box is strictly rounded inward onto grid
10^-35. The source endpoint inequalities are checked using spacing upper
and lower bounds with the correct sign for possibly negative source
indices. Thus the rational boxes are subsets of the true witness boxes.
Pruning removes only exactly contained integer boxes. The union routine
slices at all integer endpoints, recursively unions active projections,
and merges intervals in the last coordinate. It never adds overlapping
box volumes independently. Its exact ratio is
R_in=64*integer_union/10^210.

Exactly one replay of `inward_smallstep.py` was executed in this owned
directory. The only source change was ROOT, pointing to the main project's
read-only inputs. It returned 127 candidate / 96 retained boxes and the
same exact rational output as the primary result, apart from timing and
the necessarily different source hash. The complete rational R_in is in
`inward_smallstep_results.json`.

The separate E interface was inspected rather than re-executed. Its
fixed-order slab algorithm keeps the full polynomial degree, uses exact
coefficient convolution, and independently checks one internal static
union; the supplied result is degree one with exact static agreement.
Continuity of finite box-union volume extends the polynomial identity to
the chamber boundary t=1, including the endpoint tie. The source audit
above makes these exported boxes admissible witnesses independently of
the event-finding argument.

This audit independently isolated each positive root in a rational
interval of width 10^-120, checking the endpoint sixth powers with
integers. It evaluated E's two ratio-polynomial coefficient vectors at
t=1 using its own sign-aware rational monomial bounds. All 36 candidate
coefficients agree exactly with D's static answer; the constant vector
agrees with the supplied inherited round007 coefficient vector. The new
root evaluation is contained in **both** D and E's reported candidate
intervals and has width approximately 3.81216185626608e-114. Its lower
endpoint begins

```
2.0673973332673815277191179372217682427866005663975806814546176...
```

This is an exact rational enclosure, not a floating calculation or a new
union evaluation. All exact endpoints and coefficient checks are saved in
`independent_review_results.json`. No raw-arrangement event beyond the
single candidate, optimum claim, or larger block search is needed here.

## 3. Explicit finite L1 and strict-level transfer

The primary inward rounding verifies a positive uniform source-containment
margin. Its minimum margin is 10^-35, and it chooses

\[
\eta=\tfrac12 10^{-35}=5\cdot10^{-36}.
\]

Replace each selected atom by the nonnegative unit-mass bump
(2 eta)^-6 times the indicator of its translate of [-eta,eta]^6. The
strict margin ensures every witness cube contains the **entire** selected
bumps, including their support boundaries. A uniform indicator bump is
already a legitimate L1 function; a smooth-bump assumption is unnecessary.

For ell<=4 the selected relative source indices lie in [-1,2] before
reflection, and in [-2,2] after reflection. The independent checker
verified that the maximum absolute index is exactly two. Hence receiver
cells indexed by [-K,K]^6 require only source atoms indexed by
[-K-2,K+2]^6: collar two is sufficient in every direction.

Define f as the sum of those finitely many translated bumps, with

\[
K=2^{47}=140737488355328.
\]

Positivity gives the exact norm ||f||1=(2K+5)^6, regardless of any bump
overlap. Reflected inward unions within distinct receiver cells have
disjoint interiors and total measure (2K+1)^6 R_in. Each witness cube
contains m unit masses and has volume at most m, because its radius is
the lower rational sixth-root radius. Its average is therefore at least
one.

Let T=2067397333267/10^12 and set

\[
\theta=\frac{R_{in}+T}{2R_{in}}.
\]

The independently verified R_in>T>0 gives **0<theta<1**. Therefore every
selected center belongs to the **strict** superlevel set {Mf>theta},
even when the average is exactly one. The finite weak ratio is at least

\[
Q_{finite}=\theta R_{in}
\left(\frac{2K+1}{2K+5}\right)^6.
\]

The independent Fraction calculation gives

```
Q_finite = 2.067397333267014486960709041367795281764887004252858...
Q_finite - T > 1.4486960709041367795281764887004e-14.
```

The full exact fractions, including the boundary-loss factor and theta,
are saved. Thus this **one finite L1 input** already proves C6>=Q_finite>T;
no asymptotic truncation or heuristic boundary-loss estimate is required.
The function has finitely many terms given by an exact finite sum, even
though physically enumerating its very large index cube is unnecessary.

Because this nonnegative finite L1 function belongs to the general L1
class, it proves the lower bound for that class. Every witness uses a
valid positive radius in the full maximal supremum. For any dilation h>0,
f_h(x)=h^-6 f(x/h) has the same norm, witness radii h*r, threshold
h^-6 theta, and the same weak ratio. Thus the proof respects the
all-positive-scale model; it has no fixed absolute scale restriction.

## 4. Accurate comparisons and acceptance scope

Comparisons use full rational endpoints, never displayed decimals:

* `rounds/008/frozen_lower_before_candidate.json`, the complete accepted
  round007 ideal upper endpoint;
* `coordination/phase_20261025/frozen_thresholds.json`, the unchanged phase
  lower-construction upper endpoint.

The exact independent results satisfy

```
Q_finite > T > complete round007 upper > phase frozen complete upper.
```

The finite-ratio margins over these two complete upper endpoints begin
0.0004731950855014715348413993372743447444741 and
0.0004735072561335768629307575277412300067134, respectively.
Both comparisons passed with positive exact Fraction differences.

The strict target claim, ordinary-L1 transfer, collar, bump support,
normalization, common radius, ideal polynomial enclosure and threshold
comparisons are all accepted. There is no pending gap in this single
candidate's audited lower-bound chain. Upper-bound work and project-level
dual-breakthrough acceptance are separate responsibilities; this audit
changes no ledger and authorizes no new lower batch.

## 5. Provenance, dependencies, execution and stop

Dependencies: local Python **3.15.0**, standard-library Fraction, integer
arithmetic, Decimal only for display, and the explicitly hashed read-only
source/result files listed in `independent_review_results.json`. No
third-party mathematics package is required. All inequalities used as
evidence are integer/Fraction checks; Decimal displays are not proofs.

Execution records:

| execution | actual exit | program elapsed | handle |
|---|---:|---:|---|
| unique primary replay | 0 | 0.3110601902 s | exec chunk `47e370`, completed; PID not emitted |
| independent checker | 0 | 0.6903614590 s | PID 23263, exec chunk `5a06ee`, completed |

The raw primary stdout is `replay_stdout.txt`; `replay_stderr.txt` is empty.
All input hashes and checker source hash are recorded in the checker JSON;
the original/replay source hashes and sole ROOT substitution are recorded
in `registration.json`. Primary original source SHA-256:
`20a563672faed6747fc11831521b1c44f6078da3b50a74caf6191942883d7924`.
Replay source SHA-256:
`4eedde51bb32df612e22ad15f920419ccbbb734b7a71f4039319327661da3b10`.

From this owned directory, reproduction is:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B inward_smallstep.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B independent_review.py
```

The replay deliberately reads ROOT inputs from the main clone; the
independent checker also reads D/E's frozen result files without executing
their search or union programs. All audit writes are confined to this
directory. The registered caps were one candidate, one primary replay,
20 minutes; actual arithmetic took about one second total. All processes
have exited. **Stopped: no further lower computation or new round.**
