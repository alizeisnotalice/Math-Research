# C095 smooth-kernel correction: explicit global source domination

Author conclusion: the stated body32 six-dimensional construction has a complete
analytical proof of a small positive admissible amplitude. Independent review
is pending. It strictly reduces the mass of the same normalized positive-part
kernel, preserves cube domination and the whole-space positive probability
source, and does **not** establish a new numerical C6 bound or S>0.

Only upper/round095/smooth_kernel_correction was written. No body/native kernel,
source quadrature, old grid, old program, MATLAB, compilation or recursive agent
was run. No stage/commit/push occurred; main retains proof acceptance, final
cost judgment and public archival responsibility.

## Outcome and what is new

For the original unit-weight32-direction alpha1 Cauchy generator A, take any
rational certified cube floor g0>0 and c=g0/6. Let K0=(G-c)+ and use the explicit
even C-infinity product bump supported in two boxes of halfwidth1/100 around
±(6/5,1,4/5,3/5,2/5,1/5). Those boxes are outside the cube, inside{G>c}, and a
proved positive distance from true S3. Theorem.md gives explicit rational
constants and integers, including

    C0=1/c, C1=16LM/c, M=145920043,
    epsilon=c/[32L(M+1)] >0,
    b<=C0 K0, (Ab)+ dx<=C1 nu0 on all R6.

Consequently K0-epsilon b is nonnegative, even, L1 and unchanged on the cube;
its source is nu0-epsilon Ab dx-delta0, with a positive probability measure.
The unit negative atom and the old singular/interface source are fully retained.

The proof pays the substantive conditions rather than assuming them:

1. A positive coarea event argument forces three large Cauchy parameters, one
   at least|x|/(32sqrt6). It gives p_t(x)<=Abar/(R²delta⁴) at fixed distance
   delta from S2. Combining with Ct^-6 gives an explicit uniform Green tail.
2. The whole jump cylinders from the bump stay at leastdelta from S2, because
   the bump stays away from all three-direction spans. Their intersection with
   Omega={G>c} is therefore bounded by an explicitR.
3. On that bounded bad set, a single fixed body jump direction has a common
   far segment where G<=c/2. This proves the genuine source density lower bound
   c/(16L), without asserting nu0 has positive density on axes or S2.
4. On the closed complement of Omega, the positive probability approximants
   [P_tK0-K0+U_t]/t dx and the correctly directed upper-semicontinuous
   Portmanteau inequality prove measure domination, including the level set.
5. Explicit derivative bounds give ||Ab||infinity<M, while Ab is continuous
   L1 with zero integral. The stated epsilon then pays both nonnegativity and
   the entire signed source. No new truncation interface is introduced.

The original finite-source lemma §C was read in full and is stated for arbitrary
nonnegative L1 kernels with a whole-space finite positive source, not justR_c.
The dimension6 obstacle proof was also read. These supply the required general
L1/all-positive-common-radius transfer for the corrected kernel.

## Prior-work identity and limits of novelty

The08492-source registry and its scope/schema were read at the numerics clone;
its finite scope is respected.081's exact body-plus-axes duplicate seal was
checked. Killed-Laplacian cube, body-plus-axes,152-direction sparse4 and alpha3/2
routes have different generators or targets and supply no values/constants here.

The broad spatial-bump subtraction idea is old. The finite-dim U10 source already
gave three-dimensional joint-bump conditions, exterior payment and a finite
tail reduction, while leaving the bounded source region unverified. Old009 and
010c strict negative-source witnesses excluded particular larger amplitudes,
not all small perturbations. Those originals were read, not rerun. C095's
specific contribution is the six-dimensional off-S3 box geometry and explicit
global finiteC1, closing the source condition for a small amplitude. No
literature-wide originality claim is made. Full lookup fields and source hashes
are in prior_lookup.md and input_manifest.json.

## Mass and final-cost accounting

With N0=64(g0-c)=320c, the normalized mass decreases by
epsilon integral(b)/N0. It is at least

    Delta_min = [2/(100^6*9^6)]/[10240L(M+1)] >0.

Thus for a matching certified I_U, the sufficient cost slack gains at least
(g0/6)^(6/5) Delta_min. This does not prove the total slack positive. No I_U was
computed, and a reduction of an auxiliary bound is not an actual C6 improvement.

The displayed numerical-free existence instance uses the analytically proved

    gstar=1/713671068669617786572320538624000000000000000000000000000.

It comes from a positive fixed-time coarea lower bound on the cube, using the
already accepted determinant160 basis; no matrix or kernel is reevaluated. For
this intentionally weak input, the exact normalized guaranteed gain is

    1/93002857912729530577057877253245891935562180772560063112586748784117163718932807680000000000000.

All exact T,R,L,C0,C1,epsilon and tail endpoints are stored in results.json.
The bump-volume upper bound integral(b)<=2(2/100)^6 also bounds the actual
gain for this chosen epsilon by34012224 times Delta_min. Thus this explicit
existence instance is extremely small in cost effect; it is not a competitive
budget proposal. Better quantitative constants/amplitudes would need separate
work and certification. The proof itself can accept any later certified
same-model rationalg0, without mixing a new generator or stale kernel values.

## One bounded exact run

Registration was frozen before execution:
ec4c7d4ac91ce8c33c0354026abddc38f0e3d3356257bebbc467233929912721.
Input manifest:
5005f180c3727c14a51a65d63d2e44a3baadc8ba0b0a19589eb634d493839ad7.
Theorem source:
d8ff73cbc86ae25ce827308185c86b8423f55390bb6e196f77eda09213e84c16.

Three new groups ran exactly once, all PASS:

| Group | Actual verification |
|---|---|
|1|Exact rational outside-cube, absolute-gap/S3 distance and radial margins; zero plane/triple re-enumerations.|
|2|One SymPy bump first/second derivative identity and the explicit Hessian/Ab constants.|
|3|One complete rationalgstar parameter graph: short/long tail, line segment, source bound, amplitude and normalized mass gain inequalities.|

Local venv Python used SymPy1.14.0. PID14225, tool chunkbc140c, exit0. Saved
program elapsed0.331476083s and CPU0.252034s; final stdout was printed slightly
later. The tool's separate reported wall field was0.268773s and is retained as
such, not substituted for the program's monotonic measurement. CPU limit30s and
default terminating alarm179.8s covered cleanup/output; no child processes.
Zero native/body/Arb/source-integral calls, zero old replays, zero retries.
These finite checks support constants; they do not replace the analytical
all-space proof or an independent review.

## Handoff

The source, proof, prior lookup, fixed registration, pins, one-shot run marker,
exact results and execution receipt are sealed. All source pins remained
unchanged after the run. Main should arrange independent review of the full
coarea event coverage, cylinder geometry, density/Portmanteau measure argument
and finite-source transfer before accepting the theorem. There is no automatic
next batch, native authorization or change to the global C6 ledger.
