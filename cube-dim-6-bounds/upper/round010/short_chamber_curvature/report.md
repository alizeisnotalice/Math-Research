# A010: strict analytic disproof of short-potential chamber convexity

For the fixed body32 Cauchy model, the short potential
K_T=integral_0^T p_t is **not convex** on

    x1=1>x2>x3>x4>x5>x6>0

for any T>0. In particular this rules out the proposed whole-chamber
convexity shortcut at T=4/125. There must be an interior point with a
strictly negative second derivative in one of the four directions
e2-e3, e3-e4, e4-e5, e5-e6. A sufficiently short symmetric segment at
that point is a strict midpoint counterexample.

This is an analytic existence proof, not an evaluated numerical triple.
The proof is self-contained in `theorem.md`; the final review and
integration remain the parent agent's responsibility. No C6 endpoint
is changed and no global essential lower bound m_T is certified.

The mechanism has four parts:

1. Away from S3, at least four independent Cauchy parameters must be
   bounded away from zero. Solving a basis containing those parameters
   gives p_t=O(t^2) locally uniformly. The actual Cauchy first and second
   derivative bounds then prove K_T is C2 there.
2. A span of three sign vectors has at most four distinct coordinate
   absolute values. Strict ordered points have six; averaging just one
   adjacent pair leaves five. Thus all points and walls differentiated
   in the argument are off S3.
3. At c=(1,1/2,1/2,1/2,1/2,1/2), an explicit positive density estimate
   gives K_T(c)=infinity. Fatou's lemma gives divergence toward c from
   within the chamber.
4. On the explicit zero-sum spread ray with
   z=(0,51,22,2,-24,-51)/200, divergence at s=0 and finiteness at s=1
   force a negative derivative somewhere. Summation by parts forces a
   negative adjacent-coordinate gradient difference. Evenness at the
   corresponding smooth permutation wall then forces a negative second
   derivative at an interior point.

The rank-two divergence is independently proved here using a positive
six-variable density formula. It is consistent with A009's stronger
S2 barrier but does not depend on its scalar certificates. The positive
Jensen gap on B009's three points does not conflict with this result;
it only tested one segment.

## Actual finite verification

The registered exact Fraction/FLINT calculation completed successfully
in 0.001663542 seconds, PID 25147 as recorded by the output artifact.
It checked the 32 directions, AA^T=32I, the basis inverse, exact rational
path inequalities, summation-by-parts coefficients, opposite sign-pair
count, and scalar derivative polynomial identities. Its key output is

    det B=-32, C=49, t0=1/49,
    p_t(c) >= t^-2 / [1787153625614625485619200 pi^6]
                 for 0<t<=1/49.

This arithmetic supports the proof; the rank argument, dominated
differentiation, Fatou step, and mean-value arguments are analytic
statements in `theorem.md`, not claims established by a sampled test.

New G/H/p kernel passes were all zero. No finite-difference search,
Monte Carlo, spatial grid, MATLAB session, or subagent was used.

## Reproduction and necessary files

From this directory run:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python verify_structure.py

`registration.md` fixes the scope and 1200-second cap.
`input_manifest.json` locks exact inputs and SHA-256 hashes before the
check. `structural_certificate.json` records the actual exact output,
runtime, executable, and source/manifest hashes. The verifier uses
Python 3.15.0 and python-flint 0.9.0 in the recorded local environment.
The source hashes include the prior G counterexample, B009 short-potential
values, A009 barrier, and prior Cauchy derivative discussion for provenance;
the new proof does not call their evaluators or depend on their results.

All six retained files are necessary report/proof/registration/code/input/
certificate material. No scratch output or compiled binary was created.
Nothing was staged, committed, or pushed. This route is closed; there
is no reason to spend the registered point budget on finding a numerical
counterexample unless a later, separately authorized task specifically
needs explicit numerical coordinates.
