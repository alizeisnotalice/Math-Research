# C017: an entire time interval, but not yet the whole time family

At the exact B015 point

    x=(100000000000000,99999998999999,99999997999988,
       90268521812170,53714088833376,157107163435)/10^14,

the existing T0=6/125 certificate implies the strict continuous bound

    K_T(x)/T < 1/4100
       for every T in [223/5000,251/5000]=[.0446,.0502].

This excludes a below-64 bound from any fixed time in that interval.
It also excludes every probability mixture of normalized short potentials
whose time support lies in the interval: its cube-domination model cost
is at least

    4100/64=1025/16=64.0625.

This is a model obstruction. It does not lower-bound the actual cube
maximal constant, establish an all-T obstruction, or exclude probability
mixtures with time support outside this interval. It makes no inference
about other frames, beta values, or the possibility of improving C6 by
another construction.

## The continuous certificate

The B015 exact dyadic intervals give, for display only,

    K(T0)/T0 <= .00023228107620985218,
    p_T0(x)  <= .0002124475367841768.

The rigorous comparisons use their full rational endpoints. The negative
q'(T0) interval is recorded but is not used to assume a single maximum.

For clarity, the time comparison is a positive-measure statement before
it is a density statement. On R32, for s>=t each Cauchy factor satisfies
c_s(u)<= (s/t)c_t(u). Hence their product measures satisfy

    product_32 c_s du <= (s/t)^32 product_32 c_t du.

Push both nonnegative measures through the same body32 linear map A.
The resulting measures obey the same inequality. Their fixed-time
densities are continuous, so p_s(x)<=(s/t)^32 p_t(x) holds at every x.
Reversing the two times gives p_s(x)>=(s/t)^32 p_t(x) for s<=t.
Equivalently, solving six basis variables in the positive coarea formula
retains all 32 Cauchy factors and the same determinant factor. The
exponent is 32, not 6; there is no conditioning on a zero-probability
fiber and no omitted 26-factor contribution.

Integrating these comparisons on either side of T0 gives the two upper
envelopes in `theorem.md`. Dividing by T yields A/T+B T^32. The exact
check proves A>0 and B>=0, so each envelope is convex and its maximum
on its interval occurs at an endpoint. Their outer endpoint upper values
are approximately

    left:  .00024367340944126832,
    right: .00024295653734983247,

both strictly below 1/4100. Together with the center bound this proves
the entire interval, including both endpoints. No intervening time was
evaluated, and the displayed decimals are not used as proof inputs.

The mixture conclusion is valid for essential minima because the point
is off S2 and normalized K_T is jointly continuous locally in space and
on this compact time interval. Every supported probability mixture is
therefore continuous near x and has essential cube minimum at most its
value at x. The theorem also supplies the necessary domination argument
for arbitrary time supports if a genuine uniform all-T bound is proved
later; that uniform bound is currently missing.

## What was actually run

One registered Python Fraction postprocessing batch performed **11 named
exact scalar checks**, within the limit of 12. It decoded the existing
native dyadic K and p endpoints, checked them against B015's saved exact
fractions and hashes, and checked the fixed rational envelope endpoints.
The batch took 0.000848042 seconds, PID 34457, exit 0.

New kernel calls: 0. New G computations: 0. New time evaluations: 0.
No fan expansion, common denominator, root isolation, numerical search,
MATLAB invocation, or subagent was used. All previous directories stayed
read-only. The exact endpoint bounds and every named check appear in
`scalar_certificate.json`, with the input and source hashes fixed by
`input_manifest.json` before computation.

## The remaining all-T issue and a finite acceptance route

Neither q'(T0)<0 nor the previous failures at T=.032,.048,.064 proves
that q has only one peak. The two continuous power envelopes do not
control T arbitrarily close to zero or infinity. A new exact G value
would settle only the large-time tail through q(T)<=G/T.

The full exact cone formula makes p_T(x) a rational function of T^2
divided by 64 pi^6. `theorem.md` proves a complete finite all-T acceptance
criterion: choose a certified rational lower bound ell for pi, isolate
every positive root of the exact rational equation

    R(T^2)=64 ell^6/4096,

and certify c a-(u_pi/ell)^6 K_upper(b)>0 on every critical-time bracket
[a,b]. Exact derivative signs between the complete root brackets, the
zero limit at T=0, and finite G then close all T. This handles arbitrary
root multiplicities and does not presume unimodality.

That certificate has **not** been computed. The old native output already
reports 16931 cached node values, so naive common-denominator expansion
and an unrestricted Sturm calculation are not justified within a small
budget. The theorem gives a concrete future cap: one exact G pass;
exact rational aggregation with a 1200-second/1-GiB limit and reduced
degree cap 512; at most four positive critical-root brackets; and at
most four new K endpoint passes, with explicit failure and stop rules.
The algebraic-compression gate is a serious possible failure, not an
expected PASS. It should be checked for feasibility before paying for
new G/K work unless the G value is independently wanted. None of that
future work is released or run here.

## Reproduction and retained files

Run `python3 check_scalar.py` from this directory to reproduce the
registered exact scalar postprocessing; it never executes a kernel.
The proof, registration, scalar checker, input manifest, full scalar
certificate, this report and delivery receipt are the necessary retained
artifacts. There are no scratch files or binaries. Nothing was staged,
committed, or pushed. The current established C6 upper bound remains 64.
