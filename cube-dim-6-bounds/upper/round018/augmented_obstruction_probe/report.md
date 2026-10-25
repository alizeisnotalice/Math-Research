# B018: the fixed augmented candidate fails at the B015 point

**One authorized probe completed and stopped.** At the exact existing B015
point z, the augmented38-column Cauchy frame with beta=1 and T=2/45 satisfies

\[
\frac{10048510817}{10^{15}}
<K_{1,2/45}(z)
<\frac{100485108171}{10^{16}}
<\frac1{92160}=\frac{T}{4096}
<\frac1{90000}=\frac{T}{4000}.
\]

Consequently this fixed frame/cutoff cannot achieve the short-potential
strict C6<64 route, and also fails the62.5 target. The model-cost lower bound
from this point is69.1091901168893571..., **not a C6 bound**. No conclusion is
made about other beta values, cutoffs, mixtures, or the global minimizer. The
previously established operator upper bound64 remains unchanged.

## Fixed input and new-frame identity

The point was copied exactly from B015's validated second native input:

```
z = [100000000000000,99999998999999,99999997999988,
     90268521812170,53714088833376,157107163435] / 100000000000000.
```

It was not rounded, moved, perturbed or optimized. The single new evaluation
used the model

\[
\psi_1(\xi)=\sum_{v=(1,\pm1,\ldots,\pm1)}|v\cdot\xi|
             +\sum_{i=1}^6|\xi_i|,
\quad T=2/45.
\]

All38 direction weights equal1. The new model descriptor includes its
directions, weights, symbol, Fourier normalization, Lebesgue measure and
augmented fan SHA. Its independent frame SHA is
`f256ca34ba7be9417cae5d139a0e5387be1d844ae0c22ea1c82fa62b4f19f52a`.
No base-frame kernel value or old C014 numerical result was inherited. The
two earlier C012 points at this same augmented frame/cutoff were checked to
be exactly different from z; their survival did not predict the present
necessary test.

The main agent explicitly released this single probe in its assignment;
the complete authorization is preserved verbatim in `registration.md`.
No separate release file was requested or used. The batch allowed one
complete native attempt,120 seconds for its process,1200 seconds total,
and no retry, G/H pass, replacement point or second cutoff.

## Exact fan and native-domain checks before the call

The new fan was read from the C012 augmented-frame directory, with SHA
`450caf04d124334474b91f3d9fdabce4eadf3259ea47def9094ae5eb056b6c4e`.
All21 cones and126 ray occurrences were checked exactly:

- rays and stored determinants equal the audited base fan;
- each determinant was independently recomputed from its six integer rays;
- each old q equals the explicit32-column symbol at its ray;
- each new q equals old q+sum_i|r_i| and independently the explicit38-column symbol;
- rays lie in the nonnegative ordered chamber, and all38 linear forms have
  constant weak sign on every cone.

The previously certified geometric cover is unchanged. The full signed
permutation symmetry is valid for both the body32 collection and the six
unit coordinate axes. C014 reads all new q values from this fan and rebuilds
both nodes r.x/q and weights det/prod(q); its source and binary were not
modified. Thus no base-frame normalization weight is silently retained.

From the18 unique rays, all signed/permuted primitive zero forms were
generated, modulo global sign. All9308 forms are nonzero at z, with minimum
absolute integer numerator10. The conservative raw numerator bound is1.2e15,
the largest raw node denominator is1.7e16, and the comparator product bound
is2.04e31. These are inside signed int64/int128 limits. There was no exact
zero node, unsupported domain, precision retry, skipped term or input repair.

## Complete interval result and pure postprocessing

The unmodified C014 native512 engine returned all967680 logical terms,
483840 computed sign pairs, zero_node_terms=0 and confluent_terms=0. The
complete K, six gradient, p_T and Euler dyadics are retained in
`native.stdout.json`. The K display midpoint is

```
0.0000100485108170111742742895409032633580892414417577774...
```

Its exact dyadic width is below2.81e-120. Its upper endpoint is below the64
threshold by more than8.02183e-7. The opening rational enclosure and both
threshold comparisons were reconstructed with Python Fraction arithmetic,
without using a displayed floating midpoint as evidence.

`check_results.py` separately decoded the saved dyadics, checked every source
and output SHA, complete group counts, positive K and p_T, all ordered finite
intervals and their accuracy, and the native Euler interval. It independently
formed the endpoint interval of

\[
z\cdot\nabla K_T+5K_T+T p_T,
\]

which contains0. It did not import the kernel driver or execute any binary.
Only saved data were used for the supplementary scalar diagnostics

\[
K_T/T\approx0.0002260914933827514212,
\qquad (T p_T-K_T)/T^2\approx-0.0004161843593448660790.
\]

The latter interval is strictly negative at this single T; it is not a
time-interval exclusion or a claim about an optimal cutoff. Complete rational
intervals and all six derivative intervals are in `checked_results.json`.

## True augmented S2 and the essential-infimum consequence

The new frame's S2 is the union of pair spans among all38 columns. The checker
examined all703 pairs:496 body/body,192 body/axis and15 axis/axis. No old S2
membership or distance certificate was substituted for this calculation.

For a pair v,w, set a=|v|^2, b=|w|^2, c=v.w, U=z.v, V=z.w. Every pair is
independent, and its squared distance is exactly

\[
|z|^2-\frac{bU^2-2cUV+aV^2}{ab-c^2}.
\]

The minimum across the703 pairs is

\[
\operatorname{dist}(z,S_{2,\mathrm{aug}})^2
=\frac{6020805515960314567115330493}
       {40000000000000000000000000000}>0.
\]

The closest pair has indices28,31 in the explicitly recorded frame order;
its distance is approximately0.38796924865124. That pair happens to be a
body/body pair, but this was obtained by including all axis pairs as well.

The C012/F012 theory proves that the true augmented short potential is finite
and continuous off this S2. Since z is off S2 and K(z)<T/4096 strictly, some
neighborhood of z also lies below T/4096. Its intersection with the cube has
positive six-dimensional Lebesgue measure, even though z1=1. Therefore the
cube essential infimum is at most K(z), so the strict point obstruction
cannot be removed by changing a null set. The same argument excludes the
larger target1/90000. There is no convexity or coordinate-monotonicity premise.

The augmented positive Markov semigroup retains mass||K_T||1=T and the
allscale/general-L1 cube transfer C6<=T/(64m) when a global essential m is
available. This theorem is the interface being tested, not a new global m
derived from the single point. Its fixed-model cost must instead be at least
T/(64K(z))>69.10919.

## Provenance, execution and stopped state

Pinned engine source SHA:
`247a1f6c6967ad52f34320b999818214488700651b46cb8e5feb3eac7064953b`.
Pinned binary SHA:
`67119aa799041ed75b66f043444640a2a62556f58517cbdf3bd6f7cccd564f04`.
The engine was mathematically reviewed by the main agent and F014 within the
generic validated-fan domain. Both the C012/F012 frame theory and C014/F014
engine reports were read and their source hashes retained in the input
manifest. Runtime FLINT/GMP/MPFR library hashes were rechecked before launch.

Preparation: PID34305, tool chunkff3ac9, exit0, about0.0403s for exact checks.
The only native driver: PID34326, tool session49143, start chunkc1c53c and
completion chunk140a5b, exit0. Native child PID34334 exited0 after about3.910s
outer wall time. The child was waited and the session is closed. The registered
input manifest SHA is
`2d4739925ce59677d4b88830c2800c6676861fe97694446b05ec75b1e0e78861`.
The complete native stdout SHA is
`f2561cf034a25e946fffe5a1bca22ae4115e6ac94cf9b3a09a0e0e1e5b334336`.

Pure postprocessing: PID34390, tool chunkc9ba64, exit0, about0.0145s, zero
additional kernel calls. Python3.15.0, SymPy1.14.0 and python-flint0.9.0 were
actually inspected; the native ABI requires FLINT3.6.0 and512-bit precision.
Full paths, input hashes, exact fan audit rows, raw output, exit receipts and
the file allowlist are preserved in this directory. No B015 or other agent
file was modified, and no MATLAB, subagent, stage, commit or push was used.

B018 is **complete and stopped after its only allowed probe**. Any other
beta, cutoff, point, mixture or frame experiment needs a separate assignment.
