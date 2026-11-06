# r028 two-channel full-line extension audit — 2026-11-06

**Accepted.**  The separate positive/negative ray rearrangement removes the
coordinate-even hypothesis for the fixed-reflect-anchor operator.

For a fibre with strict positive/negative halo lengths `a,b`, set `D=a+b`.
The combined level-test weight has high-level measure

`H=min(a,r*b)+min(b,r*a) <= r*D`

and support measure

`U=max(a,r*b)+max(b,r*a)`.

Exactly, `U+H=(1+r)D`, hence `U>=D`.  If
`F(q)=int_0^q f*(u)du`, global Hardy--Littlewood gives

`sD <= (F(U)+F(H))/(1+r) <= (F(D)+F(rD))/(1+r)`.

The latter inequality follows from concavity: the pair `(U,H)`, at fixed
sum, is at least as imbalanced as `(D,rD)`.  This includes `r=0,1`.
Each ray output is decreasing after the independent ray rearrangements,
and their direct-integral rearrangement remains globally equimeasurable with
the original input.  The strict-level `s'>s` repair and coordinate iteration
therefore carry over verbatim.

This proves the finite-cap fixed-vector equality for arbitrary nonnegative
full-space inputs of the reflected anchored family; it remains unrelated to
point-dependent ratios or the full strong maximal operator.  Local MATLAB
random stress test of the `H,U` identities and concave-pair inequality found
only floating-point noise (maximum `8.88e-16`).
