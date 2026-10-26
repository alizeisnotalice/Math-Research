# 115: adjacent positive bands and two lawful whole-box expressions

PREPARED ONLY. Original body32 law, common Cauchy time, unit weights,
orthonormal Lebesgue/Fourier normalization are unchanged. Fix
a=1/64, b=1/32, T=1/2 and the two complete114 blocked-node boxes.
No box shrinking, new point selection or window scan is included.

## Applicability of the existing mathematical interfaces

C016 report.md states its exact fan formula for every0<a<T; the
unchanged engine source full() checks precisely positive a<T and reads
both rational endpoints from argv (lines50–53). Its arctangent/Hermite
formula and six physical derivatives use these input endpoints, and
its native Euler residual uses actual a,T. The original finite software
tests are not misrepresented as new-point independent crosschecks.

A057 theorem.md lines9 and149 allows every positive finite interval
[a,T]. Its box_certificate.py checks a>0,T>a, with no fixed-T condition.
The old C063 scheduler's T=1/2 is a scheduling restriction, not an A057
or engine theorem restriction. This pilot therefore uses a new schema
and direct general-window checker input; it does not forge an old
C063/C068 fixed-window receipt.

## Exact addition of values, gradients and Euler identities

Let L(x)=integral_a^b p_t(x)dt and H(x)=integral_b^T p_t(x)dt. Positivity
and additivity of the ordinary time integral give

    G(x)>=L(x)+H(x)=B_[a,T](x).

The intervals are bounded away from0; their positive-time densities
are smooth with matching finite derivative integrals, as in the
accepted C016/A057 formulas. Hence all six spatial derivatives add:
grad B_[a,T]=grad L+grad H. Sum the complete saved interval endpoints
as rational intervals; no displayed decimal or midpoint is used.

For the unchanged six-dimensional alpha1 model, the two exact identities
are x.grad L+5L+b p_b-a p_a=0 and
x.grad H+5H+T p_T-b p_b=0. Their sum is

    x.grad(L+H)+5(L+H)+T p_T-a p_a=0.

The internal p_b terms cancel as the same actual density, not as rounded
midpoints. Both independently saved p_b enclosures must intersect; their
intersection is recorded. The merged exact-rational Euler interval uses
only the outer endpoint densities and summed complete jets and must
contain0. Each original native Euler is independently checked using
its actual window. No additional high-band or union-band native runs.

## Separate-band versus merged-jet box certificates

For either positive band F, the accepted frame32 A057 expression on
each of32 labelled corners d is

 C_F(d)=f_d exp(b_d/f_d-|d|²/(32 a_F²)),
 b_d=sum_i min(g_i^- d_i,g_i^+ d_i),
 f_d=min(U_F,max(L_F,b_d)).

Its whole-box lower bound ell_F=min_d C_F(d) holds on the entire saved
box, including its closed order-face boundary. Thus

    G(x)>=ell_L+ell_H on that complete box.               (1)

Independently, the summed interval jet is a valid jet for B_[a,T], so
applying the same A057 formula with lower cutoffa=1/64 gives another
lawful whole-box lower bound ell_merged. Both bounds are valid even
when an individual band center is below1/65000; A057 requires a strictly
positive center lower endpoint, not a center-target pass. The two old
high-band centers are already below target, but their positive jets
remain valid inputs for the same lower-certificate formula.

The pilot computes three certificates per point: L, saved H and the
summed jet, for192 labelled corner rows total. It compares the two
outward intervals for ell_L+ell_H and ell_merged. A disjoint ordered
interval proves an ordering of these certificate expressions only;
neither expression is an upper bound on true G. If a lower bound passes
target it applies only to that complete saved local box. It is not a
full-root or globalg certificate, and does not imply a new C6 bound.

The merged jet has explicit sum provenance, not fake single-native
metadata. A057 accepts its matching model/window/center/gradient/source
fields; the C063 exact saved-corner adapter is called directly for
formula/geometry verification without pretending a merged jet was one
C016 native output. Native low-band validation uses the same accepted
512-bit physics, replacing only hard-coded Euler a,T with the actual
expected window. All old114 sources and high native outputs stay read-only.

Both114 originals have corners_checked=0 and no high box certificate:
their strictly-below-target centers stopped before that checker stage.
Thus all192 corner rows here are newly generated:64 for the first
high-from-saved-jet certificates and128 for low/merged certificates.
Only the two complete high native jets are reused; no old high corner
rows are claimed or fabricated. The original verification hashes are
explicitly bound in registration.json.
