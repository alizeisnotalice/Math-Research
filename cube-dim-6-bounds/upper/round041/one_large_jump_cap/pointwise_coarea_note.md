# Pointwise version of the line integral in theorem section1

The identity

    integral_R r_t(x-a v0) da
       =q_(v0perp),t(P_(v0perp)x)/sqrt(6)

first follows as an almost-everywhere identity from ordinary coarea and
the pushforward measure. Here is an explicit justification at every x,
including the singular directions relevant to this task.

Fix t>0. The k=0 version of the independently proved positive inequality
(5) applies to the31-parameter density r_t. All its31 one-direction
quotient laws have coercive symbols: on the corresponding quotient their
frame Gram is32I-(P_Wv0)(P_Wv0)^T>=26I. Their densities therefore have finite
peaks by Fourier inversion. Taking these finite peaks in(5) yields

    r_t(z)<=C(t)/(|z|^2+992t^2)

for a finite constant C(t). No numerical value of C(t) is used or computed
in this note or in the cap budget. The density inequality holds everywhere
by continuity, as proved for(5).

For x in any fixed compact set, the last bound at z=x-a v0 is dominated
by an integrable multiple of(1+a^2)^-1. Consequently the line integral is
a continuous function of x, by dominated convergence. The quotient density
on the right is also continuous. Their a.e. equality thus extends to every
x. This proves the pointwise form used in theorem(2) without assigning a
density to a zero-probability fiber. It introduces no new parameter choice,
scalar configuration or kernel evaluation.
