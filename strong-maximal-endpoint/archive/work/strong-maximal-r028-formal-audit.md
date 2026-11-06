# r028 formal audit — 2026-11-06

Independent, read-only review of `rounds/r028/report.md`, sections 1--8.

## Verdict: accepted

The exact theorem is correctly scoped to coordinatewise-even physical inputs
and a *single globally fixed* ratio vector.  It is not asserted for arbitrary
physical inputs, point-dependent ratios, or the full strong maximal operator.
The positive-orthant factor cancels in the score.

The fibre Hardy--Littlewood step is valid because the fixed-`x` kernel is
decreasing in the fibre variable.  The two-level global weight has the stated
rearrangement, including the endpoint cases `r=0,1`.  Applying the argument
at `s'>s` before letting `s'` decrease fixes the potential strict-level-set
plateau error.  The stated truncation domination supplies the finiteness
needed to iterate the rearrangement inequality.  The scalar test
`h_x(t)=g*(x exp(-t))` and the `L -> infinity` deep construction preserve
the cap exactly, so no artificial upper truncation of `h` is required.

For the Fourier band, the negative tilt and all signs are consistent:
`f'/f=lambda+O(1/k)`, hence `a>0`, the minimizing change is
`delta=x/(2k)(1+o(1))<0`, and the second-order gain is `p x^2 f/(8 k^2)`.
The report now explicitly retains the required
`O((lambda^2+k^-1) f)` Taylor remainder, which does imply the displayed
relative `O(sqrt(log N/k))` on the full band.  The global shifted derivative
bound accounts for the unbounded exponential component in the centered
Taylor expansion.  The report correctly labels controlled occupation as an
open issue, hence does not overstate a cumulative feedback benefit.

Section 7's 33-point third-derivative maximum is explicitly noncertifying;
the horizon test is explicitly only algebraic.  Both numerical qualifications
are adequate.  I independently ran local MATLAB to confirm
`mu=1.2784645427610735`, `p/8=0.0272264632149750`, and
`p/12=0.0181509754766500`, matching the table/coefficient claims.

## Post-rewrite check: full-line two-ray formulation

The rewrite of sections 1--4 and 8 correctly removes the evenness assumption.
For each fibre the positive and negative ray rearrangements preserve global
equimeasurability and increase the corresponding initial-interval averages.
At level `s`, with ray lengths `a,b`, the displayed measures obey
`H <= rD`, `U+H=(1+r)D`, hence `U-D=rD-H>=0`.  Since `F` is concave,
the equal-length increments satisfy `F(U)-F(D)<=F(rD)-F(H)`, exactly the
required direction.  Thus the full-space rearrangement lemma is sound.

The corollary for finitely many translated anchors is also sound: each
translated fixed-anchor operator is conjugate by a translation and has the
same individual weak constant; union bound gives the stated `K Q_n` bound.
It does not assert a bound for point-dependent anchors or ratios.
