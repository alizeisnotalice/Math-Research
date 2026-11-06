# r028 audit: two-ray rearrangement removes the even-input restriction

Date: 2026-11-06.  The proposed two-channel argument is valid and upgrades
the fixed-global-ratio anchored equality from coordinate-even inputs to all
nonnegative inputs (within the same reflected anchored operator).

For one coordinate, write the positive and negative ray fibres as \(f_+\)
and \(f_-\), and rearrange them separately decreasingly.  This preserves the
full-space distribution.  At a positive receiver radius \(x\), the average
is
\[
 \frac{\int_0^x f_+(t)dt+\int_0^{rx}f_-(t)dt}{(1+r)x},
\]
and at a negative radius the two rays exchange roles.  Separate
Hardy--Littlewood therefore increases both output rays pointwise; the outputs
are decreasing in their radii.

At a strict level, let the two ray halo lengths in a fixed transverse fibre
be \(a_z,b_z\).  After multiplying the two interior-point inequalities by
their lengths and integrating in \(z\), the unnormalised global weight has
high-level measure
\[
 H=\int[\min(a_z,rb_z)+\min(b_z,ra_z)]dz
\]
and support measure
\[
 U=\int[\max(a_z,rb_z)+\max(b_z,ra_z)]dz.
\]
With \(D=\int(a_z+b_z)dz\), exact algebra gives
\[
 U+H=(1+r)D,\qquad H\le rD\le D\le U.
\]
Wolfram 15 symbolically verified all three identities/inequalities for
\(a,b\ge0,0\le r\le1\).  Hence its global decreasing rearrangement has
weight \(2/(1+r)\) on measure \(H\), then \(1/(1+r)\) through measure \(U\).
If \(F(q)=\int_0^qf^*(y)dy\), global Hardy--Littlewood yields
\[
 sD\le\frac{F(U)+F(H)}{1+r}.
\]
Since \(F\) is concave and \((U,H)\) is more spread than
\((D,rD)\) with the same sum,
\[
 F(U)+F(H)\le F(D)+F(rD).
\]
Equivalently, use \(U-D=rD-H\) and monotonicity of \(f^*\) to compare the
two equal-length increments.  This is exactly \(sD\le DT_rf^*(D)\).

The previous strict \(s'>s\), then \(s'\downarrow s\), repair applies
unchanged.  At \(r=0\), \((H,U)=(0,D)\); at \(r=1\), the concavity step is
\(F(U)+F(H)\le2F(D)\).  Thus both endpoints are covered.

The same compact-unused-fibre induction works after tensoring: each used
physical coordinate has a two-ray \(\min(1,C/|x_i|)\) tail, while the next
unused coordinate retains compact integrable fibres.  Consequently the
global fixed reflected-ratio tensor operator has the same scalar
distribution domination and exact cap constant \(Q(n,{\bf r},A)\).  The
deep hyperbolic lower input is even and hence supplies the reverse inequality.
This remains a fixed-global-ratio reflected anchored result, not an adaptive
ratio or full strong maximal theorem.
