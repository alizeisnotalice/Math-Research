# Independent final audit: dilation obstruction and budget bridge

## Verdict

**PASS, within the explicitly qualified scope.**  The dilation note gives a
valid finite-depth obstruction to the proposed global angular-mean baseline;
it is not presented as a counterexample to the strong maximal inequality or
as an unconditional near-\(Q_N\) construction.  The budget bridge has the
correct right-slope rigidity, global source-measure identity, and
source/receiver Jacobian bookkeeping.  Its remaining cone-packing and
active-region requirements are stated as missing hypotheses rather than
claimed consequences.

## Dilation and the explicit packet

For \(f_D(x)=f(e^Dx)\), rectangles map bijectively under \(x\mapsto e^Dx\),
so both the strict receiving halo at any fixed threshold and the
\(\Phi_N\)-modular scale by \(e^{-ND}\).  Thus the endpoint
volume-to-modular ratio and the strict witness margin are exactly
preserved.

For a compact logarithmic box, the displayed formula
\[
\eta_i-\frac1{b_j}
 =\frac{b_jv_i-S_j}{b_j(b_jD+S_j)}
\]
does put the source angular support in a ball of radius \(O(D^{-1})\) in
each fixed simplex.  Its product Dirichlet mass is consequently
\(O(D^{-(N-r)})\), uniformly over the finite-depth radial coordinate.
This proves that the all-label/all-angle radial mean tends uniformly to
zero under dilation.  No limiting radial-terminal substitution is used.

The concrete cube calculation is also correct.  In one block the
Dirichlet angular probability of the dilated log box is
\[
\frac{\Gamma(b)f_{\rm IH,b}(S)}{(bD+S)^{b-1}}.
\]
Here \(f_{\rm IH,b}\le1\) and \(\Gamma(b)\le b^{b-1}\), giving the product
bound
\[
B_D(w)\le a_ND^{-(N-r)}.
\]
For \(D=2\), \(N>r\), and \(a_N=1+N^{-2}\), this is at most \(5/8<1\)
(the extremal case is \(N=2,r=1\)).  Hence the radialized terminal is
everywhere subunit and has no strict receiving halo, while the original
orthant cubes have strict margin \(N^{-2}\).

Angular averaging preserves the ordinary integral with the common
orthant factor cancelling.  Because \(B_D<1\), its modular is its ordinary
integral, so
\[
\frac{{\cal J}(f_{N,D})}{{\cal J}(B_D)}
 =\frac{\Phi_N(a_N)}{a_N}
 =(1+\log a_N)^{N-1}\longrightarrow1.
\]
This validly refutes a mechanism based only on sharp Jensen budget,
near-unit amplitude, and a forced global active baseline.  It does not
assert that this explicit cube is near-\(Q_N\); the separate near-endpoint
dilation paragraph is correctly conditional on a strict compact
near-endpoint packet/exhaustion result.

## Budget bridge

The source domain is correctly the global
\(\Omega_L=\{w:L\alpha_j+w_j>0\}\), not the local increment orthant.
With the stated unweighted coarea measure and uniform label/Dirichlet
averaging, the affine Bregman term averages to zero fibrewise.  Thus the
identity
\[
{\cal J}_L(F)-{\cal J}_L(B)
=C_L\int_{\Omega_L}e^{-\sum w_j}\prod_jv_j^{d_j}
 {\cal D}(w)\,dw
\]
has the right source measure and contains no hidden strip restriction.

For \(N\ge2\), the right derivative of \(\Phi_N\) supplies the claimed
quadratic lower Bregman bound on \(B\in[1,A]\).  Near \(B\), its curvature
is order \(N^2\Phi_N(B)\); away from that \(1/N\) neighborhood, monotonicity
of the right slope gives a linear tangent gap, which dominates the required
\(N\Phi_N(B)(x-B)^2\) on the fixed cap interval.  The linear branch
\(x\le1\), including \(B=1\), has the same bound.  Therefore the
\(1/(N\Phi_N(B))\) defect-to-\(L^2\) estimate is sound.

The finite-\(V\) response factor
\[
\prod_j\left(\frac{V_j+s_j}{V_j}\right)^{d_j/2}
\]
is the exact change-of-coordinates factor.  Squaring it and multiplying by
the receiver Jacobian \(\prod V_j^{d_j}\) gives exactly the source
Jacobian \(\prod(V_j+s_j)^{d_j}\).  The common baseline is preserved by
each conditional Markov kernel; only \(F-B\) is zero-extended.  Thus no
unweighted-simplex or label-loss assertion is being smuggled into the
local estimate.

The lower-core inequality is a source estimate: on
\(\sum w_j\le q_0-KL\), one has \(F\le\Phi_N(F)\) and
\(e^{\sum w_j}\le e^{q_0-KL}\).  The \(O(L^2/N)\) critical-strip
probability follows from width \(O(L)\) times the critical density
\(f_{N,\mu}(q_0+O(L))=O(L/N)\).  Both are correctly marked as insufficient
without source-to-receiver transport/localization.

Finally, the cone identity is an exact Fubini change of variables
\(w=u+s\).  Its multiplicity \(M_{\cal U}(w)\) is precisely the
uncontrolled upper-cone reuse factor.  A global defect budget bounds only
\(\int d\Lambda_{\rm def}\), not
\(\int M_{\cal U}\,d\Lambda_{\rm def}\).  The proposed packing inequality
is explicitly sufficient but unproved.  Accordingly, neither note makes
a hidden strip, packing, or global-\(C_N\) claim.

No numerical conclusion was needed for this audit.
