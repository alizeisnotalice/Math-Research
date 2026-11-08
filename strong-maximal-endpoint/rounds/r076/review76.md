# r076 independent review: near-corner charge and entropy chain

## Result

The three notes are mathematically consistent after one boundary-branch
repair in nearcorner76. They establish a finite post-sharing fixed-corner
Hardy module and conditional converses for a full witness cover. They do not
construct a continuous moving-anchor cover or prove the general
strong-maximal bound.

## Near-corner optimization

The containment implication and threshold direction are correct:

\[
 R_b(x)\subset R_a(x),\quad {|R_a|\over|R_b|}\le e^\tau,
 \quad H_bF>1
 \quad\Longrightarrow\quad H_aF>e^{-\tau}.
\]

For cap two, the comparison \(t^p\le\Phi_N(t)\) on \(0\le t\le2\) is valid
for \(1<p\le p_N=1+\zeta(N-1)\). Tensor Hardy and Chebyshev give

\[
 L_N(q)=N(1+q^{-1})\log(1+q)+(1+q^{-1})\tau,\qquad
 q=(p-1)^{-1}\ge {1\over\zeta(N-1)}.
\]

Direct differentiation gives

\[
 L_N'(q)={N[q-\log(1+q)]-\tau\over q^2},
\]

so the boundary/interior split at
\(\tau_{\rm crit}=N[q_{\min}-\log(1+q_{\min})]\) is correct, with the
interior branch only for \(\tau>\tau_{\rm crit}\) and equality allowed at
the cap.

There is one precision repair in the current displayed boundary calculation:
the term is \(p_N\tau=O(N\tau)\), not \(O(\tau)\). Since
\(\tau\le\tau_{\rm crit}=O(N^{-1})\) on that branch, it is \(O(1)\).
Thus the correct boundary statement is \(L_N(q_{\min})=N+O(1)\), which is
fully sufficient for, and consistent with,

\[
 \log B_N(\tau)=N+\sqrt{2N\tau}+O(\tau+1)
\]

uniformly for \(\tau=o(N)\), including \(\tau\asymp N^{-1}\). The reported
outside repair to nearcorner76 resolves this local display issue.

The optional sharper inversion is correct. With \(x=\tau/N\),

\[
 q-\log(1+q)=x
 \quad\Longrightarrow\quad
 q=\sqrt{2x}+{2\over3}x+O(x^{3/2}),
\]

and hence, on the interior branch,

\[
 \log B_N(\tau)
 =N+\sqrt{2N\tau}+{2\over3}\tau+
 O\!\left({\tau^{3/2}\over\sqrt N}\right).
\]

Adding an \(O(1)\) term makes this compatible with the constrained boundary
branch. A local Wolfram series check gives
\(q-\log(1+q)=x-\frac19x^2+O(x^{5/2})\) after substituting the two displayed
inverse terms. One cosmetic TeX repair is also needed in nearcorner76
equation (2): insert the missing backslash before \quad.

## Orientations and conditional converses

The same-anchor all-orientations assertion is correct. A rectangle having
the fixed anchor as a vertex uses only the corresponding source orthant.
Translations and reflections reduce each orientation to one-sided tensor
Hardy, and these source orthants are disjoint up to null faces. Their
modulars sum to the one global modular, so no \(2^N\) factor belongs in the
same-anchor entropy.

The witness premise in entropy76 is appropriately stated as an imported
compact cap-two actual-function witness: r019 supplies the explicit
platform construction, while r069 fixed-cap repair preserves its calibration
and strictness. Given that premise and the additional full near-corner-cover
hypothesis,

\[
 |E_N|\le M_NB_N(\tau)J_F,\qquad
 {|E_N|\over J_F}\ge e^{q_0+\alpha\log N+o(\log N)}
\]

imply the necessary lower bound on \(\log M_N\). The parallel source-reuse
conclusion likewise follows only after assuming
\(\sum_jJ(F_{N,j})\le A_NJ(F_N)\). Neither is an unconditional statement
about arbitrary policies, stopping trees, or source allocations.

For \(\tau=o(N)\), the uniform Hardy expansion yields

\[
 \log M_N,\ \log A_N\ge(\mu-1)N-o(N).
\]

This is a barrier to a complete cheap finite-corner cover, not an
impossibility theorem for a global upper bound.

## Single-corner frontier

The frontier calculation is correct. On the interior branch,
\(\log B_N(\tau_*)=q_0+\alpha\log N\) is equivalent to

\[
 \tau_*=N[q_T-\log(1+q_T)],\qquad
 q_T=q_0/N+\alpha\log N/N-1.
\]

With \(q_0=N\mu-d_N\), \(d_N^2/N=\mu\log N+o(\log N)\), Taylor expansion at
\(\mu-1\) gives

\[
 \tau_*
 =N(\mu-1-\log\mu)-p_\mu d_N+
 \left(p_\mu\alpha+{1\over2\mu}\right)\log N+o(\log N).
\]

In particular the logarithmic coefficient
\(p_\mu\alpha+1/(2\mu)\) is correct. It describes only the
single-corner method/reference frontier, not an attained net or an actual
finite-dimensional endpoint constant.

## Final integration repairs

The final nearcorner76 display retains p_N tau explicitly and bounds it by O(1) only on the capped branch. The entropy note now also states that individual inequalities F_j<=F do not imply a bound on their summed modulars; its reuse inequality remains an independent, assumed budget condition.


---

# Review: translated copies separate anchor count from budget reuse

## Verdict

**PASS.**  For every fixed dimension and fixed enlargement \(\tau\),
translated compact cap-two packets give arbitrarily many disjoint witness
components that cannot be covered by a dimension-uniform number of
near-corners.  This is an anchor-count statement only.  A packetwise source
allocation has a reuse factor independent of the number of copies, so it
does not imply a comparable budget-reuse lower bound.

## 1. Copy construction

Let \(F_0\) be a compact \([0,2]\)-valued source with a compact witness
set \(E_0\).  Translate both by multiples of \(D e_1\):
\[
 F_j(y)=F_0(y-jDe_1),\qquad E_j=E_0+jDe_1,\qquad
 F_K=\sum_{j=1}^KF_j,\quad E_K=\bigcup_{j=1}^KE_j.
\tag{1}
\]
Choose \(D\) so the source supports and witness supports are disjoint.
Then \(0\le F_K\le2\), and translation invariance plus disjointness gives
\[
 J(F_K)=KJ(F_0),\qquad |E_K|=K|E_0|.
\tag{2}
\]
Thus the score \(|E_K|/J(F_K)\) is exactly that of one packet.

Let \(L_0\) bound the length of the first-coordinate projection of one
source packet.  The projection of \(\operatorname{supp}F_K\) is a union
of total length at most \(K L_0\).

## 2. One near-corner cannot serve two sufficiently separated copies

Let a signed fixed-corner rectangle from \(a\) to \(x\) have first-side
length \(d_1=|x_1-a_1|\).  Since \(F_K\le2\), first-coordinate Fubini
followed by the trivial bound one in every other coordinate gives
\[
 H_aF_K(x)
 \le {2\,|[a_1,x_1]\cap\operatorname{proj}_1(\operatorname{supp}F_K)|
          \over d_1}
 \le {2K L_0\over d_1}.
\tag{3}
\]
This holds for every signed orientation.  Hence
\[
 H_aF_K(x)>e^{-\tau}\quad\Longrightarrow\quad
 |x_1-a_1|<2e^\tau K L_0.                                \tag{4}
\]

Let \(D_E\) be the first-coordinate diameter of \(E_0\), and choose
\[
 D-D_E>4e^\tau K L_0.                                    \tag{5}
\]
If one anchor served points in two distinct witness copies, (4) and the
triangle inequality would force their first coordinates to differ by less
than \(4e^\tau K L_0\), contradicting (5).  Consequently every near-corner
cover of all of \(E_K\) uses at least \(K\) anchors.  Since \(K\) is
arbitrary for each fixed \(N\) and \(\tau\), there is no
dimension-uniform anchor-count bound for all such translated compact
inputs.

## 3. Why this is not a source-budget reuse lower bound

Suppose one packet is covered by \(M_0\) corners.  Assign its source piece
\(F_j\) to that packet's corners.  Reusing the same \(F_j\) for the
\(M_0\) corners costs at most
\[
 \sum_{j=1}^K\sum_{m=1}^{M_0}J(F_j)
 =M_0J(F_K).
\tag{6}
\]
Thus the allocation reuse factor is \(A=M_0\), independent of \(K\),
while the anchor count is \(KM_0\).  It would be incorrect to call the
reuse factor one unless \(M_0=1\).

All orientations at a fixed anchor remain free through the disjoint-orthant
accounting of `nearcorner76.md`; no extra \(2^N\) factor is used in (3)--(6).
The copies construction is not a general endpoint-upper counterexample: it
only shows that the anchor entropy and source-budget reuse ledgers in
`entropy76.md` are logically distinct and should not be conflated.


The final elementary packet specializes the projection proof to F0=3/2 on the unit cube, E0=[1/4,3/4]^N and D=8AK+2. It gives an exact K-anchor toy and reuse one; the general high-score statement retains reuse M0.
