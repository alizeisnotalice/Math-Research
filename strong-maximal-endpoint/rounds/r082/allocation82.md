# Fixed-corner charge forces a revision of the r081 soft budget target

## Verdict

The r081 target for \(\mathcal W_{N,N^{-2}}(3)\) is impossible once the
optimized capped fixed-corner charge \(e^{q_{N,3}}\) is used.  The soft budget
\(\mathcal W\) has **not** changed.  What changes is a valid implication
from \(\mathcal W\) to the soft norm: the smaller column coefficient makes
the necessary lower bound on this unchanged \(\mathcal W\) stronger.

There is also a rigorous finite-dimensional instance at \(N=64\):
\[
 \mathcal W_{64,1/4096}(3)>
 \exp(118915/8000)>2.85\cdot10^6. \tag{1}
\]
The finite statement needs only the certified cap-two witness from r054 and
the directed cap-three root bracket from r082; it does not appeal to an
unproved asymptotic lower bound for a cap-restricted norm.

## 1. Correct cap and norm chain

Let \(C_{N,A}\) be the hard endpoint norm for \(0\leq F\leq A\), and let
\(S^+_{N,A}(\eta)\) be its upward flat-smoothed counterpart.  At
\(\eta=N^{-2}\), exact Orlicz dilation gives
\[
 C_{N,2}\leq D_N(1+\eta)S^+_{N,\,2(1+\eta)}(\eta)
 \leq D_N(1+\eta)S^+_{N,3}(\eta), \tag{2}
\]
because \(2(1+N^{-2})\leq3\).  This is a cap enlargement, not an
identification of cap two and cap three.  Also
\[
 \log D_N(1+N^{-2})=O(N^{-1}). \tag{3}
\]

For a real fixed-corner column cover in the definition of
\(\mathcal W_{N,\eta}(3)\), the fixed-corner theorem applies to every
column \(G\leq F\leq3\).  Integrating its weighted cover and then taking
the infimum/supremum gives
\[
 S^+_{N,3}(\eta)\leq e^{q_{N,3}}\mathcal W_{N,\eta}(3). \tag{4}
\]
The fixed-corner proof already sums the signed orientations through
disjoint source orthants, so (4) has no hidden \(2^N\) factor.

The r073 cap-two amplitude reduction says, in its stated large-\(N\)
regime,
\[
 C_N\leq(1+\varepsilon_N)C_{N,2},\qquad \varepsilon_N=o(1). \tag{5}
\]
Thus any established global lower witness \(C_N\geq Q_N\) yields
\[
 \mathcal W_{N,N^{-2}}(3)\geq
 {Q_N\over(1+\varepsilon_N)D_N(1+N^{-2})e^{q_{N,3}}}. \tag{6}
\]
If \(\mathcal W=+\infty\), this is still true as an extended-real
inequality.  Hence no finiteness of the missing general column budget is
being assumed.

## 2. Consequence for the asymptotic r081 target

The established asymptotic inputs are
\[
 \log Q_N=\mu N-d_\mu+{9-\mu\over24}\log N+o(\log N) \tag{7}
\]
and, at every fixed cap including \(A=3\),
\[
 q_{N,3}=N-d_1+{1\over3}\log N+o(\log N), \tag{8}
\]
where \(d_1=\sqrt N\,z_1\) and
\(\overline\Phi(z_1)=z_1/\sqrt N\).  Substitution in (6) gives
\[
 \boxed{\quad
 \log\mathcal W_{N,N^{-2}}(3)\geq
 (\mu-1)N-d_\mu+d_1+
 \left({9-\mu\over24}-{1\over3}\right)\log N
 +o(\log N).\quad} \tag{9}
\]
Mills' inequalities imply \(z_1^2=\log N+O(\log\log N)\); hence
\[
 d_1-\tfrac13\log N\gg\log N. \tag{10}
\]
The prior r081 upper target omitted precisely
\(d_1-\frac13\log N\), so it contradicts (9), rather than merely having a
different constant hidden in its \(o(\log N)\).

A scale-consistent *sufficient target* for the earlier endpoint comparison
is instead the matching upper bound
\[
 \log\mathcal W_{N,N^{-2}}(3)\leq
 (\mu-1)N-d_\mu+d_1+
 \left({9-\mu\over24}-{1\over3}\right)\log N
 +o(\log N). \tag{11}
\]
Adding (8) to (11) reproduces the desired old endpoint scale.  Equation
(11) is not proved here, and the argument does not show that
\(\mathcal W\) is uniformly finite.

## 3. Independent finite-\(N\) obstruction

The r054 directed-Decimal certificate supplies a **finite, compact,
cap-two** input \(F\) in dimension \(64\) with
\[
 {|\{M_{\rm S}F>1\}|\over J_{64}(F)}>e^{71.79}. \tag{12}
\]
This is the explicit finite-\(L=10^8\) scalar-staircase transfer recorded
in [r054/scalar54.json](../r054/scalar54.json), with the explicit proof in [scalar_proof54.md](../r054/scalar_proof54.md), not an inference that a general
cap-two norm has an unproved lower bound.

Set \(\eta=1/4096\) and \(F'=(1+\eta)F\).  Then \(F'\leq3\), and on the
hard halo of \(F\),
\[
 M_{\rm S}F'>1+\eta,\qquad
 \chi((M_{\rm S}F'-1)/\eta)=1.
\]
The exact dilation bound gives
\[
 { \int\chi((M_{\rm S}F'-1)/\eta)\over J_{64}(F')}
 >{e^{71.79}\over D_{64}(1+\eta)}
 \geq e^{71.79-1/64}. \tag{13}
\]
The 32-cell/positive-midpoint directed r082 certificate at
\((N,A)=(64,3)\) proves that the calibrated mean at \(q=56.91\) is below
one, so the unique root satisfies
\[
 q_{64,3}<56.91. \tag{14}
\]
Apply (4) to (13), then use (14):
\[
 \log\mathcal W_{64,1/4096}(3)
 >71.79-\frac1{64}-56.91
 ={118915\over8000}. \tag{15}
\]
The final comparison \(e^{118915/8000}>2.85\cdot10^6\) is reducible to a
positive Taylor lower bound for the exponential with rational input, so it
can be certified without floating-point inference.

## Scope

The smaller \(e^{q_{N,3}}\) is a better fixed-corner **upper charge**.
It does not cause the true optimum \(\mathcal W\) to change or become
larger.  It strengthens the lower requirement because (4) is stronger.
Nothing above proves a source-sharing/packing theorem, the finiteness of
\(\mathcal W\), or a new general strong-maximal upper bound.

## 第84轮后续更正（2026-11-11）

本轮提出但未证明的充分目标(11)现已被 [第84轮严格反例](../r084/obstruction84.md) 否定：保持同一预算定义，已有 cap2 输入给 $\log\mathcal W_{N,N^{-2}}(2)\ge3N\log N/2-O(N)$。本文件的必要下界(9)和有限维证书仍有效；(11)不再列为可攻关目标。
