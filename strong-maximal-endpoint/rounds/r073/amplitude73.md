# Amplitude localization audit for the wide-window module

## Verdict

There is a valid **global physical reduction** from arbitrary amplitudes
to the cap-two class.  It uses a much smaller splitting threshold than the
older \(2N^{-1/2}\) display.  This reduction is distinct from, and does
not by itself solve, the fixed-reference local mean and anchor steps in
r072.

## What is already proved

Let \(\rho_N=\sqrt{\log N/N}\), fix \(M>1\), and set
\[
 g=f{\bf1}_{\{f>e^{A\rho_N}\}}
\]
and \(s=N^{-M}\).  The archived proof with threshold \(2s\), rather than
only its displayed \(2N^{-1/2}\) specialization, gives
\[
 |\{M_Sg>2s\}|
 \le B^N e^{-(A-\sqrt{2\mu M})\sqrt{N\log N}+O_{A,M}(\log N)}
 \int\Phi_N(g).
\tag{1}
\]
For the middle range take \(p=b\sqrt{N/\log N}\).  The logarithmic
coefficient before \(-A\sqrt{N\log N}\) is
\[
 {\mu\over2b}+Mb,
\]
whose minimum is \(\sqrt{2\mu M}\), at
\(b=\sqrt{\mu/(2M)}\).  In the extreme range the replacement of
\(N^{-1/2}\) by \(N^{-M}\) adds only \(M\log N\) to the old factorial
kernel estimate, hence remains superexponentially small.  Thus (1) is a
valid extension of the existing localization argument.

The retained core
\[
\ell=f{\bf1}_{\{f\le e^{A\rho_N}\}}
\]
satisfies \(\ell/(1-2s)<2\) for all large \(N\), and
\[
 \int\Phi_N\!\left({\ell\over1-2s}\right)
 \le(1-2s)^{-N}\int\Phi_N(\ell)
 =(1+O(N^{1-M}))\int\Phi_N(\ell).                          \tag{2}
\]

## The actual r072 requirement

Let \(T_{N,2}\) denote the genuine physical weak ratio supremum restricted
to inputs bounded by two.  Subadditivity at levels \(1-2s,2s\), (1), and
(2) yield
\[
 C_N\le\max\left\{(1-2s)^{-N}T_{N,2},\,H_{N,A,M}\right\},   \tag{3}
\]
where \(H_{N,A,M}\) is the coefficient in (1).  The cap-two subclass
gives \(C_N\ge T_{N,2}\).  Since the cap-two scalar witness has
\[
 T_{N,2}\ge Q_{N,2}
 =\exp\{N\mu-\sqrt\mu\sqrt{N\log N}+o(\sqrt{N\log N})\},
\]
choose
\[
 A>\sqrt\mu+\sqrt{2\mu M}.
\]
Then \(H_{N,A,M}/T_{N,2}\to0\), while
\((1-2s)^{-N}=1+O(N^{1-M})\).  Therefore
\[
 \boxed{C_N/T_{N,2}\longrightarrow1.}                     \tag{4}
\]
This is a general physical cap reduction, not merely an r061 profile
statement.

### Explicit high-amplitude splitting bound

For completeness, fix a physical cutoff a0>1 and 0<s<1/2. Split g=f 1_{f>a0} into g_mid with log f<=sqrt N and g_ext with log f>sqrt N. The archived sharp strong-Lp estimate, with b_p=p log c_p determined by p(1−e^{−b_p/p})=1+e^{−b_p}, and the convexity in u=log f give

\[
 C_{\rm mid}=e^{Nb_p}s^{-p}\max_{u\in\{\log a_0,\sqrt N\}}
 \exp\{(p-1)u-(N-1)\log(1+u)\}.
\]

The r001 iterated truncation bound, using r=1/2 and monotonicity in u>=sqrt N, gives

\[
 C_{\rm ext}=\frac{s^{-1}4^N}{\Gamma(N)}
 \left(\frac{\sqrt N+N\log2+\log(1/s)}{1+\sqrt N}\right)^{N-1}.
\]

Consequently |{M_S g>2s}|<=H_N integral Phi_N(g), with H_N=max(C_mid,C_ext), because the two source parts are disjoint. Using C_mid+C_ext instead is a harmless conservative coefficient in the numerical scripts. No subadditivity of Phi_N on overlapping inputs is used.

### Exponentially accurate cap-two reduction

The preceding polynomially small threshold can be strengthened to an
exponentially small one.  Take
\[
 a_0={3\over2},\qquad s=e^{-N/100},
\]
split at \(a_0\), and use the fixed strong-type exponent \(p=4\) on the
middle range.  Its lower-end logarithmic rate is
\[
 r_0=b_4+{4\over100}-\log(1+\log(3/2)).                    \tag{5}
\]
It is strictly below \(\mu\).  A convenient rational verification is
\[
 b_4<{3\over2},\qquad \log(1+\log(3/2))>{1\over3},
\]
so
\[
 r_0<{181\over150}<{5\over4}<\mu.                          \tag{6}
\]
The needed elementary exponential bounds are
\[
 e^{3/8}>{16\over11},\qquad e^{2/5}<{3\over2},\qquad
 e^{1/3}<{7\over5},\qquad e^{5/4}<4.
\]
For the first implication put z=e^{−3/8}<11/16. Then 3−4z−z^4>3−4(11/16)−(11/16)^4=1743/65536>0, so the increasing b4 root equation is positive at b=3/2. Together these facts certify the displayed \(b_4\), logarithmic, and
\(\mu\) comparisons by the defining equation for \(b_4\) and
\(\mu=1+e^{-\mu}\).

At the upper endpoint \(u=\sqrt N\), the middle-range exponent is
\(-\tfrac12N\log N+O(N)\).  The iterated weak-kernel extreme range has
the same bound: replacing the old polynomial threshold by
\(e^{-N/100}\) only adds \(N/100\) to its \(O(N)\) term.  Thus the
high-amplitude coefficient obeys
\[
 H_N\le e^{r_0N+o(N)}.                                     \tag{7}
\]
Since \(T_{N,2}\ge e^{\mu N-o(N)}\), \(H_N/T_{N,2}\) is exponentially
small.  The low input satisfies
\[
 {\ell\over1-2e^{-N/100}}<2,\qquad
 (1-2e^{-N/100})^{-N}=1+O(Ne^{-N/100}).
\]
The disjoint-cost maximum argument in (3) therefore gives the sharper
physical equivalence
\[
 \boxed{\displaystyle
 1\le {C_N\over T_{N,2}}
 \le1+O(Ne^{-N/100}).}                                    \tag{8}
\]

More generally, every fixed cap \(A>1\) admits an exponential reduction.
Choose \(1<a_0<A\), then choose \(p\) so large that
\[
 b_p-\mu<{\log(1+\log a_0)\over4},
\]
and put
\[
 \kappa_A={\log(1+\log a_0)\over4p},\qquad s=e^{-\kappa_A N}.
\]
The middle lower-end rate is then strictly below \(\mu\), while the
rescaled low cap is below \(A\) for large \(N\).  The same proof gives
\[
 {C_N\over T_{N,A}}=1+O_A(Ne^{-\kappa_A N}).                \tag{9}
\]

At a retained reference receiver, r072 additionally needs
\[
 \mathbb E_c h(S)\ge1-\eta_{\rm ref,N},\qquad
 \eta_{\rm ref,N}=o(\log N/N).                              \tag{10}
\]
After the radial tail is removed, rescaling costs
\[
 (1-\eta_{\rm ref,N}-O(N^{-M}))^{-N}=N^{o(1)}.              \tag{11}
\]
An \(O(\log N/N)\) deficit instead causes an \(N^{O(1)}\) multiplier, too
large for the arbitrarily small logarithmic loss in the r069 capture
bookkeeping.

The cap reduction (4) does not give a statement of this receiverwise form.
A global maximal-tail bound cannot be converted to this
fixed-reference receiver mean statement without the missing
common-source/reference/anchor mechanism.

## Why clipping is not a valid shortcut

For a naive local clipping \(F^{(A_0)}=F\wedge A_0\), the removed mean at one receiver is
\[
 \Delta_{A_0}=\int p_c(s)
 [K_c^s(F-F^{(A_0)})_{u+s}](X,\sigma)\,ds.
\]
The global modular controls an integral of \(\Phi_N(F)\), not this local
cone quantity.  It may concentrate along one receiving cone or angular
fibre.  Hence naive local clipping does not supply
\(\Delta_{A_0}=o(\log N/N)\).  This does not contradict (4), which charges
the entire high-amplitude maximal halo separately before passing to the
cap-two supremum.

A growing cap does not fix this automatically.  The r072 tail and
rescaling constants then depend on the cap.  To retain (2), one needs an
explicit increasing-cap estimate, for example a condition controlling
\(N\eta_N\) at \(o(\log N)\) scale.  No such quantified lemma appears in
r001 or r061.

## Relation to r061

Section 7 of r061 proves a different result.  Given an *actual local
modular budget* for fixed balanced blocks and a fixed cap, a successful
terminal is \(O(\rho_N)\)-close to one in \(L^1\), and is
\(O(\rho_N^{1/2})\)-close in \(L^2\).  Its unbounded scalar feedback
extension uses a separate Fenchel tail argument.  Neither statement turns
an arbitrary physical common source into a fixed-cap source with the
receiverwise loss (1).

The long angular-band example in r061 likewise shows why a global modular
budget cannot simply be substituted for an actual local cone budget.

## Remaining issues are not amplitude issues

Equations (8)--(9) complete the general amplitude reduction.  The
remaining r072 tasks are genuinely different: reference-mean comparison
for variable-action successful receivers, common-source/anchor
compatibility, and packing.  They should not be recorded as an amplitude
localization remainder.
