# Physical fixed-corner defect ledger

Let \(K=Q_{\Gamma,N}^{[A]}\) be the exact cap-\(A\) physical weak norm of
one signed fixed corner, as identified by r028 and `exact83.md`.  Consider a
finite or countable physical column family
\[
 E_i=\{H_{a_i}G_i>1\},\qquad 0\le G_i\le F\le A,
 \qquad w_i\ge0.
\]
Let \(\rho\ge0\) be an integrable target density on physical Lebesgue
space, and set
\[
 \Gamma=\sum_iw_i{\bf1}_{E_i},\qquad
 B=\sum_iw_iJ_N(G_i),\qquad V=\int\rho,
 \qquad D=\int(\rho-\Gamma)_+ .                         \tag{1}
\]
All statements below assume these quantities are finite.  In particular,
the fixed-corner inequality \(|E_i|\le KJ_N(G_i)\) and Tonelli imply
\(\int\Gamma<\infty\) whenever \(B<\infty\).

## Exact identity

The following is an equality, not an estimate:
\[
 \boxed{
 KB+D-V
 =\sum_iw_i\bigl(KJ_N(G_i)-|E_i|\bigr)
   +\int(\Gamma-\rho)_+ .}                              \tag{2}
\]
Indeed, the pointwise elementary identity
\[
 (\rho-\Gamma)_+-\rho=(\Gamma-\rho)_+-\Gamma
\]
and Tonelli give
\[
 KB+D-V
 =KB-\int\Gamma+\int(\Gamma-\rho)_+,
\]
which is exactly (2).  Every term on the right is nonnegative by the exact
fixed-corner norm bound.

## Consequences of a near-multiplicative budget

Suppose the proposed column budget is nearly multiplicative in the sense
that
\[
 T:=KB+D\le(1+\delta)V.                                 \tag{3}
\]
Then (2) gives the two necessary bounds
\[
 \int(\Gamma-\rho)_+\le\delta V,                       \tag{4}
\]
and, for every \(\varepsilon>0\),
\[
 \sum_{i:\,|E_i|\le(1-\varepsilon)KJ_N(G_i)}
 w_iJ_N(G_i)
 \le {\delta V\over\varepsilon K}.                    \tag{5}
\]
The second follows because every listed low-efficiency column contributes at
least \(\varepsilon K w_iJ_N(G_i)\) to the first term on the right of (2).

Thus any genuinely near-optimal physical column construction must place
almost all of its source budget on individually near-extremal fixed-corner
columns and must have negligible overcoverage.  This is only a necessary
condition for a near-multiplicative budget.  The source-sharing LP of r082
can still be degenerate, and an \(o(\log N)\) target error does not by itself
imply \(\delta\to0\).  Therefore this ledger does not close a general upper
bound.

There is an additional logical requirement before applying (5) to a real
covered main body.  Condition (3) alone permits \(B=0,D=V\), hence
\(T=V\).  One must separately establish an actual coverage-deficit estimate
\[
 D\le\varepsilon_RV.                                    \tag{6}
\]
Then the integrated covering inequality
\(\rho\le\Gamma+(\rho-\Gamma)_+\), together with
\(\int\Gamma\le KB\), yields
\[
 KB\ge V-D\ge(1-\varepsilon_R)V,
 \qquad KB\le(1+\delta)V.                              \tag{7}
\]
This shows that a nontrivial amount of source budget is present, but it
still does not remove the cross-corner/source-sharing degeneracy of the LP.
The defect ledger must not be used as a substitute for (6).

The bookkeeping is explicitly on physical Lebesgue measure.  It cannot be
replaced by a finite list of atomic receiver samples without separately
justifying the measure identity and its limiting passage.  Its purpose is to
focus the remaining work on cross-corner reuse of sources; the one-corner
constant itself has already been reduced exactly.
