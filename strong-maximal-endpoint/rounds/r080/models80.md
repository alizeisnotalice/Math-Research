# Sharp square-root scale and rational step-source persistence

This note complements `no_budget80.md`.  Its budget is the cost of a
fixed-corner source-column cover, not the weak-type constant of the strong
maximal operator.

## 1. Matching construction on a truncated halo patch

Retain the source \(F\), the parameter \(u\), and the patch coordinates of
`no_budget80.md`.  Thus
\[
 {2\over3}\le u\le{3\over4},\quad
 X={u+u^{-1}\over2},\quad
 Y=u+{1\over4u},\quad
 (x,y)=(X,Y-\delta).                                    \tag{1}
\]
For \(0<\varepsilon<1/96\), write
\[
 E_\varepsilon=\{(X,Y-\delta):2/3\le u\le3/4,
                       \ \varepsilon\le\delta<1/96\}. \tag{2}
\]
Choose a uniform rational \(\alpha\)-mesh of \([2/3,3/4]\) with
\[
 M=\left\lceil{1\over12\sqrt\varepsilon}\right\rceil,
 \qquad h={1\over12M}\le\sqrt\varepsilon,              \tag{3}
\]
and put one rational anchor at
\[
 a=\alpha,\qquad b={1\over2\alpha}                     \tag{4}
\]
for every mesh point.  Every column uses the unmodified source \(F\).

For the rectangle joining \((a,b)\) to \((X,Y-\delta)\), its average of
\(F\) exceeds one exactly when
\[
 P=(1-a^2)(1-b^2)-2(X-a)(Y-\delta-b)>0.                 \tag{5}
\]
Substitution of (1) and (4) gives the exact identity
\[
 P=2\delta(X-\alpha)-(u-\alpha)^2R(\alpha,u),
 \qquad
 R(\alpha,u)={4\alpha^2u^2-2\alpha u+1\over4\alpha^2u^2}.
\tag{6}
\]
For reference, after multiplying the \(\delta=0\) part by
\(4\alpha^2u^2\), (6) is the elementary factorization
\[
 4\alpha^2u^2P\big|_{\delta=0}
 =-(u-\alpha)^2(4\alpha^2u^2-2\alpha u+1).
\]
This verifies the identity without an optimization or computer algebra.

On the displayed parameter square,
\[
 R(\alpha,u)\le{73\over64},\qquad X-\alpha\ge{7\over24}. \tag{7}
\]
Indeed, with \(z=\alpha u\in[4/9,9/16]\),
\(R=1-(2z)^{-1}+(4z^2)^{-1}\) is decreasing, and its value at
\(z=4/9\) is \(73/64\).  A nearest mesh point has
\(|u-\alpha|\le h/2\).  Therefore, if \(\delta\ge\varepsilon\),
\[
 P\ge 2\varepsilon{7\over24}
       -{\varepsilon\over4}{73\over64}
 ={229\over768}\varepsilon>0.                           \tag{8}
\]
The \(M+1\) rational anchors consequently cover every point of
\(E_\varepsilon\) strictly.

The source cost is
\[
 J(F)={1+5\log2\over8};                                 \tag{9}
\]
for example, integrate \(2st\) on the region \(2st\le1\) and
\(2st(1+\log(2st))\) on \(2st>1\).  Thus the exhibited cover has cost
\((M+1)J(F)=O(\varepsilon^{-1/2})\).  This matches the lower square-root
scale from `no_budget80.md`: exact coverage of \(E_\varepsilon\) requires
\[
 B\ge {1\over73728\sqrt\varepsilon}.                   \tag{10}
\]
The discarded strip has measure \(|I|\varepsilon\).  Hence the construction
also gives a sequence of finite covers with budgets \(B_\varepsilon\) and
coverage deficits at most \(|I|\varepsilon=O(B_\varepsilon^{-2})\).
No uniform finite budget follows as \(\varepsilon\downarrow0\).

## 2. A rational compact step family with diverging column budget

For integers \(m\ge1\), let
\[
 F_m(s,t)=2\,{\lfloor ms\rfloor\over m}
                {\lfloor mt\rfloor\over m}\,
                {\bf1}_{[0,1]^2}(s,t).                  \tag{11}
\]
It is compactly supported, step-constant on a rational grid, and has
rational heights.  Moreover
\[
 0\le F_m\le F,\qquad \|F-F_m\|_\infty\le{4\over m},
 \qquad J(F_m)\le J(F).                                 \tag{12}
\]
The middle inequality follows from
\(st-\lfloor ms\rfloor\lfloor mt\rfloor/m^2\le2/m\) on the unit square.

On the patch in `no_budget80.md`, the exact formula
\(M_{\rm S}F=2u(x)u(Y(x)-\delta)\), together with
\(u(x)\ge2/3\) and \(|u'|\ge8/5\) on the relevant interval, gives
\[
 M_{\rm S}F(x,Y(x)-\delta)-1\ge{32\over15}\delta.      \tag{13}
\]
The maximal operator is \(L^\infty\)-Lipschitz.  Thus, for
\(\delta\ge4/m\),
\[
 M_{\rm S}F_m\ge M_{\rm S}F-{4\over m}>1.               \tag{14}
\]
When \(m>384\), the set \(E_{4/m}\) is nonempty and lies in the true halo
of \(F_m\).

Any source used in a fixed-corner column cover of that halo satisfies
\(G\le F_m\le F\).  Applying (10) with \(\varepsilon=4/m\) therefore
shows that its full true-halo column budget is at least
\[
 {\sqrt m\over147456}.                                  \tag{15}
\]
Since \(J(F_m)\le J(F)<\infty\), the normalized fixed-corner column budget
diverges at least on the order of \(\sqrt m\), even within rational compact
step inputs in fixed dimension two.

For each individual \(m\), `step_halo79.md` still gives a finite signed
corner representation of the halo of the actual step source \(F_m\), and
therefore a finite (but grid-dependent) fixed-corner cover using copies of
\(F_m\).  The divergence in (15) says precisely that these finite budgets
cannot be made uniform over this rational-step family.  It is not a claim
that the strong-maximal weak constant, or \(C_N\), diverges.

## 3. Actual modular bounds and bounded true halo volume

The logarithmic objective has not been replaced by an \(L^p\) proxy.
Since \(\Phi_2\) is Lipschitz with constant less than three on \([0,2]\),
(12) also gives
\[
 J(F)-12/m\le J(F_m)\le J(F).
\]
The lower bound is useful for the large meshes used in the certificate.
The upper bound for the actual halo follows from \(F_m\le2\,1_{[0,1]^2}\):
\[
 |\{M_{\rm S}F_m>1\}|
 \le 1+8\log2.
\]
Indeed, for \(A>1\) the constant-square halo has area
\(1+4A\log A\), by splitting its interior, four side strips and four
external corners. Thus the true halo/modular ratios of this family stay
bounded while its fixed-corner column-cover budget divided by the actual
modular tends to infinity. No divergence of an operator endpoint constant
is asserted.

## 4. Exact evaluation without forming the huge step array

Put \(f_m(s)=\lfloor ms\rfloor/m\) on the unit interval. For \(x>1\),
the step-source endpoint lemma reduces \(Mf_m(x)\) to anchors \(i/m\).
The exact tail integral is
\[
 \int_{i/m}^1 f_m(s)\,ds
 =\frac{m(m-1)-i(i-1)}{2m^2}.
\]
Thus its anchored average is this rational number divided by \(x-i/m\).
The rational interpolant of these grid values,
\[
 Q_m(a,x)=\frac{1-1/m-a^2+a/m}{2(x-a)},
\]
has derivative with the sign of
\(a^2-2xa+1+(x-1)/m\). This sign polynomial is strictly decreasing
on \(0\le a\le1\), positive at zero and negative at one. Hence the sampled
sequence is unimodal; exact checks at the two nodes bracketing its unique
continuous maximum, or equivalent neighbor checks, determine the entire
grid maximum. Finally
\[
 M_{\rm S}F_m(x,y)=2Mf_m(x)Mf_m(y).
\]
The large-m tests use these exact one-dimensional prefix sums, not an
instantiated \(m^2\)-element array.
