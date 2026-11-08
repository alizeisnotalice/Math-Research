# Smooth endpoint norms: exact scaling comparison and its limits

Let
\[
 \Phi_N(t)=t(1+\log^+t)^{N-1},\qquad
 J_N(F)=\int\Phi_N(F),\qquad
 C_N=\sup_F { |\{M_{\rm S}F>1\}|\over J_N(F)}.         \tag{1}
\]
The suprema below range over nonzero nonnegative inputs of finite modular.

## 1. A flat smooth cutoff

Define
\[
 \chi(z)=
 \begin{cases}
 0,&z\le0,\\
 \displaystyle {e^{-1/z}\over e^{-1/z}+e^{-1/(1-z)}},&0<z<1,\\
 1,&z\ge1.
 \end{cases}                                             \tag{2}
\]
For \(0<z<1\),
\[
 \log {\chi(z)\over1-\chi(z)}=-{1\over z}+{1\over1-z},
\]
whose derivative is \(z^{-2}+(1-z)^{-2}>0\).  Thus \(\chi\) is
nondecreasing (strictly increasing on \((0,1)\)).

For completeness, \(r(z)=e^{-1/z}\) for \(z>0\), extended by zero for
\(z\le0\), is \(C^\infty\) and flat at zero: every derivative is a
polynomial in \(1/z\) times \(e^{-1/z}\), hence tends to zero as
\(z\downarrow0\).  Near zero, (2) is \(r(z)\) divided by a smooth function
bounded away from zero.  It is therefore flat at zero; the symmetry
\(\chi(1-z)=1-\chi(z)\) proves the same assertion at one.  Thus (2) is a
nondecreasing \(C^\infty\) cutoff on the real line.

For \(\eta>0\), set
\[
 S_N^+(\eta)=\sup_F
 {\displaystyle\int\chi\!\left({M_{\rm S}F-1\over\eta}\right)
  \over J_N(F)}.                                        \tag{3}
\]
Pointwise,
\[
 {\bf1}_{\{M_{\rm S}F>1+\eta\}}
 \le\chi\!\left({M_{\rm S}F-1\over\eta}\right)
 \le {\bf1}_{\{M_{\rm S}F>1\}}.                       \tag{4}
\]

## 2. Exact Orlicz dilation and the upper smoothing sandwich

For \(a\ge1\), the exact dilation constant of \(\Phi_N\) is
\[
 D_N(a)=a(1+\log a)^{N-1}
   =\sup_{t>0}{\Phi_N(at)\over\Phi_N(t)}.               \tag{5}
\]
Indeed, for \(t\ge1\), divide the two factors in (5) and use
\((1+\log(at))/(1+\log t)\le1+\log a\).  For
\(a^{-1}<t<1\), the denominator is \(t\) and the same bound is immediate;
for \(t\le a^{-1}\), the ratio is \(a\).  Equality in (5) occurs at
\(t=1\).

The right inequality in (4) gives \(S_N^+(\eta)\le C_N\).  Conversely,
apply (4) to \(F'=(1+\eta)F\).  If \(M_{\rm S}F>1\), then
\(M_{\rm S}F'>1+\eta\), so the smooth integrand for \(F'\) equals one.
By (5),
\[
 |\{M_{\rm S}F>1\}|
 \le S_N^+(\eta)J_N((1+\eta)F)
 \le D_N(1+\eta)S_N^+(\eta)J_N(F).
\]
Taking the supremum proves the exact comparison
\[
 \boxed{\quad S_N^+(\eta)\le C_N
 \le D_N(1+\eta)S_N^+(\eta).\quad}                     \tag{6}
\]

For a cap \(A>0\), attach a subscript \(A\) to restrict the inputs to
\(0\le F\le A\).  The cap-sensitive statement is
\[
 \boxed{\quad S_{N,A}^+(\eta)\le C_{N,A}
 \le D_N(1+\eta)S_{N,(1+\eta)A}^+(\eta).\quad}         \tag{7}
\]
The enlarged cap on the right is necessary because the proof scales the
source up.  In particular, (7) does not identify cap two with cap
\(2(1+\eta)\).

## 3. Downward threshold smoothing

For \(0<\eta<1\), define
\[
 S_N^-(\eta)=\sup_F
 {\displaystyle\int\chi\!\left({M_{\rm S}F-(1-\eta)\over\eta}\right)
  \over J_N(F)}.                                        \tag{8}
\]
This cutoff is one whenever \(M_{\rm S}F\ge1\), and vanishes whenever
\(M_{\rm S}F\le1-\eta\).  Hence
\[
 C_N\le S_N^-(\eta)
 \le D_N((1-\eta)^{-1})C_N.                             \tag{9}
\]
For the second inequality, bound the integrand by
\({\bf1}_{\{M_{\rm S}F>1-\eta\}}\), then apply the definition of \(C_N\)
to \(F/(1-\eta)\).  At cap \(A\), this reads
\[
 \boxed{\quad C_{N,A}\le S_{N,A}^-(\eta)
 \le D_N((1-\eta)^{-1})C_{N,A/(1-\eta)}.\quad}         \tag{10}
\]
Thus the downward version has the opposite cap enlargement, again essential.

## 4. The \(\eta=N^{-2}\) scale and what it does not solve

For \(N\ge2\), take \(\eta=N^{-2}\).  Since
\(\log(1+\eta)\le\eta\),
\[
 \log D_N(1+\eta)
 =\log(1+\eta)+(N-1)\log(1+\log(1+\eta))
 \le N\eta={1\over N}.                                 \tag{11}
\]
For the lower threshold, \(-\log(1-\eta)\le\eta/(1-\eta)\le4\eta/3\),
so \(\log D_N((1-\eta)^{-1})\le4/(3N)\).  A local Python calculation
checked the corresponding \(\log D\) values at \(N=2,10,100,1000\);
these numerical checks are only sanity checks, while (11) is the proof.

Consequently smoothing at this width has multiplicative loss
\(e^{O(1/N)}\).  It does not by itself provide a new sharp upper bound:
one would still need a noncircular estimate of \(S_N^+(N^{-2})\) (or of the
downward functional) that is sharper than applying the definition of
\(C_N\) back to it.  In the bounded-input setting one additionally needs
such an estimate at the enlarged cap \(2(1+N^{-2})\), or a separately
quantified return from that cap to two.  The r073 cap-two amplitude
reduction must not be silently treated as exact cap preservation here.

Thus (6)--(10) are useful norm equivalences and regularizations, not a
closure of the general endpoint problem.
