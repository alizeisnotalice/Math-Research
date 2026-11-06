# r027: exact block likelihood log-concavity and safe Schur consequence

Date: 2026-11-06.  This is a strict auxiliary theorem for equal block
parameters.  It does not prove the proposed sharp curved-boundary asymptotic.

Let \(Z=\operatorname{Exp}(1)+c\operatorname{Bern}((1+e^c)^{-1})\), let
\(f_b\) be the density of a sum of \(b\) copies, and set \(R_b(t)=e^tf_b(t)\).
Positive constants and the scale \(c\) do not affect log-concavity.  At
\(c=1\), on the piece \(j\le x\le j+1\),
\[
 R_b(x)\doteq\sum_{k=0}^j\binom bk(x-k)^{b-1},\tag{1}
\]
with the last piece \(x\ge b\).

Wolfram 15 executed exact `Resolve[ForAll[x, dom => R R''-R'^2<=0],Reals]`
for every piece of \(b=4,5,6,7\), and returned `True` in all cases.  This is
a symbolic semialgebraic sign decision, not a grid test.  For \(b\ge3\),
the new truncated-power term has continuous first derivative at each knot.
Thus the piecewise log derivatives join continuously and the four checks give
global log-concavity of \(R_4,R_5,R_6,R_7\).

Since
\[
 R_{p+q}=R_p*R_q,\tag{2}
\]
Prékopa's convolution theorem and the four consecutive bases imply
\[
 \boxed{R_b\ \hbox{is log-concave on }[0,\infty)\quad\hbox{for every }b\ge4.}\tag{3}
\]
Use \(b=4,5,6,7\) as the bases modulo four and iterate (2).  Consequently
\(f_b=e^{-t}R_b\) is also log-concave.

## Safe implication for the symmetric KKT receiver

Take two identical blocks of size \(b\ge4\).  The capped KKT profile has
the form \(H_\lambda(t_1,t_2)=F(R_b(t_1)R_b(t_2))\) for a nondecreasing
scalar function \(F\).  At fixed total \(r=t_1+t_2\), concavity of
\(\log R_b\) makes \(R_b((r+z)/2)R_b((r-z)/2)\) even and nonincreasing in
\(|z|\).  Hence the same is true of \(H_\lambda\) along that line.

Conditioned on \(T_1+T_2=r\), the block difference has density proportional
to
\[
 f_b((r+z)/2)f_b((r-z)/2),
\]
which is symmetric decreasing because \(f_b\) is log-concave.  The
convolution of two nonnegative symmetric decreasing one-dimensional functions
is symmetric decreasing.  Conditioning on the total and then mixing proves
\[
 B(d,v)=\mathbb EH_\lambda(d+v+T_1,d-v+T_2)
 \quad\hbox{is even and nonincreasing in }|v|.\tag{4}
\]
Together with coordinate monotonicity in \(d\), every boundary selected by
\(B(d(v),v)=1\) is nondecreasing in \(|v|\).  This is useful global shape
control.  It does not supply the needed quantitative lower growth
\(d(v)-d(0)\sim v^2/(\mu n)\), so it does not by itself establish uniform
integrability or the \(e^{-1/2}\) score conjecture.
