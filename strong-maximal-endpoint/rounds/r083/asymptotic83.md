# Actual corner norm versus radial Fenchel charge

## Verdict

For every fixed cap \(A>1\),
\[
 \log Q_{\Gamma,N}^{[A]}
 =N-\sqrt N z+\frac {z^2}{3}-\frac16+o(1),
 \qquad \overline\Phi(z)=z/\sqrt N. \tag{1}
\]
The r082 proof-class charge is
\[
 q_{N,A}=N-\sqrt N z+\frac {z^2}{3}+\frac43+o(1).
 \tag{2}
\]
Therefore
\[
 q_{N,A}-\log Q_{\Gamma,N}^{[A]}\longrightarrow\frac32. \tag{3}
\]
The rough finite-dimensional impression of a gap near \(2\) is not the
limit.  This is a comparison between the r082 radial
Fenchel--Fubini--Markov class and the actual *fixed-corner* norm.  It says
nothing about moving corners or a general strong-maximal upper bound.

## Exact scalar KKT problem

The r028 equality at ratio \(r=0\) identifies the actual fixed-corner
cap-\(A\) norm with
\[
 Q_{\Gamma,N}^{[A]}=(\min J_N(h))^{-1},\quad
 J_N(h)=\int_0^\infty e^{-s}\Phi_N(h(s))\,ds,\quad
 \mathbb E_\Gamma h(S)=1. \tag{4}
\]
Let \(b\) be the zero-to-one activation point.  The KKT equation is
\(\Phi_N'(h(s))=(s/b)^{N-1}\) in the middle, with the expected zero,
one, and cap-\(A\) pieces.  Thus, with
\(Q_\Gamma(x)=\mathbb P(\Gamma(N,1)\geq x)\),
\[
 Q_\Gamma(b)+\int_1^A
 Q_\Gamma\left(b[\Phi_N'(t)]^{1/(N-1)}\right)\,dt=1, \tag{5}
\]
and
\[
 e^bJ_N(h_b)=1+\int_1^A\Phi_N'(t)
 e^{-b([\Phi_N'(t)]^{1/(N-1)}-1)}\,dt. \tag{6}
\]

## Expansion of the mean root

Set
\[
 x=z^2,\quad r=\sqrt N\phi(z),\quad q_0=N-\sqrt Nz,\quad b=q_0+u.
\]
Uniformly for \(u=O(x)\), the Gamma expansion used in r082 gives
\[
 Q_\Gamma(q_0+u)=1-\frac z{\sqrt N}
 -\frac rN\left(u+\frac{1-x}{3}\right)+o(x/N). \tag{7}
\]
For the integral in (5), put \(t=1+v/\sqrt N\), \(L=\log N\).  In its
contributing range,
\[
 \begin{aligned}
 \log\Phi_N'(t)&=L+\sqrt Nv-v^2+O((v+v^3)/\sqrt N),\\
 [\Phi_N'(t)]^{1/(N-1)}
 &=1+\frac v{\sqrt N}+\frac{L-v^2/2}{N}
   +O((Lv+v^3)/N^{3/2}),\\
 \frac{b[\Phi_N'(t)]^{1/(N-1)}-N}{\sqrt N}
 &=v-z+\frac{u+L-v^2/2-zv}{\sqrt N}+o(N^{-1/2}).
 \end{aligned} \tag{8}
\]
The survival expansion
\[
 Q_\Gamma(N+\sqrt Ny)=\overline\Phi(y)
 +\frac{y^2-1}{3\sqrt N}\phi(y)+\text{ integrated }o(x/N) \tag{9}
\]
and \(y=v-z\) now give
\[
 \int_1^A Q_\Gamma\left(b[\Phi_N'(t)]^{1/(N-1)}\right)dt
 =\frac z{\sqrt N}+\frac{-u+\frac32x-L+o(x)}N. \tag{10}
\]
This uses the three elementary truncated Gaussian moments
\[
 \int_{-z}^\infty\phi=\Phi(z),\quad
 \int_{-z}^\infty y\phi=\phi(z),\quad
 \int_{-z}^\infty y^2\phi=\Phi(z)-z\phi(z).
\]
Gamma Chernoff bounds and monotonicity of \(\log\Phi_N'\) make the omitted
tails \(o(x/N)\).

For the word “integrated” in (9), an \(L^1\) Edgeworth remainder by itself
should not be quoted without its tail bookkeeping.  The needed statement is
obtained by cutting at \(|y|\leq N^{1/12}\), where the direct log-density
Taylor remainder has weighted \(L^1\) norm \(O(N^{-1})\), and using the
Gamma Chernoff bound outside.  Fubini then converts the survival remainder
over the \(v\)-integral into that weighted density remainder, followed by
the prefactor \(N^{-1/2}\).  This is \(o(x/N)\).  Thus the claimed precision
does not follow from an unweighted Berry--Esseen estimate.

Combining (7) and (10),
\[
 m_N(q_0+u)-1=
 \frac{-(r+1)u+r(x-1)/3+\frac32x-L+o(x)}N. \tag{11}
\]
The differentiated expansion has slope
\(-(r+1)N^{-1}(1+o(1))\), hence
\[
 u=\frac{r(x-1)/3+\frac32x-L}{r+1}+o(1). \tag{12}
\]
Mills' relation from \(\overline\Phi(z)=z/\sqrt N\) is
\[
 r=x+1+O(x^{-1}),\qquad L=x+4\log z+\log(2\pi)+O(x^{-1}).
\]
Equation (12) consequently becomes
\[
 b=q_0+\frac x3-\frac16+o(1). \tag{13}
\]
The \(o(x/N)\) precision in (11) is essential: it becomes \(o(1)\) after
division by the mean slope.

## Modular correction with an explicit bound

Let \(h_b'(s)\) denote the ordinary derivative in the middle branch and
zero on the constant branches; the activation jump is accounted separately.
Integration by parts in the physical exponential measure gives
\[
 e^bJ(h_b)-1=\lambda_b e^b\,\mathbb E_\Gamma h_b'(S)
 ={\mathbb E_\Gamma h_b'(S)\over f_N(b)}.
\tag{14}
\]
In the middle branch, implicit differentiation of the KKT equation gives
\[
 h_b'(s)={h_b(s)(1+\log h_b(s))(N+\log h_b(s))
 \over s(N-1+\log h_b(s))}
 \le {2A(1+\log A)\over s}\le {C_A\over b}.
\tag{15}
\]
The density expansion at the calibrated root gives
\(f_N(b)\sim r/N\), while \(b\sim N\). Thus
\[
 0\le e^bJ(h_b)-1=O(r^{-1})=O(z^{-2})=o(1).
\tag{16}
\]
In fact, on a Gamma event of probability tending to one, the profile is in
its middle branch, \(h_b(S)\to1\) and \(S/N\to1\), so
\(Nh_b'(S)\to1\) in probability. Bound (15) permits dominated
convergence, hence \(e^bJ(h_b)-1\sim1/r\). No unproved saddle-tail
estimate is needed for (16).

Consequently
\[
 \log Q_{\Gamma,N}^{[A]}=b-\log(e^bJ(h_b))
 =N-\sqrt N z+\frac{z^2}{3}-\frac16+o(1),
\]
which proves (1)--(3). This is the actual fixed-physical-corner norm, by
the already established r028 reduction; it does not compare variable
physical corners or close the general strong-maximal upper bound.
