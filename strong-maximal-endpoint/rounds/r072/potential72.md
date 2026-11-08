# The Cartesian kernel is not a subordinator potential

Let
\[
 k_N(z)={\Gamma(N)\over(\sum_{i=1}^Nz_i)^{N-1}}
 {\bf1}_{\{z_i\ge0\}},
\]
with the inessential normalization chosen so that its multivariate Laplace
transform is
\[
 F_N(a)=\mathbb E_\Theta {1\over a\cdot\Theta},\qquad
 \Theta\sim {\rm Dirichlet}(1,\ldots,1).
\]
Put \(\phi_N=F_N^{-1}\).  If \(k_N\) were the potential density, up to a
positive scalar, of a killed or unkilled positive convolution semigroup
of subordinators, then \(\phi_N\) would be its Laplace exponent.  In
particular it would be a multivariate Bernstein function: every
\(\partial_i\phi_N\) would be completely monotone.  A necessary condition
is therefore
\[
 \partial_i\partial_j\phi_N(a)\le0\qquad(i\ne j).           \tag{1}
\]

This condition fails for every \(N\ge2\).  Write \(a_i=1+x_i\) and
\[
 L=\sum_i\Theta_i x_i.
\]
Near \(x=0\),
\[
 F_N(1+x)=1-\mathbb EL+\mathbb EL^2+O(|x|^3),
\]
and hence
\[
 \phi_N(1+x)
 =1+\mathbb EL+(\mathbb EL)^2-\mathbb EL^2+O(|x|^3).        \tag{2}
\]
The Dirichlet moments are
\[
 \mathbb E\Theta_i={1\over N},\qquad
 \mathbb E(\Theta_i\Theta_j)={1\over N(N+1)}\quad(i\ne j).
\]
Thus the coefficient of \(x_ix_j\), \(i\ne j\), in (2) is
\[
 {2\over N^2}-{2\over N(N+1)}
 ={2\over N^2(N+1)}>0.
\]
Equivalently,
\[
 \boxed{\displaystyle
 \partial_i\partial_j\phi_N(1,\ldots,1)
 ={2\over N^2(N+1)}>0\qquad(i\ne j).}                      \tag{3}
\]
This contradicts (1).

For \(N=2\), the same obstruction is visible in the explicit formula
\[
 \phi_2(a,b)={b-a\over\log(b/a)},
\]
the logarithmic mean: its mixed derivative at \((1,1)\) is \(1/6>0\).

Consequently the Cartesian kernel has no representation as the potential
of an ordinary positive convolution semigroup of \(\mathbb R_+^N\)
subordinators when \(N\ge2\).  The one-dimensional exception is
\(F_1(a)=a^{-1}\), \(\phi_1(a)=a\), the pure-drift case.  This only rules
out the stated potential representation; it does not address non-Lévy,
signed, or nonconvolution auxiliary representations.

## The same obstruction for every finite \(c\), with full labels

Sum the complete labelled kernel over labels, equivalently take the
trivial character of the label group.  After removal of the common
exponential factor, its scalar density on \(z_i\ge0\) is
\[
 k_c(z)={\prod_{i=1}^N(1+{\bf1}_{\{z_i\ge c\}})\over H_{N,c}(\sum z_i)}.
\]
Its radial shell mass is exactly \(ds\).  If \(F_c\) denotes its Laplace
transform, then on the diagonal \(a=A{\bf1}\),
\[
 F_c(A{\bf1})={1\over A},\qquad
 \partial_iF_c(A{\bf1})=-{1\over NA^2}.                    \tag{4}
\]

Conditional on \(S=s\), the coordinates are exchangeable, have mean
\(s/N\), and satisfy
\(\operatorname {Cov}(Z_i,Z_j\mid S)=-\operatorname {Var}(Z_i\mid
S)/(N-1)\) for \(i\ne j\).  Hence
\[
 \partial_{ij}F_c(A{\bf1})
 ={2\over N^2A^3}
 -{1\over N-1}\int_0^\infty e^{-As}
   \operatorname {Var}(Z_i\mid S=s)\,ds.                   \tag{5}
\]
For \(\phi_c=F_c^{-1}\), substitution of (4)--(5) into the inverse
derivative identity gives the exact positive formula
\[
 \boxed{\displaystyle
 \partial_{ij}\phi_c(A{\bf1})
 ={A^2\over N-1}\int_0^\infty e^{-As}
 \operatorname {Var}(Z_i\mid S=s)\,ds>0,\qquad i\ne j.}    \tag{6}
\]
The strict inequality holds for every finite \(c>0\): on \(0<s<c\) the
conditional law is the uniform simplex law, whose coordinate variance is
\[
 {s^2(N-1)\over N^2(N+1)}.
\]
In particular,
\[
 \partial_{ij}\phi_c(A{\bf1})
 \ge {2-e^{-Ac}\{(Ac)^2+2Ac+2\}\over N^2(N+1)A}>0.          \tag{7}
\]
This violates the necessary Bernstein sign
\(\partial_{ij}\phi_c\le0\), so no finite-\(c\) scalar projection is a
subordinator potential either.

If a proposed labelled convolution semigroup uses the natural joint label
group and is invariant under that group convolution, its trivial-character
projection is a positive scalar convolution semigroup.  Its potential
would be precisely the summed-label kernel above, contradicting (6).
Thus retaining the complete labels cannot rescue this particular
convolution-semigroup representation.

There is also a direct probabilistic obstruction.  The radial projection
\(ds\) would force the total-coordinate Laplace exponent of a hypothetical
subordinator to be \(A\), hence
\(\sum_i Z_i(t)=t\) almost surely.  For \(m\) independent increments,
\[
 0\le Z_i(t/m)\le t/m,\qquad
 \operatorname {Var}Z_i(t)
 =m\operatorname {Var}Z_i(t/m)\le {t^2\over m}.
\]
Letting \(m\to\infty\) makes every coordinate deterministic, hence fixes
the angular direction.  This is incompatible with the nonzero uniform
simplex angular variance on the small-\(s\) shells.  It independently
explains the sign obstruction in (6).
