# Asymptotics of the capped fixed-corner coefficient

## The result

Fix \(A>1\). Let \(q_{N,A}\) solve the capped mean equation of
`corner82.md`. Let \(z=z_N\) be defined by
\[
 \overline\Phi(z)={z\over\sqrt N},\qquad r=\sqrt N\phi(z).
\tag{1}
\]
Then
\[
 q_{N,A}=N-\sqrt N z+{z^2\over3}+{4\over3}+o(1).        \tag{2}
\]
For the exact threshold root, uniformly for \(0\le\tau\le N^{-2}\),
\[
 q_{N,A}(e^{-\tau})-q_{N,A}(1)
 ={N\tau\over r+1}(1+o(1)).                             \tag{3}
\]
The cap is invisible at this precision for every fixed \(A>1\), but it is
essential for the fixed-corner theorem itself.

## 1. Gamma input without a black-box Edgeworth theorem

For \(S\sim\operatorname{Gamma}(N,1)\) and \(X=(S-N)/\sqrt N\), Stirling
and direct Taylor expansion of the density give
\[
 g_N(x)=\phi(x)\left(1+{x^3/3-x\over\sqrt N}\right)+R_N(x),
 \qquad \int|R_N|=O(N^{-1}).                            \tag{4}
\]
To prove (4), expand the logarithm on \(|x|\le N^{1/12}\); the integrated
Taylor remainder has the pointwise bound \(C N^{-1}(1+|x|^6)\phi(x)\), hence its integral is \(O(N^{-1})\). The gamma Chernoff bound and Gaussian
tail make the complement exponentially small. Integrating the correction
in (4) gives, uniformly for \(u=O(z^2)\),
\[
 \mathbb P\{S\le N-\sqrt Nz+u\}
 ={z\over\sqrt N}+{r\over N}\left(u+{1-z^2\over3}\right)
 +o(z^2/N).                                              \tag{5}
\]
Elementary Mills expansion in (1) gives
\[
 r=z^2+1+O(z^{-2}),\qquad z^2\asymp\log N.              \tag{6}
\]

## 2. Expansion of the capped mean

Put \(q=N-\sqrt Nz+u\), \(u=O(z^2)\), and \(V=(S-q)_+\). The lower branch
of the profile contributes exactly
\[
 L(q)=\mathbb E[e^{S-q}{\bf1}_{\{S\le q\}}]
 ={e^{-q}q^N\over N!}={r\over N}+o(z^2/N).              \tag{7}
\]
The last equality follows from \(L=(q/N)f_S(q)\) and (4).

For the middle branch, direct inversion at one gives, uniformly for
\(0\le v\le N^{2/3}\),
\[
 \Phi_N^{-1}(e^v)=1+{v\over N}
 +\left(1-{1\over2N}\right){v^2\over N^2}
 +O_A(v^3/N^3).                                         \tag{8}
\]
The cap begins at \(v=\log\Phi_N(A)=c_A N+O(1)\), where
\(c_A=\log(1+\log A)>0\). Gamma Chernoff bounds show that
\(V>N^{2/3}\), including the capped region, contributes \(o(z^2/N)\) to
the mean. Hence (8) can be integrated.

Let \(h=(N-1)/q-1\asymp z/\sqrt N\). Since
\(f_S(s)\le f_S(q)e^{-h(q-s)}\) for \(s\le q\),
\[
 \mathbb E(q-S)_+=O(1),\qquad
 \mathbb E(q-S)_+^2=O(\sqrt N/z).
\]
It follows that
\[
 {\mathbb EV\over N}={N-q\over N}+O(N^{-1}),\qquad
 {\mathbb EV^2\over N^2}={N+(N-q)^2\over N^2}+o(z^2/N).
\tag{9}
\]
Combining (5), (7)--(9) in \(\mathbb Et_q(S)=1\) gives
\[
 \mathbb Et_q(S)-1
 ={-(r+1)u+r(z^2-1)/3+r+z^2\over N}+o(z^2/N).           \tag{10}
\]
Terms \(O(N^{-1})\), including the final \(+1\) from the second moment,
are absorbed because \(N^{-1}=o(z^2/N)\).

## 3. Root slope and conclusion

The exact derivative formula in `corner82.md` and the same estimates yield
uniformly for \(u=O(z^2)\),
\[
 -{d\over dq}\mathbb Et_q(S)={r+1\over N}(1+o(1)).      \tag{11}
\]
The lower branch contributes (7), the middle branch contributes
\(N^{-1}(1+o(1))\), and the cap tail is negligible. Equation (10) first
localizes the root to \(u=O(z^2)\), then (10)--(11) give
\[
 u={r(z^2-1)/3+r+z^2\over r+1}+o(1)
 ={z^2\over3}+{4\over3}+o(1),                           \tag{12}
\]
by (6). The \(o(z^2/N)\) error in (10) is essential: divided by the root
slope it is \(o(1)\), so (2) is a genuine constant-order expansion.

For \(\sigma=e^{-\tau}\), the derivative formula remains uniform in the
resulting \(o(1)\)-neighborhood of the root. Since
\(e^{-\tau}-1=-\tau+O(\tau^2)\), integration of the inverse derivative
gives (3).

These results concern only the fixed-corner tensor-Hardy coefficient. They
do not provide variable-corner source sharing or a new general
strong-maximal endpoint bound.
