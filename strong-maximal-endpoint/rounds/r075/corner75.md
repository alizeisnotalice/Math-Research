# A positive global charge for the affine fixed-corner policy family

## Result

The qualified affine policy
\[
 c_i(x)=\log{x_i\over r_a},\qquad r_a>0,
\tag{1}
\]
does not threaten the endpoint upper bound, even though its canonical
reanchor has the large weighted distortion identified in
`efficient_review75.md`.  Its actual rectangles are
\[
 [-e^{-c_i(x)}x_i,x_i]=[-r_a,x_i],                         \tag{2}
\]
so all receivers use one fixed negative corner.  They are controlled by a
translated one-sided tensor Hardy operator.  For every bounded compact
\(0\le F\le2\),
\[
 \bigl|\{x:H_aF(x)>1\}\bigr|
 \le \exp\!\left\{N+{1\over2\zeta}+O(N^{-1})\right\}
 \int_{\mathbb R^N}\Phi_N(F),
\tag{3}
\]
where
\[
 \zeta={\log(1+\log2)\over\log2}>0.
\tag{4}
\]
Since \(Q_{\rm ref}=\exp\{\mu N-o(N)\}\) and \(\mu>1\), this is
exponentially smaller than the critical endpoint scale.  Thus the entire
affine-corner family can be globally charged, independently of local
success or the local cost screen.

This is an application to the r075 bad-geometry family, not a new
fixed-corner theorem; the fixed-anchor theory of r028 is stronger.

## 1. Reduction to tensor Hardy

Put \(X_i=x_i+r_a\) and \(g(z)=F(z-r_a\mathbf1)\) on the positive
orthant.  Then
\[
 H_aF(x):={1\over\prod_iX_i}
 \int_{\prod_i[-r_a,x_i]}F(y)\,dy
 ={1\over\prod_iX_i}\int_{[0,X_1]\times\cdots\times[0,X_N]}g(z)\,dz
 =H_1\cdots H_Ng(X),                                     \tag{5}
\]
where \(H_i\) is the one-dimensional Hardy average in coordinate \(i\).
The translated positive-orthant modular is at most the full physical modular.
If (1) is additionally required
to have \(c_i\ge0\), then \(x_i\ge r_a\), which only restricts the
receiver set and leaves the upper bound unchanged.

## 2. The cap-two Orlicz-to-\(L^p\) comparison

For \(N\ge2\), set
\[
 p_N=1+\zeta(N-1)>1.                                     \tag{6}
\]
For \(0\le t\le1\), \(t^{p_N}\le t=\Phi_N(t)\).  For
\(1\le t\le2\), write \(u=\log t\in[0,\log2]\).  The desired bound
\(t^{p_N}\le\Phi_N(t)\) is equivalent to
\[
 \zeta u\le\log(1+u).                                   \tag{7}
\]
The difference \(\log(1+u)-\zeta u\) is concave and vanishes at both
endpoints of this interval, by (4).  It is therefore nonnegative.  Hence
\[
 \|g\|_{p_N}^{p_N}\le\int_{\mathbb R_+^N}\Phi_N(g)
 \le\int_{\mathbb R^N}\Phi_N(F).                         \tag{8}
\]

## 3. Hardy norm and weak charge

The one-dimensional Hardy inequality has norm \(p_N'=p_N/(p_N-1)\) on
\(L^{p_N}(\mathbb R_+)\).  Iterating it in the \(N\) coordinates and
using Chebyshev at height one gives
\[
 \begin{aligned}
 |\{H_aF>1\}|
 &\le\|H_1\cdots H_Ng\|_{p_N}^{p_N}\\
 &\le(p_N')^{Np_N}\|g\|_{p_N}^{p_N}\\
 &\le\left({p_N\over p_N-1}\right)^{Np_N}
       \int\Phi_N(F).
 \end{aligned}                                           \tag{9}
\]
With \(q=\zeta(N-1)\),
\[
 \log\left(\left({p_N\over p_N-1}\right)^{Np_N}\right)
 =N(q+1)\log(1+q^{-1})
 =N+{1\over2\zeta}+O(N^{-1}),                            \tag{10}
\]
which proves (3).  The calculation is valid for bounded compact inputs;
more generally it only needs the cap-two modular to be finite.

## 4. Finite corner families and scope

For \(M_N\) prescribed affine corners, the union bound multiplies (3) by
\(M_N\).  Hence a constant number of corners, or \(M_N=e^{o(N)}\), remains
\(\exp\{N+o(N)\}\) and is negligible versus \(Q_{\rm ref}\).  The same
holds as long as
\[
 \log M_N\le(\mu-1-\delta)N
\tag{11}
\]
for some fixed \(\delta>0\), at the leading exponential scale.

This does not discretize a general continuum of moving anchors: no automatic
net or union count is available for that class.  It only removes the
specific fixed-negative-corner affine policy family from the list of
possible qualified geometric obstructions.

## One fixed corner with all orientations

The same constant works if receivers lie in all 2^N orthants around one fixed corner a: reflect each source/receiver orthant to the positive orthant, apply tensor Hardy there, and add their L^p integrals. The source orthants are disjoint up to null coordinate hyperplanes, so their modulars sum to the one original modular. No factor 2^N is introduced. This concerns boxes with a as one vertex, not all boxes merely containing a.

For completeness, the one-dimensional Hardy estimate used above follows by integration by parts for a compact nonnegative f: with U(x)=integral_0^x f, integral (U/x)^p = p' integral (U/x)^(p−1) f. Holder gives norm Hf<=p' norm f. Boundary terms vanish for compact f, and truncation extends the inequality to L^p. Fubini then iterates it in the N coordinates.
