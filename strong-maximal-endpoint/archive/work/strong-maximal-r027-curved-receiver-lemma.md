# r027: a strict curved-receiver band for symmetric two-block profiles

Date: 2026-11-06.  This is an actual anchored-block rectangle statement.  It
does not bound the full strong maximal halo or prove that a block profile
beats \(Q_n\).

Let \(n=2b\), take identical coordinate parameter \(c\) in both blocks, and
let \(f_b\) be the density of a block sum \(T_b\).  Let
\(H:[0,\infty)^2\to[0,A]\) be symmetric under coordinate interchange and
coordinatewise nondecreasing.  Define the actual block-kernel average
\[
 B(u_1,u_2)=\mathbb E H(u_1+T_1,u_2+T_2),
\]
where \(T_1,T_2\) are independent copies of \(T_b\).  In coordinates
\(u=(d+v,d-v)\), write \(B(d,v)=B(d+v,d-v)\).

## Gamma-mixture translation bound

The block density is a probability mixture of translates of
\(g_b(y)=y^{b-1}e^{-y}/\Gamma(b)\).  Hence its first two \(L^1\) derivative
bounds are no larger than those of \(g_b\).  For \(b>2\), putting
\(U=(b-1)/Y-1\), \(Y\sim\Gamma(b,1)\),
\[
 \|g_b'\|_1\le(\mathbb EU^2)^{1/2}=\frac1{\sqrt{b-2}},
\]
\[
 \|g_b''\|_1
 \le \mathbb EU^2+\mathbb E\frac{b-1}{Y^2}=\frac2{b-2}.\tag{1}
\]
Wolfram 15 independently evaluated both expectations as \(1/(b-2)\).
For the product translation kernel, differentiation along \((v,-v)\) gives
\[
 \|\partial_v^2[f_b(\cdot-v)f_b(\cdot+v)]\|_1
 \le2\|f_b''\|_1+2\|f_b'\|_1^2\le\frac6{b-2}.\tag{2}
\]
Mixture shifts do not alter these bounds.

Symmetry of \(H\) makes \(B(d,v)\) even in \(v\).  Taylor's theorem in
convolution form and (2) therefore yield the uniform curvature estimate
\[
 |B(d,v)-B(d,0)|\le \frac{3A}{b-2}v^2.\tag{3}
\]

## Actual curved halo consequence

Suppose for some \(d\),
\[
 B(d,0)\ge1+\eta\qquad(0<\eta\le1).\tag{4}
\]
Set \(R=\sqrt{\eta(b-2)/(6A)}\).  Equation (3) gives
\(B(d,v)>1\) for \(|v|<R\).  Coordinate monotonicity then gives the real
anchored-rectangle witness region
\[
 \{(u_1,u_2)=(d'+v,d'-v):d'\ge d,\ |v|<R\}.
\]
Its exact reference measure is
\[
 \int 1_{\rm region}e^{-u_1-u_2}\,du_1du_2
 =2R e^{-2d}
 \ge c\,e^{-2d}\sqrt{\eta n/A}.\tag{5}
\]
Thus a symmetric cap-\(A\) block profile with a fixed central strict margin
has a true curved receiver band that contributes a \(\sqrt n\) transverse
volume factor, rather than only the corner orthant volume \(e^{-2d}\).

The r026 KKT profile at equal blocks and equal \(c\) is symmetric and
coordinatewise nondecreasing, so it satisfies the hypotheses (with \(A=2\)).
Because its kernel tends to the cap as \(d\to\infty\), a point satisfying
(4) exists for every fixed \(\eta<1\).  What is *not* established here is
that this point has \(d=d_0+O(1)\), where \(B(d_0,0)=1\).  That radial
location determines whether (5) exactly cancels the r026 corner loss and
whether the complete curved receiver score improves or merely matches the
scalar benchmark.  The present lemma is therefore a rigorous nonlinear
receiver mechanism, not a claim of an improvement over \(Q_n\).

It must not be used cheaply at a near-root margin \(\eta\to0\): its width is
only \(\sqrt{\eta n}\), and it gives no estimate for the radial displacement
needed to create that margin.

## Conditional Gaussian bridge: exact algebra and conditional conclusion

This paragraph is a quantitative target, not a proved Gamma--Bernoulli
theorem.  For equal blocks \(b=n/2\), condition at a lower-moderate total
value \(S=s\) and use the common negative tilt that centers it.  The normal
bridge has
\[
 W=\frac{2(T_1-s/2)}{\sqrt{\mu n}},\qquad
 k_s(T_1)=\sqrt{\frac2{\pi\mu}}n^{-1/2}e^{-W^2/2}(1+o(1)).\tag{6}
\]
The exact likelihood identity \(R_2(t,s-t)=R_n(s)k_s(t)\) then gives
\[
 \log R_2=\log R_n(s)-\frac12\log n+
 \log\sqrt{\frac2{\pi\mu}}-\frac12W^2+o(1).\tag{7}
\]
If the capped KKT response is uniformly linearizable throughout its active
strip, the first-order scalar/block budget matching shifts its log threshold
by \(\alpha=\mathbb E(W^2/2)=1/2\).  Under the additional uniform tail and
cost hypotheses needed to integrate this statement, one obtains
\[
 \frac{I_2}{I_1}\sim e^{1/2}\sqrt{\frac{\pi\mu n}{2}},
 \qquad
 V_2\sim\sqrt{\frac{\pi\mu n}{2}}\,V_{\rm corner},
 \qquad
 \frac{V_2/I_2}{Q^{\rm cap}_{n,c}}\to e^{-1/2}.\tag{8}
\]
Here \(Q^{\rm cap}_{n,c}\) denotes only the same-\(c\), height-capped
scalar benchmark, not \(Q_n\).  The derivation of (8) still requires a
uniform tilted local limit in the KKT active strip, response linearization at
the plateau/soft interfaces, and global conditional-tail control.  None is
claimed in this note; (6)--(8) merely state the precise route supported by
the current numerical values.

For comparison, the quadratic surrogate
\(\ell_b(t)=t-(t-b\mu)^2/(2\mu b)\) is completely solvable.  If the
response is an increasing function of
\(\ell_b(b\mu+d+v)+\ell_b(b\mu+d-v)\), a central boundary \(d_0=o(n)\)
has exact equation
\[
 d-\frac{d^2}{\mu n}=d_0-\frac{d_0^2}{\mu n}+\frac{v^2}{\mu n}.
\]
Thus \(d(v)=d_0+v^2/(\mu n)+o(1)\) for \(v=O(\sqrt n)\), and its receiver
volume is asymptotic to \(e^{-2d_0}\sqrt{\pi\mu n/2}\).  Wolfram 15 directly
integrated this Gaussian bridge and returned the normalized constant
\(\sqrt{\pi\mu/2}=1.4171123483\ldots\).  This exact surrogate is not a
replacement for the missing real-kernel estimates.
