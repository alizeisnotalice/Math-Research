# r028: KKT threshold/soft decomposition and the precise bridge gap

Date: 2026-11-06.  This note isolates a rigorous part of the proposed
equal-block curved-receiver asymptotic.  It does not prove the claimed
\(e^{-1/2}\) limit for the real Gamma--Bernoulli kernel.

## Exact cap-two KKT calculus

Put \(P_n(h)=\Phi_n'(h)=(1+\log h)^{n-2}(n+\log h)\) for \(h>1\).  If
\(x=\log(\lambda R)\), the height-capped KKT response is
\[
 \Psi_n(x)=
 \begin{cases}
 0,&x<0,\\
 1,&0\le x\le\log n,\\
 P_n^{-1}(e^x),&\log n<x<\log P_n(2),\\
 2,&x\ge\log P_n(2).
 \end{cases}\tag{1}
\]
Write \(\Psi_n=1_{[0,\infty)}+G_n\).  Then \(G_n\) is globally
Lipschitz with constant \(C/n\), uniformly in \(n\ge4\).  Indeed Wolfram
15 differentiates to
\[
 \frac{d}{dh}\log P_n(h)=
 \frac{(n-1)(n-1+\log h)}{h(1+\log h)(n+\log h)}\ge c n
 \quad(1\le h\le2).\tag{2}
\]
Thus, for arbitrary real \(X,a\),
\[
 \left|\Psi_n(X+a)-\Psi_n(X)-
 \bigl(1_{\{X+a\ge0\}}-1_{\{X\ge0\}}\bigr)\right|
 \le C|a|/n.\tag{3}
\]
This proves rigorously that a bounded angular likelihood perturbation changes
the soft KKT response by only \(O(1/n)\); the leading response comes from the
moving indicator threshold.  Since the scalar root has
\(f_S(q)\asymp(\log n)/n\), the soft term is lower order by \(1/\log n\).

## What remains to turn the bridge into a theorem

For equal blocks \(b=n/2\), the exact identity
\[
 R_2(t,s-t)=R_n(s)k_s(t)\tag{4}
\]
has a rigorously available conditional local limit.  At the scalar KKT
threshold \(q\), uniformly for \(s-q=O(\log n)\) and
\(|W|\le C\sqrt{\log n}\), where
\(W=2(t-s/2)/\sqrt{\mu n}\), common tilting gives
\[
 k_s(t)=\sqrt{\frac2{\pi n v_s}}
 \exp\left(-\frac{\mu W^2}{2v_s}\right)
 \left(1+O_C\left(\frac{1+|W|^3}{\sqrt n}\right)\right),
\]
where \(v_s=\mu+o(1)\).  In particular,
\[
 \log k_s(t)=-\tfrac12\log n+
 \log\sqrt{\frac2{\pi\mu}}-\tfrac12W^2+o(1),\tag{5}
\]
The proof is the r016 Fourier calculation after common tilting: uniform
fourth-order characteristic expansion on the central arc and the Exp-factor
bound \(|\varphi(t)|\le(1+(t/(1+\lambda))^2)^{-1/2}\) off it.  Thus a
standalone conditional LLT is not the remaining gap.

The missing proposition is the uniform moving-boundary estimate for the
indicator part of (3):
\[
 \mathbb E[\Psi_n(\log R_2-\tau_2)]-
 \mathbb E[\Psi_n(\log R_n-\tau_1)]
 =f_S(q)\{\alpha-\tfrac12\}+o(f_S(q)),\tag{6}
\]
where \(\alpha=\tau_1-\tau_2-\tfrac12\log n+
\log\sqrt{2/(\pi\mu)}\).
Equation (3) supplies the soft error in (6); the unproved portion is the
uniform moving-boundary estimate, together with its cost and curve-root
counterparts.

More explicitly, (6) follows from (5) once the following standard tilted
density estimates are written with uniform constants:
\[
 (\log R_n)'(s)=1+o(1),\qquad f_S(s)=f_S(q)(1+o(1))
 \quad(|s-q|\le C\log n).\tag{7}
\]
They follow by differentiating the tilted density at its mean:
\(f_S=e^{\lambda s+K}g\), \(g'(s)/g(s)=O(1/n)\), and
\(\lambda\log n=o(1)\).  Then the indicator in (3), conditional on \(W\),
has threshold
\[
 s\ge q-\alpha+W^2/2+o(1).
\]
Its tail probability minus \(\mathbb P(S\ge q)\) is uniformly
\(f_S(q)(\alpha-W^2/2)+o(f_S(q))\) for \(|W|\le C\sqrt{\log n}\).
The tilted Bernstein bound for one block, divided by the total tilted local
density, makes the remaining conditional \(W\)-tail smaller than
\(o(f_S(q))\) after taking \(C\) sufficiently large.  Gaussian integration
then gives (6), since \(\mathbb EW^2=1+o(1)\).

This reduces the real proof burden further: the moving *mean* constraint is
a local tilted-density calculation.  The still-unclosed steps are the
matching asymptotic for the Orlicz cost and the receiver boundary uniformly
out to \(|v|\asymp\sqrt{n\log n}\), required for the volume/UI conclusion.

If (5)--(6), their analogous cost statement, and the uniform boundary form
\(d(\sqrt n y)=y^2/\mu+o(1)\) are proved, then
\[
 \alpha\to\tfrac12,qquad
 I_2/I_1\sim e^{1/2}\sqrt{\pi\mu n/2},
 \qquad V_2\sim\sqrt{\pi\mu n/2}\,V_{\rm corner},
\]
and the same-\(c\), cap-two full-receiver score ratio tends to \(e^{-1/2}\).
The log-concavity theorem from r027 gives monotonicity of \(d(v)\), useful
for the remaining uniform-integrability tail step, but it does not supply
(5)--(6) or the quadratic lower growth by itself.
