# 014 — \(\alpha=3/2\) transfer interface and isotropic-family audit

**Scope.** This is one independent, theory-only audit of the alternative
\(\alpha=3/2\) transfer chain.  It reads (but does not modify)
`alternative_alpha/transfer.md`; it does not treat any QMC value of \(g\) or
\(I\) as a certificate, and it does not import an \(\alpha=1\) small-time
local estimate.  The final part separately checks the algebraic
isotropic-family exclusion.  It is not a statement about the actual \(C_5\).

## Inputs and finite check

| input | SHA-256 |
|---|---|
| `Math-Research-cube5-lower/.../alternative_alpha/transfer.md` | `37688656ae9613578fa38cd2f7ec5d43896d1b083dbe995257569fafaa48c2da` |
| `Math-Research-cube-dim-5/.../round004/isotropic_exclusion.md` | `a26856aa716fa9806d16bcc3b37729b0a20cd87185b472a5520065207441b87f` |

One finite Python/Fraction check was run in the existing mathematical venv.
It verified the sixth-power comparison using the rational lower bound
\(\pi>333/106\):

```text
pi_lower_bound_validates_B2_gt_32= True
exact_margin_numerator=39582329299713887752981216536367920753269
exact_margin_denominator=8439763710605736088843734482944
```

No QMC, heat-kernel, or inherited \(\alpha=1\) local-time computation was
run or used.

## \(\alpha=3/2\) interface audit

Let \(V=\{(1,\pm1,\ldots,\pm1)\}\), with 16 elements, and

\[
 L=-\sum_{v\in V}|v\mathbin\cdot D|^{3/2},\qquad
 \widehat{p_t}(\xi)=e^{-t\psi(\xi)},\quad
 \psi(\xi)=\sum_{v\in V}|v\mathbin\cdot\xi|^{3/2}.
\]

The normalization is coherent: \(L\) has Fourier multiplier \(-\psi\),
and spanning of \(V\), compactness on the unit sphere, and homogeneity give
\(\psi(\xi)\asymp |\xi|^{3/2}\).  Therefore the operator domain is
\(H^{3/2}(\mathbb R^5)\) and its closed quadratic-form domain is
\(H^{3/4}(\mathbb R^5)\), with the usual inhomogeneous \(L^2\) addition.
This is the correct distinction; neither occurrence may be replaced by the
other.  The Dirichlet-form sign and the directional second-difference formula
are consistent with \(L\le0\).

For a cutoff \(\chi_R(x)=\chi(x/R)\), its exact homogeneous scaling gives

\[
 \|L\chi_R\|_\infty=R^{-3/2}\|L\chi\|_\infty.
\]

Likewise \(\|L\chi_\varepsilon\|_1=
\varepsilon^{5-3/2}\|L\chi\|_1=
\varepsilon^{7/2}\|L\chi\|_1\).  These are the two required nonlocal
boundary estimates.  They suffice respectively for recovering the obstacle
mass after \(R\to\infty\), and for proving \(\mu(\{0\})=0\) after
\(\varepsilon\downarrow0\).  Thus the mass-one positive measure used for
the all-scale step is supported away from the origin in the needed sense.

The only small repair needed here is a stated cutoff choice: the sentence
using \(\chi_R\uparrow1\) needs \(\chi\) to be *radial nonincreasing*,
with \(0\le\chi\le1\) and \(\chi(0)=1\), not merely smooth and radial.
Alternatively, one may choose an explicit nested cutoff family.  This is a
specification repair, not a changed estimate.  When compactness of the
obstacle is invoked, it should also say explicitly that
\((V_t-\kappa t)_+\) is compactly supported because \(V_t\in C_0\).

The global stable scaling, not a local small-time expansion, yields
\(p_t(0)=t^{-10/3}p_1(0)\).  Combining it with the \(L^1\to L^\infty\)
bound gives

\[
 V_t\le\min\{F,\,C m t^{-10/3}\},\qquad
 \int_0^\infty\min\{F,Cm t^{-10/3}\}\,dt
 ={10\over7}(Cm)^{3/10}F^{7/10}.
\]

Consequently the stated obstacle energy, frequency-tail, and mass-recovery
arguments have the correct exponents.  This part remains conditional only
on the explicitly supplied positive-source/obstacle construction, not on a
numeric candidate.

For the Green kernel, homogeneity is
\(G(rx)=r^{-7/2}G(x)\), so its critical weak exponent is
\(q=5/(7/2)=10/7\).  Hence \((G-c)_+\in L^1\) for every \(c>0\), which is
the needed finite-source fact.  The face-layer calculation gives

\[
 \| (G-c)_+\|_1={7\over3}\,32 I\,c^{-3/7},
 \qquad
 M(c)={7\over3}\,{I c^{-3/7}\over g-c}.
\]

Differentiating \(\log M=-\tfrac37\log c-\log(g-c)+\mathrm{const}\)
gives the unique minimum \(c=3g/10\).  Substitution gives, exactly,

\[
 M_{\min}=\left({10\over3}\right)^{10/7}I g^{-10/7}.
\]

Thus the requested \((10/3)^{10/7}\) factor is correct (and is not
\((10/3)^{3/7}\) or a factor inherited from the \(\alpha=1\) formula).

Finally, the Fourier scaling identity

\[
 K_r*Lu={\mu_r*u-u\over N r^{3/2}}
\]
has the correct sign and power: evaluate
\(\widehat\mu-1=-N\psi\widehat K\) at \(r\xi\).
On \(D^c\), positivity of \(\mu\) then gives the fixed-scale weak
estimate.  Taking a common null set over positive rational scales, using
strong \(L^1\) continuity of dilations for bounded compactly supported
sources, and then monotone truncation of arbitrary \(L^1\) sources is a
valid route to all positive scales.  This establishes the general-input,
all-scale interface conditional on the two numerical/analytic inputs
\(g>0\) and finite \(I\); it supplies neither one.

## Isotropic stable positive-part family

The second input is correct as the explicitly scoped algebraic exclusion.
For \(0<\alpha\le2\), the radial Green model has
\(q=5/(5-\alpha)\).  Radial projection of each signed cube-face orthant
partitions \(S^4\) into \(10\cdot16=160\) congruent pieces.  Since
\(|S^4|=8\pi^2/3\), the exact shape integral is

\[
 J=\int_{[0,1]^4}(1+|y|^2)^{-5/2}\,dy={\pi^2\over60}.
\]

For the already-displayed positive-part cost formula its minimizer is
\(c=\alpha g/5\); substituting the best essential cube bound produces

\[
 B(\alpha)=\left({5\over\alpha}\right)^{5/(5-\alpha)}
             {5^{5/2}\pi^2\over60}.
\]

Its logarithmic derivative is
\[
 {5\over(5-\alpha)^2}
 \left[\log(5/\alpha)-(5/\alpha-1)\right]<0,
\]
so the minimum on \((0,2]\) is \(B(2)\).  The Fraction run above proves
strictly \(B(2)>32\) by comparing sixth powers with \(\pi>333/106\).

This last result only excludes improvement through that isotropic Green,
positive-part, finite-source **cost formula**.  It does not establish a
lower bound for \(C_5\), exclude anisotropic models, or validate a general
\(\alpha\) transfer theorem beyond the displayed algebra.

## Audit disposition

**PASS subject to the explicit nested-cutoff specification above.**  The
\(\alpha=3/2\) transfer contains the required generator, Sobolev-domain,
truncation, finite-source, and all-scale interfaces.  Any eventual numeric
upper bound remains conditional until \(g\) and \(I\) receive their own
rigorous certificates; QMC is not one.  No unproved \(\alpha=1\)
small-time local exponent was used in this audit.
