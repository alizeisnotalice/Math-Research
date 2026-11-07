# Arbitrary capped terminals: budget barriers replace monotonicity

## Statement and scope

Use the exact notation, independent Gamma–binomial blocks, canonical d_N,
and parameter box of r062/log_blocks62.md. Let N=rb, 1≤r≤C log N,
L=log N, ρ=√(L/N), |c_i−μ|≤K√(L/b), and

\[
|s-(N\mu-d_N)|\le C_0rL.
\]

Let H:R_+^r→[0,A] be **any fixed measurable terminal**, with fixed A>1,
independent of the actions, and assume

\[
J(H)=\int e^{-S(t)}\Phi_N(H(t))dt\le e^{-s},\qquad
E_{c_*}H\ge1
\]

at one vector c_* in the box. No monotonicity, KKT condition, angular
smoothness, or finite-layer hypothesis is imposed on H. Then, uniformly
throughout the box,

\[
\boxed{\nabla_c^2G(c)\preceq-(p_\mu-o(1))\frac{L}{r}I_r,\qquad
G(c)=E_cH.}                                                     \tag{1}
\]

The fixed-r common-location response also remains
G_d=(r+o(1))L/N. In particular the two-block common second derivative is
−p_μ L(1+o(1)). This is a theorem about the stated terminal model, not a
physical-space truncation or a representation theorem for arbitrary
noncentered rectangles.

## 1. A small exceptional set from the actual modular budget

Write F_c(t)=∏_i f_{b,c_i}(t_i). A Gamma(b,1) density has supremum at
most one; each block density is a convex mixture of its translates.
Therefore F_c≤1 for every parameter in the box (in fact for every finite
nonnegative action vector).

Fix K_big≥20, increasing it once if necessary to absorb all desired
polynomial errors. The sets below depend on H and s but not on c. Define

\[
B=\{H>1:\Phi_N(H)>e^{S-s+K_{\rm big}L}\}.
\]

Since ∫_B e^(−S)Φ_N(H)>e^(−s)N^K_big |B|, the budget gives

\[
|B|\le N^{-K_{\rm big}},\qquad
\sup_cP_c(B)\le N^{-K_{\rm big}}.                              \tag{2}
\]

The first expression is Lebesgue measure, the second uses F_c≤1. This is
not a pointwise envelope for H on all of its domain.

Set I=1_{S≤s−K_big L}. Since Φ_N(h)≥h for h≥0,

\[
E_c[HI]\le\int_{S\le s-K_{\rm big}L}Hdt
\le e^{s-K_{\rm big}L}J(H)\le N^{-K_{\rm big}}.                 \tag{3}
\]

Thus the lower-core mass estimate also needs no pointwise ordering.

## 2. The scalar envelope and its moments

Let Ψ_N denote the inverse of the strictly increasing Φ_N on [0,∞), and
put

\[
U(S)=\min\{A-1,[\Psi_N(e^{S-s+K_{\rm big}L})-1]_+\}.
                                                                    \tag{4}
\]

U is a bounded nondecreasing scalar function. For P=(H−1)_+, off B we
have P≤U; on B both functions lie in [0,A−1]. Therefore

\[
E(P-U)_+\le(A-1)N^{-K_{\rm big}},\qquad
EP^2\le EU^2+(A-1)^2N^{-K_{\rm big}}.                         \tag{5}
\]

Let ε_N=C[L^−1+rL^(3/2)/√N+√rρ]. The same central box as r062,
max_i |T_i−bm(c_i)|≤B_0√(bL), gives

\[
EU\le\frac{d_N}{N}(1+\varepsilon_N),\qquad
EU^2\le\left(\frac{d_N}{N}\right)^2(1+\varepsilon_N).           \tag{6}
\]

Here are the moment details, which do not use H. Set
δ=(S−s+K_big L)_+ and a=ES−s+K_big L=d_N+O(rL)>0. On the central box,
δ/N=O(√rρ)=o(1). For h=1+y>1,
log Φ_N(1+y)≥(N−1)y/(1+y), so (4) implies
U≤δ/(N−1−δ). Off the box use U≤A−1 and its arbitrarily high polynomial
tail. Since Var S=O(N),

\[
E\delta\le a+\operatorname{Var}(S)/a,\qquad
E\delta^2\le a^2+\operatorname{Var}(S).
\]

The first bound follows from (−x−a)_+≤x²/a for centered x; the second
uses δ²≤(S−s+K_big L)². Expanding the central denominator yields (6);
N/d_N²=O(1/L), and the slack rL/d_N is covered by ε_N. The cap handles
the complement, not an extrapolation of the denominator bound there.

The already-proved triangular-array moderate deviation lemma of r062 gives

\[
F_0:=P_c(I=1)=\frac{d_N}{N}(1+O(\varepsilon_N)).                \tag{7}
\]

## 3. Near reception forces both positive and negative parts

Suppose temporarily G(c)≥1−η, where η=o(ρ). With D=(1−H)_+,
ED≤EP+η, while (3) gives

\[
ED\ge F_0-N^{-K_{\rm big}},\qquad
E[DI]\ge F_0-N^{-K_{\rm big}}.
\]

Combining (5)–(7) proves

\[
\|D-I\|_1=O(\rho\varepsilon_N+\eta),\qquad
\|H-1\|_2=O(\sqrt\rho).                                      \tag{8}
\]

We also have EP≥ED−η≥F_0−N^−K_big−η. The identity

\[
E|P-U|=EU-EP+2E(P-U)_+
\]

therefore gives the crucial replacement for monotonicity:

\[
\boxed{\|P-U\|_1=O(\rho\varepsilon_N+\eta).}                  \tag{9}
\]

No assertion is made that P itself is monotone, or that P=U pointwise.

## 4. Lower coordinate location response

Let ℓ_i=f'_{b,c_i}/f_{b,c_i}. The ordinary Gamma mixture gives
Eℓ_i²≤1/(b−2), and on the per-block central box |ℓ_i|≤C√rρ.
The translated expectation is differentiable against bounded H by the
C²(L1)/W¹,¹ mixture lemma in r062. Its location response is

\[
G_{d_i}=E[D\ell_i]-E[P\ell_i].
\]

Only U needs an order property. Its scalar monotonicity and the vanishing
axis trace give E[Uℓ_i]≤0 by Stieltjes integration by parts. Hence

\[
G_{d_i}\ge E[I\ell_i]
-E[|D-I||\ell_i|]-E[|P-U||\ell_i|].                            \tag{10}
\]

On the box use (8)–(9) and the score supremum. Off the box use the cap,
Fisher Cauchy–Schwarz and the arbitrary polynomial tail. The exact identity
E[Iℓ_i]=f_{S,c}(s−K_big L) then yields

\[
G_{d_i}\ge\frac LN\left[1-o(1)
-O\left(\frac{\sqrt r\eta}{\rho}\right)\right],              \tag{11}
\]

where the other relative errors are bounded by
C[√r/L+r^(3/2)L^(3/2)/√N+rρ]=o(1) for r≤C log N. Equation (10) is the
only growing-block step in the former proof that used H monotonicity.

## 5. Continuation and Hessian

Equations (5)–(8) give E|H−1|=O(ρ) in the near-receiving regime. The
ordinary parameter score bounds from r062 thus still give
|∂_{c_i}G|≤CρL. The same first-exit argument along each segment from c_*
propagates reception to the entire box with

\[
G(c)\ge1-O(r^{3/2}L^2/N).
\]

The final η makes √rη/ρ=O(r²L^(3/2)/√N)=o(1). Thus (11) is uniform
throughout the box.

Finally the ordinary kernel identities of r062 are valid for any bounded
measurable H:

\[
G_{c_ic_i}=bm''(c_i)G_{d_i}+E[(H-1)e_i],\qquad
G_{c_ic_j}=E[(H-1)s_is_j]\quad(i\ne j).
\]

The L2 bounds ∥e_i∥₂≤CL², ∥s_i∥₂≤CL, independence, and (8) bound every
error entry by C√ρ L². Since m''(c_i)=−p_μ+o(1), Gershgorin gives (1)
with relative spectral error Cr²√ρ L=o(1). For fixed r the earlier
Cauchy bound on EPℓ gives the full, two-sided common-location asymptotic,
so the two-block common-curvature equality also survives.

## What is still missing

The terminal must remain one fixed action-independent H with the stated
actual block modular budget. This theorem removes its monotonicity
assumption; it does not prove that an arbitrary physical rectangle family
has such a representation, nor control the extra angular-kernel transport
terms when the effective conditional terminal changes with c. The general
C_N coefficient bracket is unchanged.
