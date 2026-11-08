# r074: nested physical support obstruction

## Result

This is a valid and stronger obstruction to a *pointwise* transfer from an
actual, locally efficient action to a deterministic reference action.  It
uses neither a posterior moderate-deviation event nor a conditional-label
comparison.  It does not contradict a theorem with a global packing or an
explicit common-source allocation condition.

Let `L=log N`, let mu solve mu=1+exp(−mu), set c_r=mu−sqrt(mu)/2 sqrt(L/N), and choose critical actions

\[
 c_a=c_r+\Delta c,\qquad \Delta c={L^{1/4}\over\sqrt N}>0,       \tag{1}
\]

and work at the receiver `x=(1,...,1)`.  Put

\[
 R_c=[-e^{-c},1]^N,\qquad V_c=|R_c|=(1+e^{-c})^N.                  \tag{2}
\]

Then `R_{c_a}` is a strict subset of `R_{c_r}`.  With
`S(y)=-sum_i log|y_i|` away from coordinate hyperplanes, take a cap-`3/2`
scalar KKT profile `h_a=h_{*,3/2,c_a}` and define the single physical source

\[
 f(y)=h_a(S(y))\,\mathbf1_{R_{c_a}}(y).                         \tag{3}
\]

The same function is used for both actions.  In sign/log coordinates,
`R_c` means `z_i>=c` on a negative coordinate and `z_i>=0` on a positive
one.  Thus at `c_a` the conditional action kernel sees exactly the scalar
profile `h_a(s)` on every allowed label, with no label averaging error.

## Exact rectangle identity

For the normalized physical action average at this receiver,

\[
 A_c f:={1\over V_c}\int_{R_c}f(y)\,dy,
\]

the scalar normalization gives `A_{c_a}f=E_{c_a}h_a(S)=1`.  Since the
support in (3) is contained in `R_{c_a}`, nesting yields the exact identity

\[
 A_{c_r}f={1\over V_{c_r}}\int_{R_{c_a}}f
 ={V_{c_a}\over V_{c_r}}A_{c_a}f.                            \tag{4}
\]

No conditional Fourier comparison is being used in (4).  Expanding the
volume ratio gives, with `p_c=(1+e^c)^{-1}`,

\[
 \log {V_{c_a}\over V_{c_r}}
 =-Np_{c_r}\Delta c+O(N\Delta c^2)
 =-p_{c_r}\sqrt N\,L^{1/4}+O(L^{3/4})\longrightarrow-\infty. \tag{5}
\]

Consequently the reference average tends to zero, not merely to
`1-O(L/N)`.  This is legitimate despite the two actions having different
allowed sign supports: the physical source is fixed, and the reference
rectangle simply contains the source rectangle.

## Local scalar cost and strict compact version

For the local cost operator at the actual action, coarea/conditional
disintegration gives directly

\[
 j_{c_a}\Phi_N(f)(1)=\int_0^\infty e^{-s}\Phi_N(h_a(s))\,ds
 =J_{*,3/2}(c_a).                                             \tag{6}
\]

The uniform scalar calibration in the critical compact action window gives

\[
 \log {J_{*,3/2}(c_a)\over J_{*,2}(c_r)}=o(L),              \tag{7}
\]

which is all that is needed for `J_{*,3/2}(c_a)=J_{*,2}(c_r)N^{o(1)}`.
The quadratic model term is `O(D^2L)=O(sqrt L)` for
`D=Delta c/sqrt(L/N)=L^{-1/4}`, but the available uniform statement should
only be quoted as (7), because its residual is merely `o(L)`.

The scalar calibration used in (7) is the compact-parameter statement

\[
 \log J_{*,A}(\mu+t\sqrt{L/N})=-q_0-\kappa(t)L+o(L),\qquad
 \kappa(t)=\alpha-\frac p2(t-t_*)^2,
 \quad t_*=-\sqrt\mu/2,
\]

uniform for t in a fixed compact set and each fixed A>1. The uniform platform lower estimate is r019 (4), with its bounded-t qualification; the fixed-action specialization of r022 (12),(14) gives the matching upper, and r069's fixed-cap repair is uniform on this shrinking action window. At t_a=t_*+L^{-1/4}, the explicit difference is (p/2)sqrt L plus o(L), hence (7). Here p=(mu−1)/mu, alpha=(9−mu)/24, and q_0 is the calibrated critical location Nmu−sqrt(mu N)z_N from r069. This imported equality is only at logarithmic precision. It is not an assertion of an O(sqrt L) bound for the full remainder.

Use instead the r072 wide sum window `W_N`.  Put
\(\eta_a=E_{c_a}[h_a(S)\mathbf1_{S\notin W_N}]\).  With A_0=3/2, its uniform Cramer
bound gives `eta_a<=2A_0N^{-6}`.  The strict compact source is

\[
 F(y)=\tau h_a(S(y))\mathbf1_{W_N}(S(y))\mathbf1_{R_{c_a}}(y),
 \qquad \tau={1+N^{-3}\over1-\eta_a}=1+O(N^{-3}).        \tag{8}
\]

It has actual average exactly `1+N^{-3}`, cap below `2`, and
`tau^N=1+o(1)`, so the actual local cost is at most `(1+o(1))J_{*,3/2}(c_a)`.  Since `S in W_N`, it is
compact away from the axes.  Identity (4) survives exactly, because it uses
only support containment in `R_{c_a}`.

There is a stronger conditional-cost identity.  On `W_N`, r065 section 3
gives `|partial_c log p_c(s)|<=CL`; hence

\[
 {H_{N,c_a}(s)\over H_{N,c_r}(s)}
 ={Z_{c_a}\over Z_{c_r}}{p_{c_a}(s)\over p_{c_r}(s)}
 ={V_{c_a}\over V_{c_r}}\,(1+o(1)),                       \tag{9}
\]

uniformly, because `Delta c L=o(1)`.  The numerator of the reference
conditional kernel is the `c_a` numerator, while its denominator is
`H_{N,c_r}`.  Therefore, for every band contained in `W_N`,

\[
 j^{\rm band}_{c_r}\Phi_N(F)(1)
 ={V_{c_a}\over V_{c_r}}(1+o(1))
   j^{\rm band}_{c_a}\Phi_N(F)(1).                        \tag{10}
\]

Thus not only the reference average but also its local/band cost is
exponentially (indeed superpolynomially) suppressed; the actual captured
band is preserved up to a factor 1+o(1).  This makes clear that no posterior-tail lemma is needed.

## Scope

The example refutes an inference of the form

> at a receiver, success plus `j_{c_a}Phi(f)<=J_*N^epsilon` at its actual
> critical action forces near-success for one predetermined `c_r`.

It does not refute a global result that additionally controls how such
action-matched nested source packets can occur across receivers.  In
particular it supplies no global `C_N` counterexample by itself.  Its lesson
is that a proof requiring a common deterministic reference needs a
common-source/action allocation or packing statement; local scalar cost
qualification alone cannot furnish that transfer.
