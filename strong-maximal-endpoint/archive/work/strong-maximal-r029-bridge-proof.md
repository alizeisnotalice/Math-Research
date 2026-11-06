# r029 Poisson-corrector: conditional bridge interface

Date: 2026-11-06.  This proves the bridge component requested for the
corrector, conditional on the standard uniform local Fourier lemma stated
below.  It does not prove the separate soft-terminal or supersolution
interfaces.

## Set-up and tilt

Write `p_c=(1+exp(c))^-1`, `m(c)=1+c p_c`, and let `P=K_{c0}` where
`c0=mu-(sigma/2)sqrt(L/N)`.  Let `m0=m(c0)` and let `f_j` be the density
of a sum of `j` baseline increments.  Suppose

`q=N m0-sigma sqrt(N L)+O(L)`.

Use the signed exponential tilt

`dP_eta/dP=exp(eta Z)/M(eta)`

and choose `eta=eta_N` with `m_eta=q/N`.  Uniformly here,

`eta=-(1/sigma)sqrt(L/N)+O(L/N)`, `v_eta -> mu`.

Indeed the tilted increment is `Exp(1-eta)` plus a bounded Bernoulli
shift with success probability

`p_{c0,eta}=p_{c0}exp(eta c0)/(1-p_{c0}+p_{c0}exp(eta c0))`.

For `c0` near `mu` and `|eta|=o(1)`, its variance is uniformly positive,
all fixed moments are uniform, and its centred characteristic function is
bounded by the exponential factor

`(1+u^2/(1-eta)^2)^(-1/2)`.

Consequently the usual small-/large-Fourier-interval argument is uniform
in this triangular family:

`sqrt(j) f_j^eta(j m_eta+sqrt(j) w) = phi_{v_eta}(w)+O(j^-1/2)`

uniformly for bounded `w` and `j -> infinity`.  The same assertion at
`N-j` applies below.  This is the signed-tilt analogue of the r028 Fourier
lemma; no lattice qualification is needed because of the exponential
component.

## Exact conditional density and its limit

Let `j in [L^4,T]`, `T=N/(C1 L)`, and write

`r=j m0+z sqrt(j)`, with `z` in a fixed compact set.  Set

`d=r-jm_eta=sqrt(j)(z+s_j)`,
`s_j=sqrt(j)(m0-m_eta)`.

Here `|s_j|<=C/sqrt(C1)`.  Change of measure cancels exactly in the
conditional density:

`k_j(r) = f_j^eta(r) f_{N-j}^eta(q-r)/f_N^eta(q)`.

Since `q=N m_eta`, the second numerator is evaluated at its tilted mean
minus `d`.  The uniform local Fourier estimate gives

`f_{N-j}^eta(q-r)/f_N^eta(q)
 = sqrt(N/(N-j))*exp(-d^2/(2(N-j)v_eta))*(1+o(1))
 = 1+O(j/N+N^-1/2)`,

uniformly in this range (a weaker `o(1)` bound is already enough).  Thus
the conditional density of `z` is, uniformly on the source band,

`sqrt(j) k_j(jm0+z sqrt(j))
 = phi_{v_eta}(z+s_j)+O(j^-1/2+o(1)).`

This is a density statement; conditioning on the continuous value `S_N=q`
is therefore unambiguous via the regular conditional density.

## Cutoff second moment and the exact corrector identity

Let `psi=theta(2-theta)`, with a fixed smooth cutoff supported on
`c<=|z|<=C`.  First choose `c` small and `C` large.  Continuity of Gaussian
integrals then gives, for all sufficiently small `|s|` and `|v-mu|`,

`integral z^2 psi(z) phi_v(z+s) dz >= (1-eps) mu`.

Next choose `C1` large to make `|s_j|` uniformly small, and then take `N`
large.  Integrating the displayed conditional density over the conditional
law of `S_j` gives

`E[ ((S_j-jm0)^2/j) psi((S_j-jm0)/sqrt(j)) | S_N=q]
 >= (1-2eps) mu`.

For the source

`g_j(r)=alpha f_j(r) ((r-jm0)^2/j^2) psi((r-jm0)/sqrt(j))`,

and the exact recurrence `U_t=P U_{t-1}+g_t`, convolution and Bayes give
the exact identity

`U_N(q)/f_N(q)
 = alpha sum_{j=L^4}^T (1/j)
 E[ z_j^2 psi(z_j) | S_N=q]`.

Since

`sum_{j=L^4}^{N/(C1 L)}1/j = L-5 log L-log C1+o(1)`,

we obtain

`U_N(q) >= alpha (1-2eps) mu f_N(q) L (1-o(1))`.

With `alpha` chosen below `p/8` to leave the one-step and corrector-error
reserve, this is the desired bridge lower bound.  In particular its leading
candidate coefficient is `mu p/8=(mu-1)/8`.

## Audit qualifications

The notation in the penultimate displayed expectation is shorthand for the
conditional expectation over `S_j`; it is *not* conditioning simultaneously
on a prescribed value of `S_j`.  The required theorem is the uniform
signed-tilt local Fourier lemma above.  The ratio estimate need only be
`1+o(1)`; its stated sharper rate is not needed for the harmonic-sum result.
