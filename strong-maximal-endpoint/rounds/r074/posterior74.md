# Conditional posterior moderate deviations and action likelihood ratio

## Result and scope

Put \(L=\log N\), let \(\mu>1\) solve
\(\mu=1+e^{-\mu}\), and write
\[
 p=p_\mu={1\over1+e^\mu}={\mu-1\over\mu},\qquad
 \tau^2={p\over\mu^2}.
\tag{1}
\]
For fixed \(C\), uniformly for
\[
 |s-N\mu|\le C\sqrt N,\qquad
 c\in[\mu-\varepsilon,\mu+\varepsilon],
\tag{2}
\]
with \(\varepsilon>0\) sufficiently small, consider the posterior
probability mass function
\[
 \pi_{c,s}(k)={w_k(c,s)\over\sum_jw_j(c,s)},\qquad
 w_k(c,s)={N\choose k}(s-ck)_+^{N-1}.
\tag{3}
\]
The following is the needed conditional result.  It is a fixed-central-slab
statement; it makes no assertion uniform up to a radial endpoint.

Let
\[
 c_r=\mu-{\sqrt\mu\over2}\sqrt{L/N},\qquad
 c_a=c_r+{L^{1/4}\over\sqrt N},\qquad
 \delta_N={1\over NL^2}.
\tag{4}
\]
For every \(s\) in (2), there is an upper-tail event \(A_s\) (with a
fraction of one lattice atom allowed) such that
\[
 \pi_{c_a,s}(A_s)=\delta_N,
\tag{5}
\]
and
\[
 \pi_{c_r,s}(A_s)
 \ge \delta_N
 \exp\!\left\{\sqrt{2p}\,L^{3/4}-C_C\sqrt L\right\}.
\tag{6}
\]
The constant in the exponential is positive.  In particular, (6) is the
posterior comparison required by the deletion construction in r073.  The
proof below is by the exact finite sum (3), Stirling's formula, and Taylor
estimates; no general triangular-array moderate-deviation theorem is being
invoked.

The same calculation gives the smooth posterior-score estimate needed for
the radial comparison.  If \(\sigma=s/N\), let \(x=x(c,\sigma)\) be the
unique solution of
\[
 \log{1-x\over x}={c\over\sigma-cx},
\qquad \gamma={1\over\sigma-cx}.
\tag{7}
\]
Thus \(x=p_{\gamma c}\).  Uniformly on any fixed compact subregion of the
interior parameter set in (2),
\[
 {\mathbb E}_{c,s}{K\over s-cK}=\gamma p_{\gamma c}+O(N^{-1}),
\qquad
 (N-1){\mathbb E}_{c,s}{K\over s-cK}
   =N\gamma p_{\gamma c}+O(1).
\tag{8}
\]
This is an expectation statement for the full labelled posterior after the
label count has been summed into (3); it does not replace a common-source
or anchor argument.

## 1. A uniform discrete Laplace lemma

Set
\[
 \ell_{c,\sigma}(x)
 =-x\log x-(1-x)\log(1-x)+\log(\sigma-cx).
\tag{9}
\]
Shrinking \(\varepsilon\) if needed, the solution in (7) stays in a
closed interval strictly inside both \((0,1)\) and \((0,\sigma/c)\).  It
is the unique maximizer of (9), because
\[
 \ell'(x)=\log{1-x\over x}-{c\over\sigma-cx},\qquad
 \ell''(x)=-{1\over x(1-x)}-{c^2\over(\sigma-cx)^2}<0.
\tag{10}
\]
All derivatives through order four are therefore uniformly bounded on a
fixed neighborhood of the maximizers, and \(-\ell''\) is uniformly bounded
away from zero.

For completeness, the elementary discrete Laplace calculation is recorded
here.  Stirling's formula with its usual bounded remainder gives, uniformly
on that neighborhood,
\[
 \log w_k=C_N+N\ell_{c,\sigma}(k/N)+a_{c,\sigma}(k/N)+O(N^{-1}),
\tag{11}
\]
where \(C_N\) is independent of \(k\), and
\[
 a_{c,\sigma}(x)=-\tfrac12\log[x(1-x)]-\log(\sigma-cx).
\tag{12}
\]
The harmless constant \(-\frac12\log(2\pi N)\) has been absorbed into
\(C_N\).  Taylor expansion at the maximizer, followed by comparison of a
lattice Gaussian sum with its integral, yields
\[
 \pi_{c,s}(k)
 ={e^{-z^2/2}\over\tau(c,\sigma)\sqrt{2\pi N}}
 \left(1+O\!\left({1+|z|^3\over\sqrt N}\right)\right),
\quad
 z={k-Nx(c,\sigma)\over\tau(c,\sigma)\sqrt N},
\tag{13}
\]
whenever \(|z|\le C_1\sqrt L\).  Here
\[
 \tau(c,\sigma)^2={1\over-\ell''_{c,\sigma}(x(c,\sigma))}.
\tag{14}
\]

Here are details that make the uniformity in (13) explicit.  On
\(|k/N-x|\le N^{-2/5}\), Taylor's remainder in \(N\ell\) is bounded by
\(C|z|^3/\sqrt N\).  The first-order variation of \(a\) is
\(O(|z|/\sqrt N)\).  On the complement of this neighborhood strict
concavity gives a contribution at most \(e^{-c_0N^{1/5}}\) relative to the
central sum.  Thus the normalizing sum is the Gaussian Riemann sum with a
uniform relative error \(O(N^{-1/2})\), and the displayed point-mass
formula follows.  These estimates use only constants from the compact
parameter set.

The same argument, with a one-sided Gaussian sum, gives the logarithmic
tail estimate
\[
 \log \pi_{c,s}\{K\ge Nx+\tau\sqrt N z\}
 =-{z^2\over2}-\log z+O_C\!\left(1+{z^3\over\sqrt N}\right)
\tag{15}
\]
for \(1\le z\le C_1\sqrt L\), with the usual rounding of the threshold.
Indeed, apply (13) up to \(z+1\), use the successive-ratio form of (11)
to bound the remaining terms by a geometric integral, and use
\[
 {e^{-z^2/2}\over\sqrt{2\pi}(z+z^{-1})}
 \le \int_z^\infty {e^{-u^2/2}\over\sqrt{2\pi}}\,du
 \le {e^{-z^2/2}\over\sqrt{2\pi}z}.
\tag{16}
\]
For terms beyond the Taylor window, (10) makes the ratio decrease by a
fixed amount per macroscopic displacement, so the same upper bound closes.
This proves (15), rather than merely quoting a normal-tail theorem.

Finally, (11) also proves a uniform smooth-test rule.  For a function
\(q\) with two bounded derivatives near the compact set of saddles,
symmetry of the leading Gaussian and the cubic correction in (11) give
\[
 {\mathbb E}_{c,s}q(K/N)=q(x(c,\sigma))+O(N^{-1}).
\tag{17}
\]
More explicitly, the odd \(N^{-1/2}\)-term in the local density can only
produce \(O(N^{-1})\) after multiplying the linear Taylor term of \(q\);
the quadratic term has expectation \(O(N^{-1})\).  The rational function
\(q(u)=u/(\sigma-cu)\) needs a separate endpoint argument: it is not
legitimate to call it polynomially bounded when \(s-cK\) can be arbitrarily
small.  Instead use its exact unnormalised numerator
\[
 w_k(c,s){k\over s-ck}
 ={N\choose k}\,k\,(s-ck)_+^{N-2}.\tag{17a}
\]
For \(N\ge3\), the same Stirling--Laplace comparison applies to (17a).
On the complement of a fixed saddle neighborhood, its exponent has a
strictly smaller compact maximum than the interior one (and tends to
\(-\infty\) at \(s-ck=0\)); the extra factor \(k\) and the change
\(N-1\) to \(N-2\) affect only subexponential terms.  Hence that whole
endpoint contribution is exponentially small relative to the normalizer.
The smooth cutoff argument is then valid on the remaining compact
neighborhood.  Since
\[
 {K\over s-cK}={K/N\over\sigma-cK/N},\qquad
 {x\over\sigma-cx}=\gamma p_{\gamma c},
\tag{18}
\]
(17) gives (8).

## 2. Location of the actual \(\delta_N\)-tail

At \((c,\sigma)=(\mu,\mu)\), (7) gives \(x=p\) and
\(\sigma-cx=1\).  Consequently
\[
 -\ell''_{\mu,\mu}(p)
 ={1\over p(1-p)}+\mu^2={\mu^2\over p},
\qquad \tau(\mu,\mu)^2={p\over\mu^2}.
\tag{19}
\]
Apply (15) with \((c,s)=(c_a,s)\).  The upper-tail quantile having mass
\(\delta_N\), with a fractional boundary mass if necessary, has the form
\[
 k_s=Nx(c_a,s/N)+\tau_a\sqrt N\,
       \sqrt{2L+O_C(\log L)}+O(1).
\tag{20}
\]
Here \(\tau_a=\tau(c_a,s/N)=\tau+O_C(\sqrt{L/N})\).  Hence using the
base value \(\tau\) instead would give the same displayed precision after
altering the \(O_C(\log L)\) term.
To see this without a circular normal approximation, insert the two trial
values \(z^2=2L\pm C_2\log L\) into (15) and choose \(C_2\) large.  Their
tail probabilities respectively straddle \(N^{-1}L^{-2}\).  Monotonicity
of the tail then proves (20).  This also shows that the cutoff lies a fixed
positive fraction of \(N\) away from both posterior endpoints.

## 3. Exact monotonicity and the likelihood-ratio gain

The change of saddle under the action displacement in (4) has the required
sign and size.  Differentiating (7) at \((\mu,\mu)\) gives
\[
 \partial_cx\big|_{(\mu,\mu)}=-{p\over\mu}.
\tag{21}
\]
For reference, if \(F(c,x)=\log((1-x)/x)-c/(\sigma-cx)\), then at the
base point \(F_c=-\mu\) and \(F_x=-\mu^2/p\), which proves (21).
Thus, uniformly for the central \(s\)-slab,
\[
 x(c_r,s/N)-x(c_a,s/N)
 ={p\over\mu}{L^{1/4}\over\sqrt N}
 +O_C\!\left({L^{3/4}\over N}\right).
\tag{22}
\]
The error in (22) is more than sufficient below.  Combined with (19), it
says that the two posterior centers are separated by
\[
 {N[x(c_r,s/N)-x(c_a,s/N)]\over\tau_a\sqrt N}
 =\sqrt p\,L^{1/4}+o(1)
\tag{23}
\]
actual posterior standard deviations.  It also recovers
\(\partial_c\mathbb E(K\mid S=s)=-Np/\mu+o(N)\) on this slab.

There is an exact monotonicity which avoids any issue caused by the
fractional boundary in (5).  For every integer in the common support,
\[
 {d\over dk}\log{w_k(c_r,s)\over w_k(c_a,s)}
 =(N-1)\left\{-{c_r\over s-c_rk}+{c_a\over s-c_ak}\right\}>0,
\tag{24}
\]
because after clearing denominators the numerator is
\(s(c_a-c_r)>0\).  Division by the two normalizing constants preserves
this monotonicity.

For clarity, the size of the ratio at the cutoff can also be obtained from
the same Taylor calculation rather than from a formal LAN citation.  Put
\(d=x(c_r,s/N)-x(c_a,s/N)\), and write
\(k=Nx(c_a,s/N)+\tau_a\sqrt N z\), with \(z=O(\sqrt L)\).  Subtracting the
two versions of (11), including their Laplace normalizers, gives
\[
 \log {\pi_{c_r,s}(k)\over\pi_{c_a,s}(k)}
 ={d\over\tau_a^2}\,\tau_a\sqrt N z
  -{Nd^2\over2\tau_a^2}+O_C(1)
 =\sqrt p\,L^{1/4}z+O_C(\sqrt L).
\tag{25}
\]
The first equality has a uniform \(O(1)\) Taylor remainder: cubic local
terms are \(O(L^{3/2}/\sqrt N)\), while the change in curvature across
\(c_a-c_r\) is negligible at this scale.  The second, negative quadratic
term is only \(O(\sqrt L)\).  Formula (22) absorbs the remaining compact
parameter variation.

At the boundary point of the event in (5), (20) and (25) imply
\[
 \log {\pi_{c_r,s}(k_s)\over\pi_{c_a,s}(k_s)}
 \ge \sqrt{2p}\,L^{3/4}-C_C\sqrt L.
\tag{26}
\]
Define \(A_s\) as the strict upper tail above the first boundary atom and
include the required fraction of that atom to make its actual mass exactly
\(\delta_N\).  Every included atom has likelihood ratio at least the one
in (26), by (24).  Multiplication by the actual mass proves (6).

## Limitation

The result is conditional on the radial sum and compares two fixed actions.
It supplies the local posterior-tail loss and the score bias used in the
reference-baseline calculation.  It does not make the deleted sets for
different receivers arise from one common physical source, and hence does
not settle the separate common-source, anchor, or packing interfaces.
