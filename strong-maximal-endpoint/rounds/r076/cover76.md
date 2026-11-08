# r076: local corner nets, global source reuse, and necessary entropy

## Conclusion

At one fixed receiver scale, no rho^{-N} action net is needed: one geometric
cell and one enlarged negative corner cover every vector action in the
critical window, at an exactly computable threshold loss. This is not a
global finite-cover theorem. Across scales, the same source region occurs in
infinitely many enlarged corner boxes, so an unweighted Fubini sum of source
budgets diverges even for one bounded compact packet. Moreover, the known
cap-two witness lower bound forces exponential anchor entropy for any corner
cover which purports to charge the full witness halo.

Let

\[
 c_-=c_r-K\rho,\quad c_+=c_r+K\rho,\quad
 r_+=e^{-c_-},\quad r_-=e^{-c_+},\quad
 \rho=\sqrt{\log N/N}.                                   \tag{1}
\]

## One-parameter geometric grid at a fixed scale

Fix a>1 and, for m in Z^N, set

\[
 C_m=\prod_{i=1}^N[a^{m_i},a^{m_i+1}),\qquad
 A_{m,i}=-r_+a^{m_i+1}.                                   \tag{2}
\]

If x is in C_m and c_i is in [c_-,c_+], then the actual lower face
-e^{-c_i}x_i is no smaller than A_{m,i}, so

\[
 R_{\mathbf c}(x)\subset\prod_i[A_{m,i},x_i].             \tag{3}
\]

The enlarged rectangle has one fixed negative corner. Its volume relative
to the actual rectangle is bounded coordinatewise by

\[
 {x_i+r_+a^{m_i+1}\over(1+e^{-c_i})x_i}
 \le C(a,\rho):={1+a r_+\over1+r_-}.                     \tag{4}
\]

Hence actual success of a nonnegative source implies that the enlarged
fixed-corner average is at least

\[
 C(a,\rho)^{-N}=e^{-\tau_N},\qquad
 \tau_N=N\log C(a,\rho).                                  \tag{5}
\]

This is only a geometric reduction: the near-corner theorem must be used at
the threshold e^{-tau_N}, not at level one. Since

\[
 \log C(a,\rho)=
 \log{1+a e^{-c_r}\over1+e^{-c_r}}+O(\rho),               \tag{6}
\]

any permitted linear threshold tau_N<beta N allows a fixed a>1 sufficiently
close to one and then all sufficiently large N. Thus the whole critical
vector-action window fits one corner locally; local action entropy is not
the obstruction.

For the sharp leading cheap threshold
\(\beta_*=\mu-1-\log\mu\), this condition can be written exactly as

\[
 1<a<a_{\max,N}:=
 {e^{\beta_*}(1+r_-)-1\over r_+}.                         \tag{6a}
\]

Indeed, it is precisely \(C(a,\rho)<e^{\beta_*}\).  As \(N\) tends to
infinity, \(r_+,r_-\) tend to \(\mu-1\), and

\[
 a_{\max,N}\longrightarrow
 {e^{\mu-1}-1\over\mu-1}>1.                               \tag{6b}
\]

For a finite N it is possible that \(a_{\max,N}\le1\); then this particular
leading-cheap grid is unavailable at that N.  The assertion here is
asymptotic and does not blur that finite-N qualification.

At one anchor, all 2^N orientations are free: partition the source into its
disjoint orthants relative to that anchor, reflect each to the positive
orthant, and sum the modulars. Only genuinely different anchor locations
count below.

## Why an unweighted sum over cells fails

Let

\[
 B_m=\prod_i[-r_+a^{m_i+1},\,a^{m_i+1}] .                 \tag{7}
\]

Localizing a fixed-corner charge for C_m naturally replaces a source by
F 1_{B_m}. The desired unweighted allocation inequality would be

\[
 \sum_{m\in\mathbb Z^N}\mathcal J(F\mathbf1_{B_m})
 \ \lesssim\ \mathcal J(F).                               \tag{8}
\]

It is false. Take any nonzero bounded packet F=t 1_Q, 0<t<=2, supported in
a small cube Q about the origin. For every m with all coordinates
sufficiently large, Q is a subset of B_m, and

\[
 \mathcal J(F\mathbf1_{B_m})=\mathcal J(F)>0.
\]

Already along m=(M,...,M), infinitely many terms are equal to the same
positive number. The left side of (8) diverges although the global source
budget is finite. The same ancestor reuse persists for a disjoint multiscale
family with summable total modular: every nonzero packet belongs to
infinitely many larger corner boxes.

This is not an example of one packet producing successful receivers at all
large scales; a compact packet cannot do that at a fixed average level. It
only proves that unweighted source reuse pays a packet before it has used
the scale decay of its averages.

For a completely explicit instance, take \(a=2\) and
\(F=\frac12\mathbf1_{[0,1]^N}\).  Then
\(\mathcal J(F)=\frac12\), and for every diagonal ancestor
\(m=(M,\ldots,M)\), \(M\ge0\), one has \([0,1]^N\subset B_m\).
The first K such ancestors already contribute \(K/2\) to the left side of
(8).  On the other hand any actual receiver in that scale cell has
\[
 H_{\mathbf c}F(x)\le\tfrac12\,2^{-MN}<1.
\]
Thus the repeated-budget calculation is explicitly not a repeated-success
example.

## Necessary entropy

The full-halo anchor count and independently assumed source-budget reuse bounds are proved in [entropy76.md](entropy76.md). For sublinear enlargement each is at least exp((mu−1)N−o(N)) when used to charge an entire efficient witness.

## Qualitative finite exact-corner cover

There is no topological obstruction to a finite cover for one fixed compact
strict halo. Let F be bounded and compactly supported, and let

\[
 E\Subset\{M_{\rm S}F>1\}
\]

be compact.  For every x in E choose a rectangle R containing x with
average strictly above one. By an arbitrarily small enlargement, strictness
allows x to be put in the interior of R while retaining that inequality.
Cut R by the coordinate hyperplanes through x into its \(2^N\) vertex
subrectangles. The average on R is the volume-weighted average of the
subrectangle averages, so at least one subrectangle has average above one.
It has x as a vertex and an opposite fixed corner a_x, with all side lengths
positive.

For the corresponding orientation, the fixed-corner average based at a_x
is continuous in the opposite vertex near x: the numerator is a bounded
integral over moving rectangles and the denominator stays positive. Hence
there is an open neighborhood U_x on which that same fixed-corner average
is above one. The U_x cover E, and compactness gives a finite subcover.
This is a finite exact-corner cover with \(\tau=0\).

For the known cap-two compact witness, take a strict compact subset of its
finite strict halo by inner regularity. Keeping, for example, a fixed
positive fraction of the halo volume loses only O(1) in its logarithmic
score, so it still has score
\[
 \exp\{q_0+\alpha\log N+o(\log N)\}.
\]
The conditional entropy converse above is therefore not vacuous: every
finite exact-corner cover of such a retained witness must obey
\[
 \log M_N\ge q_0+\alpha\log N-\log B_N(0)+o(\log N)
 =(\mu-1)N-d_N+\alpha\log N+o(\log N).                   \tag{12a}
\]

The set of cardinalities of finite exact-corner covers is nonempty, so an
integer minimal cardinality exists. This makes no claim of a quantitative
net, an algorithm for the minimizer, a source-budget allocation, or any
upper bound for general halos.

## Meaningful reduced target

A scale-aware cover would need cell weights w_m with

\[
 \sup_y\sum_{m:\,y\in B_m}w_m\le C                         \tag{13}
\]

and a compatible per-cell halo estimate charging
w_m J(F 1_{B_m}), or an equivalent tree stopping functional. Uniformly
positive weights violate (13) by the packet above. Decaying weights can
make ancestor sums finite, but their reciprocal loss must be recovered from
a real radial/scale decay in the successful averages. Establishing that
recovery is the remaining stopping/packing problem; no general upper theorem
is claimed here.
