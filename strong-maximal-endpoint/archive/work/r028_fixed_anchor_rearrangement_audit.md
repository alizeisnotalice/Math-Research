# r028 audit: fixed-ratio anchored rearrangement

Audit date: 2026-11-06.  The proof in
`fixed_anchor_rearrangement_proof.md` is valid with the stated
coordinate-even/symmetrized-operator scope.

For a single positive-axis fibre,
\[
 A_rf(x)=\frac1{(1+r)x}\left(\int_0^xf(t)dt+\int_0^{rx}f(t)dt\right)
 =\int_0^1k_r(v)f(vx)dv,
\]
and the kernel in the physical fibre variable is decreasing.  Fibrewise
Hardy--Littlewood therefore has the correct direction
\(A_rf\le A_rf^\#\).  At a strict output level, the fibre halo is an initial
interval.  Integrating the interior-point inequality produces the global
two-level weight
\[
 W=(1_{t<a_z}+1_{t<ra_z})/(1+r).
\]
Its global decreasing rearrangement is exactly
\(k_r(y/D)1_{y<D}\): it has value \(2/(1+r)\) on measure \(rD\) and value
\(1/(1+r)\) on the remainder to measure \(D\).  Global
Hardy--Littlewood thus proves the claimed distribution inequality.  The
threshold \(s'>s\), then \(s'\downarrow s\), is the correct repair for
plateaux.  The endpoint cases \(r=0,1\) reduce respectively to a one-level
weight and the same one-level weight, so introduce no exception.

Iterating the one-axis distribution bound and the order-preserving scalar
operators gives the fixed-vector tensor bound.  For bounded compact input,
the induction is legitimate: a used coordinate has a \(\min(1,C/x)\) tail,
and an unused coordinate retains compact fibres, providing the required
fibre \(L^1\), boundedness and finite level measures.

The conversion to \(Q(n,{\bf r},A)\) is exact: for an interior rearranged
halo point \(x\), \(h(t)=f^*(xe^{-t})\) satisfies
\[
 J(h)=x^{-1}\int_0^x\Phi_n(f^*(y))dy,
\]
and the fixed-vector product average is \(\mathbb Eh(S_{\bf r})>1\).
The deep hyperbolic even input supplies the reverse inequality after a
strict feasible truncation.  Thus the cap and uncapped equalities follow.

Scope warning: the physical tensor identity uses coordinatewise-even inputs
and reflected intervals.  The result is not an upper bound for arbitrary
non-even inputs, adaptive ratio vectors, or the full strong maximal operator.
Within the coordinate-even fixed-r anchored family it is an exact finite-n
barrier, hence rules out all same-fixed-r microblock and multi-receiver
improvements, not merely the single-corner candidate.
