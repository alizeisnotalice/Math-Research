# Necessary anchor entropy for a full near-corner cover

## Conditional converse

This is a necessary counting consequence of the cap-two near-corner Hardy
charge.  It does not supply a quantitative optimal net. The compact strict receiver case has a qualitative finite exact cover by the separate compactness argument in [cover76.md](cover76.md); arbitrary unbounded halos are not asserted to have a finite cover.

Let \(F_N\) be a compact cap-two r019/r069 witness with modular
\(J_F=\int\Phi_N(F_N)\) and witness set \(E_N\) satisfying
\[
 {|E_N|\over J_F}\ge
 \exp\{q_0+\alpha\log N+o(\log N)\}.
\tag{1}
\]
Suppose, as an additional geometric hypothesis, that all of \(E_N\) is
covered by \(M_N\) post-sharing fixed-corner families, each of which has
actual-to-corner volume enlargement at most \(e^\tau\).  (All signed
orientations at one anchor are free, as in `nearcorner76.md`; \(M_N\)
counts only anchors/source partitions which cannot be shared.)  The
near-corner upper bound gives
\[
 |E_N|\le M_NB_N(\tau)J_F.
\tag{2}
\]
Combining (1)--(2) yields the necessary entropy inequality
\[
 \boxed{\ \log M_N\ge q_0+\alpha\log N-
                  \log B_N(\tau)+o(\log N).\ }
\tag{3}
\]

## Sublinear enlargement

For \(\tau=o(N)\), the sharp Hardy optimization from `nearcorner76.md`
gives
\[
 \log B_N(\tau)=N+\sqrt{2N\tau}+O(\tau+1).
\tag{4}
\]
Using \(q_0=N\mu-d_N\) in (3),
\[
 \log M_N
 \ge (\mu-1)N-d_N-\sqrt{2N\tau}-O(\tau+1)
       +\alpha\log N+o(\log N).
\tag{5}
\]
In particular,
\[
 {\log M_N\over N}\ge\mu-1+o(1).                        \tag{6}
\]

Thus a complete cover of a near-extremal halo by sublinearly enlarged,
cheap fixed-corner charges cannot have subexponentially many genuinely
distinct anchors/source partitions.  A finite-corner Hardy module can still
discard a subset of receivers, or a proof may carry the necessary
exponential anchor entropy by another budget mechanism; (5) does not rule
out either possibility.

The conclusion is only an upper-cover converse.  It neither constructs an
actual optimal anchor net nor claims a lower bound for all possible packing,
stopping, or action-indexed source-allocation arguments.

## Necessary source-budget reuse

There is a parallel converse if the cover is paid for by a countable family of separately assigned
source pieces rather than by an explicit anchor count.  Suppose
\(E_N\subset\bigcup_jE_{N,j}\), and for every \(j\) there is a cap-two
source \(F_{N,j}\) and one fixed anchor \(a_{N,j}\), independent of the receiver in \(E_{N,j}\), such that every assigned receiver succeeds through a near-corner rectangle based at that anchor, with enlargement at most \(e^\tau\). Signed orientations at this same anchor may vary.  Assume the proposed source
allocation has reuse factor \(A_N\):
\[
 \sum_jJ(F_{N,j})\le A_NJ(F_N).
\tag{7}
\]
For example, pointwise subpieces \(F_{N,j}\le F_N\) would be an additional structural restriction. It does not by itself imply (7); the summed reuse inequality is a separate hypothesis. No such allocation is constructed here.

Applying the one-corner charge to each \(F_{N,j}\), summing, and using
(1), gives
\[
 |E_N|\le B_N(\tau)\sum_jJ(F_{N,j})
 \le B_N(\tau)A_NJ(F_N),
\]
and hence
\[
 \boxed{\ \log A_N\ge q_0+\alpha\log N-
                  \log B_N(\tau)+o(\log N).\ }
\tag{8}
\]
For \(\tau=o(N)\), its expansion is the same as (5), with \(\log M_N\)
replaced by \(\log A_N\).  Thus a complete cheap near-corner charge of the
whole witness halo cannot simultaneously have subexponential distinct-anchor
entropy and subexponential source-budget reuse.  This is only a necessary
condition for an assumed allocation/cover; it does not prove that an actual
policy admits one.
