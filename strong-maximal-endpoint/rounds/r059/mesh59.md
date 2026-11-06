# Geometric mesh: conditional rank, actual cost, and outer escape

This note concerns one explicit mesh family only.  It does not assert a general endpoint theorem.

## Construction and conditional rank

Fix `N,M >= 2`, let
\[
 \Phi_N(t)=t(1+\log^+t)^{N-1},\quad Q=1+NMR,\quad
 r_0=0,\quad r_j=Q^{j-1}\ (1\leq j\leq M),\quad h_j=r_j-r_{j-1}.
\]
On outer slab `i`, start with the inner fibre \(1_{[0,r_M]}\), lower
cell \([r_{i-1},r_i]\) by \(D_i/h_i\), where \(D_i=MRr_i\), and, if
`i<M`, raise the next cell by \(D_i/h_{i+1}\).  The levels are
\[
 1-MR\ (i=1),\qquad 1-MR-1/N\ (i\geq2),\qquad 1+1/N.
\]
They are nonnegative if \(MR+1/N<1\).

At prefix \([0,r_j]\), fibre `j` averages \(1-MR\), and all other
fibres average one.  Thus all `M` fibres average \(1-R\), while omitting
any fibre `j` lets the same prefix make every retained fibre equal one.
Between grid radii, the total loss is still exactly \(MRs\) at prefix
length `s`: on each cell the current downshift minus preceding recovery
has rate \(MR\).  The family therefore has rank `M` and conditional loss
`R` at every grid prefix.

## Exact actual modular cost

After division by outer length \(M\) and inner length \(r_M\), direct
cell integration gives
\[
 \boxed{\mathcal B_N=1-R+\frac{1-1/r_M}{M}
 \left[\Phi_N(1+1/N)-(1+1/N)\right].} \tag{1}
\]
All lower values are at most one, hence contribute exactly the total
linear decrease \(-R\); the recovery cells, of total inner length
\(r_M-1\), give the displayed nonlinear term.  As
\(\Phi_N(1+1/N)\leq(1+1/N)^N\leq e\),
\[
 \mathcal B_N\leq1+(e-1)/M.
\]
Consequently \(M=\lfloor N^{1/4}\rfloor\),
\(R=\lceil(\log N)^2\rceil/N\) has `M -> infinity`, critical-sized
`R`, and \(\mathcal B_N\to1\); all amplitudes also tend to one.
Diagonal rescaling sends the support to a unit cube without changing
normalized averages or normalized modular cost.  The large quantity is
the internal mesh span \((M-1)\log Q\), not support aspect ratio.

## Global outer escape: the exact correction

Rescale outer support to \([0,1]\), put the receiver at \(1+\eta\), and
keep the inner receiver coordinate zero.  Write \(w=(M-1)/M\).  Full
outer averaging has value
\[
 S_{\rm full}(\eta)=\frac{1-R}{1+\eta}.
\]
A proper suffix retaining `s` slabs can omit a deficient fibre and use a
common maximizing inner prefix; its value is
\[
 \frac{s/M}{s/M+\eta},\qquad 1\leq s<M.
\]
Therefore the exact global value is
\[
 \boxed{S^*(\eta)=\max\left\{\frac{1-R}{1+\eta},\frac{w}{w+\eta}\right\}.} \tag{2}
\]
There is no hidden fractional-outer escape: with an inner prefix fixed,
the average is fractional-linear as the outer left endpoint traverses a
cell, so its maximum is at an endpoint.  The identical anchored-prefix
argument in the inner coordinate reduces the remaining candidates to the
grid; then the list above is exhaustive.  In particular \(\eta=0\)
has a proper suffix of value one, so conditional rank alone gives no
global gap at the boundary.

For adaptive benchmark \(A^*(\eta)=1/(1+\eta)\), this yields
\[
 A^*(\eta)-S^*(\eta)=
 \min\left\{\frac{R}{1+\eta},
 \frac{\eta}{M(1+\eta)(w+\eta)}\right\}. \tag{3}
\]
The two branches meet exactly at
\[
 \eta_* =\frac{(M-1)R}{1-MR}. \tag{4}
\]
To retain a fixed positive fraction of the full loss `R`, (3) requires
\(\eta\asymp MR\) when \(MR=o(1)\).  At the crossing,
\(\eta_*\sim(M-1)R\); before outer rescaling to unit length the
physical separation is \(M\eta_*\asymp M^2R\).

## Closing escape forces large actual cost

At (4), scale the input by
\[
 b_*=1+\eta_*=\frac{1-R}{1-MR}.
\]
Then \(b_*A^*=1\) and \(b_*S^*=1-R\) (a further \(1+o(1)\) factor makes
the inequality strict).  The recovery cell of fibre `M-1` in the largest
inner cell has normalized measure \((1-1/Q)/M\) and original height
\(1+1/N\).  If `Q>=2` (as in the critical-scale sequence below), monotonicity gives
\[
 \mathcal B_N(b_*f)\geq\frac{1-1/Q}{M}\Phi_N(b_*)
 \geq\frac1{2M}\Phi_N(1+\eta_*). \tag{5}
\]
For \(0\leq u\leq1/2\),
\(\log\Phi_N(1+u)\geq\frac49(N-1)u\).  If \(MR\leq1/4\), then
\(\eta_*\leq1/3\), and so
\[
 \log\mathcal B_N(b_*f)\geq
 \frac49(N-1)(M-1)R-\log(2M). \tag{6}
\]
For \(R\geq\kappa(\log N)^2/N\), a weaker uniform consequence is
\[
 \log\mathcal B_N(b_*f)\geq
 \frac{2\kappa}{9}(M-1)(\log N)^2-\log(2M).
\]
So this is a low-cost, growing-rank *conditional* example, but not a
cheap globally near-threshold example: a suffix escapes, or closing it
has exponential-in-\(NR\) modular cost.

## Wolfram verification

`mesh59.wl` directly integrates cells for `N=2,4,8,16`, `M=2,3,5`, and
\(R=1/(4M)\).  All 12 exact rows verify (1), grid-prefix identities, and
proper-subset witnesses.  It also stores rational diagnostics for
`N=64,256,1024,4096,16384,65536`, with
\(M=\lfloor N^{1/4}\rfloor\), \(L=\lceil(\log N)^2\rceil\),
\(R=L/N\), and integer \(Q=1+ML\), in `mesh59.json`.  At (4), all six
rows verify equality of the two branches in (2).  Raw costs are
approximately \(1.5271,1.3049,1.2951,1.1976,1.1504,1.1055\); the logs of
the one-recovery-cell scaled lower bound are
\(25.17,108.13,205.74,497.72,953.69,1845.94\), respectively.
