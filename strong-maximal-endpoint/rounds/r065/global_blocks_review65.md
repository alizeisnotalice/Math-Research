# Final audit — global fixed-block halo

**PASS.**  `global_blocks65.md` validly globalizes the fixed-block result
to the full weighted block-model halo with receiver-dependent action
vectors in \(I_K^r\):
\[
 \nu\{u:\sup_{\mathbf c\in I_K^r}E_{\mathbf c}H(u+T)>1\}
 \le(1+o(1))Q_NJ(H).
\]

The angular grid cut loses only \(A N^{-K_0}\), its compensation costs
\(\beta^N=1+o(1)\), and random grid shifts give the exact bounded-overlap
budget factor \((1+2R/W)^{r-1}=1+o(1)\).  Cellwise normalization plus the
r026 point-budget bound rules out the negative radial tail.  The positive
tail is first compared after normalization to \(Q_Ne^{-q_0}\), then scaled
back by the cell budget; this is what makes summation over infinitely many
cells legitimate.  The remaining strip is covered by the universal-action
result and the unconditional fixed-action block transfer.

Dependencies: r026 deep block construction and point budget; r028
fixed-global-ratio all-input weak bound; r062--r063 envelope, location,
and vector-curvature estimates.  The result remains a fixed-\(r\), bounded
cap, equal-block, near-\(\mu\) **block-radial model** theorem.  It does not
settle action-dependent within-block angular transport, arbitrary physical
rectangles, or general \(C_N\).
