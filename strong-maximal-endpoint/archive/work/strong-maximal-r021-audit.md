# r021 independent theory audit — 2026-11-06

Read-only audit of `rounds/r021/report.md`; no repository/Git edit.

## Verdict

Accepted. The nonlinear shifted Gaussian barrier supplies the stated
order-λogarithmic general upper bound; no hidden (L) loss or barrier-closure
gap was found.

## Key checks

* The exact chain-rule residual (5) is correct. In the quadratic core it
  reduces to (6): the leading part is
  `-K(1+z^2)f/t`, while the two nonlinear errors have relative size
  `O(K/sqrt(C0)+K^2/C0)`. The logarithmic-in-time part is controlled by
  `log(t/t0)/sqrt(t)`, so it brings no hidden `sqrt(L)` loss.
* The C3 cutoff is correctly matched: its derivative factorization and values
  at `v=0,1` agree with the quadratic core and constant cap (also verified in
  Wolfram). It supplies global residual bounds for the escaped region.
* In the spatial moment comparison, `W_xx/W_x` and `W_xxx/f` have the stated
  scales after `C0>>K^2`. Completion of the mean/variance term produces only
  `C(1+z^2)f/t`, with a coefficient independent of K. Thus the choice order
  `C_*`, then K, then C0 is noncircular.
* Duhamel has the correct sign. On `|U|<=1`, the shifted density and
  `1+z^2` are uniformly comparable without an L-dependent constant; the
  in-wide-band remainder is nonnegative because `LW<=0`. Only Gaussian escape
  from the wider band uses the global bound, and it is exponentially small in
  `tL`. This yields the needed positive discrete gain.
* The capped outer bands and signed initial function close exactly as in r020.
  In particular, the right band uses `H<=x/4` after selecting C0 and gets a
  normal plus exponential-tail error below delta.
* The CDF-to-excess integration is finite-range plus a positive-mgf tail.
  With mean shift `O(L)`, it gives `d+O(L)`. The shorter middle interval then
  gives the `O(L^(3/2)/N)` accounting error. The final normal shift gains
  `L^2/N`, which strictly dominates it.

No blockers found in Duhamel, derivative uniformity, the capped domain, or
the final probability comparison.
