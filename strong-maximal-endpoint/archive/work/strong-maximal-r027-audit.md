# r027 common-c block log-concavity audit — 2026-11-06

Audit based on the supplied r027 statement and its exact external certificate;
formal r027 report was not yet on disk. No repository/Git change.

## Verdict

Accepted within its stated common-c, two equal-block, cap-2 scope.

* The Fraction script was actually rerun. It certifies `P'^2-P P''>=0` on all
  pieces and tails for b=4,5,6,7; the b=3 right-piece defect at one is exactly
  `-4`. Knot continuity supplies global log-concavity. Convolution closure
  then proves it for every common-c block order b>=4. Endpoints c=0 and
  c=infinity reduce directly to Gamma.
* `R_b=e^t f_b` is log-concave and nondecreasing. Hence all superlevels of
  the capped dual profile, which are product-R superlevels, are convex and
  exchange-symmetric.
* For each level set, log-concavity of the two independent shifted block
  densities times the convex-set indicator survives marginalization. Thus
  `B_a(d,v)` is even and log-concave in v, hence nonincreasing in `|v|`.
  Layer-cake preserves only this monotonicity; it does not assert that B is
  log-concave. The acceptance boundary therefore moves outward with `|v|`.
* In `(u1,u2)=(d+v,d-v)`, the Jacobian 2 cancels the radial integral, giving
  the exact receiver volume `integral exp(-2 d(v)) dv`. The support/height
  Markov argument gives the stated linear lower tails for d(v), hence total
  omitted tail at most `exp(-2R)`.

No claim of a real e^{-1/2} asymptotic is justified by these facts alone, and
none is made. The proof must retain common c; arbitrary blockwise c would not
follow from the certified convolution closure.

## Formal report final audit

Accepted. The exact report preserves all essential scope restrictions. The
portable Fraction `--verify` and the Wolfram all-piece check were actually
rerun successfully. The selected-kernel halo reduction is valid: the Markov
lower bound on each u-coordinate prevents negative-end mass escape, while
the exponential measure controls the positive end; compact exhaustion and
dominated convergence then yield (6). Levelwise Prékopa is used correctly
and the layer-cake conclusion is only evenness/imbalance monotonicity. The
boundary-volume identity and explicit `e^{-2R}` tail bound follow with the
correct Jacobian. §5 is clearly noncertified, and §7 labels the e^{-1/2}
calculation only as a pending conditional-limit target.
