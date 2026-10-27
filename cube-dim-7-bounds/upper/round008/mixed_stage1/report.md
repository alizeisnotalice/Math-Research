# R008 signed mixed moments — Stage I stopping result

Stage I completed, but its first acceptance gate failed: the rigorously enclosed surrogate is

    0.90127456571435855668015 < L_H < 0.90127456571435855676623 < 1.

Thus the proposed lower surrogate does not prove a positive μ=Eν(H−1) lower bound. The true H is at least 1 by determinant Minkowski; the result above exposes slack in this particular lower surrogate, not a contradiction and not evidence that the actual matrix gain is absent. **Stage II was not started.** No new strict I, g, or Cn bound follows.

The four complete mixed-expectation enclosures have approximate values I₀≈0.03472141017593843, I₁≈0.03529716982053700, I₂≈0.03568713219857787, I₃≈0.03588397417371473. With exact positive binomial weights, these give L_H=8[2(I₀⁷+7I₁⁷+21I₂⁷+35I₃⁷)]^(1/7)/q. The complete R004 q(1/8) interval was copied with its manifest and reused; no old integral or shape sample was rerun. All dyadic enclosures are in certificate.json.

The requested formal calculation M_pair,lower/√L_H,upper is approximately 0.0853623313, above the declared 31/500 target. However, because L_H,upper<1, this substitution corresponds to a negative candidate μ, outside the nonnegative-moment domain of the Radau correction. It is recorded as a formal expression, not asserted to be a valid optimistic bound over feasible Radau moments or a lower bound on actual M. The already known H≥1 supplies only the zero-correction baseline from this surrogate; it does not certify improvement over the pair budget. The negative μ gate alone is decisive and requires this batch to stop.

All 256 signed mixed Mellin panels finished: 64 per j=0,1,2,3, including the first four per j as pilots. The implementation used Γ(z)F_j(z) times the coefficientwise conjugate of Γ(z)F_(7−j)(z). Every one of the 65 integrated even signed-real Taylor terms per panel and each panel's signed result is retained; no absolute-value integrand, panel clipping, or imaginary-zero assertion was used. Each group had 64 constant imaginary intervals excluding zero, and hundreds of strictly negative Taylor terms. The real mixed expectations, though positive, were assembled with the maximum of the two crossed Cauchy endpoint products and **both signs** of the full absolute tail. Cauchy raw errors are about 2.76–3.53e-28; normalised absolute tails about 8.29–8.31e-25. The 512-bit intervals and retained signed even terms support independent reconstruction.

The run took 6.305 s wall and 27.974 s summed panel CPU, about 4.44 core-equivalents over wall; one worker-phase snapshot observed 15 processes at 664.8% summed CPU and 422608 KiB summed RSS, which may include shared pages. Maximum individual RSS was 38,567,936 bytes and system swap use sampled at 3006.06 MiB. These are observations on a shared machine, not sustained occupancy claims.

The portable checker verifies 256 unique terminal panels, all signed terms, nonzero-imaginary examples, reused q hash, crossed radius envelopes and two-sided tails, and independently reconstructs the strict L_H<1 gate. It does not recompute the 256 series; run.py does. The separately hashed theory establishes the expectation identity and Cauchy bound. No Stage II, Monte Carlo, MATLAB, commit or push was performed.
