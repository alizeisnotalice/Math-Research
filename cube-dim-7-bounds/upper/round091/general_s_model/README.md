# R091 general-s finite numerical model

Scope: this is a single preregistered floating-point batch for c=1.1 using the supplied coarse general-s constants. The R090 proof and constants are not frozen here; no rational certificate, all-s proof, C_n result or global optimization claim is made. The fixed tau=n exp(-delta L) cost uses

    1/n+exp(-delta L)+exp(-sigma L)+L exp(-(sigma-delta)L),
    tau/k <= rhi L² exp[-(delta-chi(1-sigma))L].

The registered 15 equal-length subintervals partition [5,64] up to endpoints. Each worker used differential evolution with at most 1722 evaluations, BLAS restricted to one thread, and no artificial memory allocation. The actual foreground process ended before the 1-second tool yield returned a live session: parent PID 32563 and workers 32572–32586 all completed with exit code 0. There were 22,554 evaluations; 13 intervals returned model candidates and 2 returned none in their bounded evaluations. Rejected evaluation categories were 2781 nonpositive full-cost surplus, 312 secondary-budget failures and 96 domain failures; these are algorithm counts, not exhaustive impossibility proofs. The two candidate-free runs ended early due optimizer convergence at constant penalty.

The parent sampled `ps` five times at 0.1-second intervals. The sampled aggregate worker CPU maximum was 986.4%, summed RSS peak 1106208 KiB, wall 0.764 s; individual worker CPU times are `process_time` end minus worker-start baseline, totaling 3.369 s. These are genuine samples, not claims about unsampled peaks; summed RSS can count shared pages repeatedly. No run was repeated or prolonged to increase utilization.

Top three float diagnostics (raw full-precision parameter fields are in `summary.json`):

| worker | s | chi | sigma | delta | qmax | gamma | modeled log2 L |
|---|---:|---:|---:|---:|---:|---:|---:|
| 9 | 42.57987 | 12.90459 | 0.92808239 | 0.92807606 | 0.99999352 | 3.7331913 | 32.19323 |
| 10 | 44.37905 | 13.48821 | 0.93097957 | 0.93096998 | 0.99999188 | 3.7743236 | 32.19566 |
| 8 | 40.36881 | 12.19160 | 0.92419529 | 0.92418602 | 0.99999259 | 3.6805187 | 32.19742 |

These narrow margins are sensitive to parameter rounding. The modeled thresholds use ordinary floating exponentials and a log-domain N16 surrogate, not S16 rational enclosures, and do not certify every geometric finite gate. Any future candidate must fix fresh rational parameters, register them, verify all guards including the exact finite-history geometry, and wait for the actual R090 proof/constant audit.

The model implements the supplied coarse s-dependent constants: |J-1|<=2s/rho; |d log J|<=(4s²+6s+1)/rho; J1<=s/rho, |dJ1|<=5s/rho and J2<=s/rho; baseline 4s²n/rho²+4s/rho and variance 16s²n/rho²; source pin 1 and variance <=2s²x. The implied log transport is bounded by 2(4s²+6s+1)xh and profile transport by 10s xh. It enforces e<1/[2(s+1)], source rho>4(s+1), and uses target-cap fraction 1-1/(4rhi), not the earlier fixed .95 bound. Actual availability of the stated bounds for all target and source faces remains an external hypothesis.

`registration.json` preceded execution; `run_batch.py`, worker receipts and `summary.json` preserve the single batch. Rerunning the script overwrites PID-bearing receipts, so preserve originals first. No old round, proof, paper or Git state was modified or submitted.
