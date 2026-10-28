# R078 fixed-tuple Bernstein diagnostic

Scope: numerical and exact arithmetic diagnosis conditional on the R073 full-chain hypotheses and the source-face independence and boundedness interface. No new C_n theorem is certified here. This is general-dimensional and does not concern n=7.

The fixed R073 geometry and all tolerances remain unchanged. On each captured signed source face, suppose the n−1 free-coordinate terms phi_i=-log ell_R are independent and lie in [0, log(6/5)] subset [0, 183/1000]; the pinned term is deterministic. The centered summands are then bounded above by 183/1000, while their total variance is at most V=(531/100000)x. With mean at most p+mx, p=183/1000, m=2211/[50000(1-2^-20)], and threshold chi x, the gap is u=(chi-m)x-p. These facewise statements must be proved for the actual inverse-capture product law before applying Bernstein, and the per-face results mixed only afterwards; no law is asserted for the actual increment measure. R073 states the needed product law and relative guard, but independent verification is outside this arithmetic run.

Bernstein gives failure at most exp[-u²/{2(V+bu/3)}], b=183/1000. At x=7392, u=36.7674482681353, V=39.25152, exponent=16.2895160690. The positive degree-64 exponential series, evaluated exactly as Fractions, gives exp(-E)<=1/S64(E)<8.424674e-8<1e-6. For x>=7392, u/x=chi-m-p/x is positive and increasing, and x (u/x)^2/[2(v+b(u/x)/3)] is increasing. Hence the same uniform per-face bound follows conditional on the interface. Bennett's numerical display 8.23015e-8 uses Decimal log/exp and is not an interval certificate; no Bennett improvement is used.

All reported deficit figures include the R073 target mass lower bound, Ulo/Uhi, radial correction, source and tau payments, and (1-eD), while retaining the same price bound 257/2500. They are conditional on unchanged R073 guards and chain:

| retention choice | (1-eD) Gamma | after price upper |
|---|---:|---:|
| R073 9717/10000 | 0.29108217809647 | 0.18828217809647 |
| Cantelli at base (unrounded) | 0.29110726933885 | 0.18830726933885 |
| Bernstein 999999/1000000 | 0.29955943113786 | 0.19675943113786 |

A separately registered exact Fraction guard verifies gamma=599/2000 below the latter budget by 0.0000594311 and gamma-257/2500-393/2000=1/5000. R073's existing scalar remainder was strictly less than 1/5000 at this same geometry; it was not rerun. This specifies an improved *conditional candidate* c=393/2000, not an accepted bound.

Reproduce with `python3 diagnostic.py` and `python3 candidate_guard.py`. Registration and amendment preceded successful execution. Initial attempted mpmath import failed before any calculations or output; Python standard-library Fraction/Decimal was used instead. The sole successful batch ran PID 25707, the candidate guard PID 25735. Full source input hashes are in registration.json. Git clone checkout failed on missing source Git objects, so no usable clone was established; this directory is within the incomplete clone, and the locked R070/R072/R073 working-tree input files were read from the original without modification. No old scans or 53-gate verifier were rerun, no MATLAB or GPU used, no commit/push performed.

## Final registered candidate

The subsequent fixed candidate gamma=599/2000, c=983/5000 has exact scalar allowance 1/10000. The separate `final_candidate_guard.py` recomputes the exact R073 positive-prefix scalar remainder at L=8192 and verifies it is strictly below this smaller allowance, together with the complete deficit gate. It uses the same previously proved tail monotonicity, conditional on the unchanged chain. The coordinator selected the earlier 393/2000 candidate; this extra guard remains as a separately recorded diagnostic.
