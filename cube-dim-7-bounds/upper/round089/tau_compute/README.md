# R089: independent finite new-tau search and conditional candidate

Scope: evaluate the proposed fixed tau=n exp(-delta L), conditional on R085 geometry/probability/history and on the new tau mechanism. Neither this diagnostic nor the exact arithmetic checks certify a C_n theorem. Under the proposed estimates the paid full scalar remainder is

    1/n + exp(-delta L) + exp(-sigma L) + L exp(-(sigma-delta)L),

and tau/k <= rhi L² exp(-[delta-chi(1-sigma)]L). The original global captured threshold and entire cover factor are assumed to remain valid and require independent theory audit.

A preregistered finite parallel batch divided s in [4,5] into 15 equal, non-overlapping interiors. Each worker optimized six coordinates in a bounded differential-evolution run, up to 1302 evaluations, with BLAS limited to one thread. This was a floating log-domain diagnostic, not exhaustive optimization or interval proof. Its process parent was 32075, children 32086–32100, and all child exit codes were zero. It used 17,094 evaluations; 13 workers found model-feasible candidates, while the first two low-s intervals returned no candidate in their bounded evaluations. Across workers 2,580 trials had nonpositive cost surplus. Those counts are observations, not proofs that a continuous interval is infeasible. Best floating display was log2 L≈25.9616. Parameters in the worker receipts are unrounded floats and must not be used as certificates.

The full batch wall time reported by the runner is 0.907 s. The peak over two `ps` samples was aggregate worker CPU 783.6% and summed RSS 894720 KiB (summed RSS includes shared pages, and two samples may miss a higher peak). The summary's summed `worker_cpu_seconds` field is NOT a validated batch CPU total: time.process_time() was not baselined at worker start and individual entries exceed the measured worker wall time. Do not cite the 6.40561 value as actual net CPU seconds. The first attempted background launch PID 32041 did not persist after its shell exited, produced an empty run.log, and did not run the batch. The later foreground process 32075 produced the results. No process was artificially prolonged and no memory was filled to increase load.

After the batch, a fresh explicit rational candidate was registered: s=5, chi=201/200, sigma=5013/10000, delta=12532/25000, qmax=49999/50000, gamma=3193/2000, c=9/10; trim=1/10000000, e=1/100000, Ret=999/1000, m=4007/2500, m_face=247/250, and price upper=3479/5000. The scalar allowance gamma-price-c=7/10000 is exact. At L0=2^26, the actual Fraction check PID 32148 failed exactly one of 56 finite arithmetic gates: the source-face Chebyshev bound was 0.0010339126 > 1-Ret=0.001. This failed receipt is retained. The same tuple was separately registered at L0=2^27 before PID 32161, and all 56 gates passed:

| Payment | Exact-envelope decimal display | Required |
|---|---:|---:|
| source-face Chebyshev | 0.000516955743 | <=0.001 |
| weighted source cap | 3.47023e-21 | <=e |
| new tau/k | 1.72738e-35 | <=e |
| full new scalar remainder | 3.84058e-34 | <0.0007 |
| complete (1-e)Gamma | 1.6010050502 | >1.5965 |

The exact log enclosure gives mu(5)<0.935046258<chi and d_mu=(chi-mu)/4>0.01748843; (mu+2d_mu)/(1-e)<m_face<chi. Strict primary conditions, original target mass, 5 target-event probabilities, full matching/count bounds, individual J, weighted cap, radial calibration, positive B-gamma n, complete deficit and source/target history separation are checked using N=L^16/16!, 96-term atanh log enclosures and S16 exponential prefixes. Candidate files contain exact Fraction bounds and all gate booleans; decimal values here only orient the reader.

For all L>=2^27, N grows as L^16 and x as L: source radius N/(sx) grows; polynomial(L)/N probability bounds, epsilon, log-product, profile transport, count and radial errors decrease. The source-face 10x/(ax-p)^2 decreases for ax>p. The weighted-cap exponent (nu-1/N)eta x/(1+e) increases. For the tau payment, yS16'(y)-2S16(y)>0 for y>=4 follows from the stronger positive lower bound -5/2-3y/2-y²/2+y³/12+y⁴/16. The exact quadratic coefficient for yS16'-(5/2)S16 is -1/4; the displayed -1/2 is conservative. Hence L²/S16(tau-rate L) decreases. For the new scalar remainder, 1/N and the two exponential prefixes decrease, while L/S16((sigma-delta)L) decreases once yS16'(y)-S16(y)>0, which holds for y>=2 by the positive quadratic term. The B>=nL/4>gamma n and a_alg>=8n^4 size/displacement inequalities persist. These are tail comparisons for the assumed envelopes, not proofs of their geometric or operator origin.

All files are under this owned directory. `registration.json`, `run_batch.py`, 15 worker receipts and `summary.json` preserve the batch. `candidate_registration.json`, `candidate_verify.py` and `candidate_result.json` preserve the 2^26 failure; their `candidate2_...` counterparts preserve the 2^27 pass. Exact scripts can be run directly but overwrite PID-bearing receipts; preserve originals first. Original repo, earlier rounds and paper were not changed.
