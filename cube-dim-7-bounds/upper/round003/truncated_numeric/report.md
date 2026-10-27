# Truncated kernel finite diagnostics

R003 tested the already frozen 8192 original Dirichlet shapes, not a new shape sample. The independently proved formula is Hε=KD E[d Q⁻³ I_{Q/(Q+ε²)}(3,65/2)] with Q=xᵀA⁻¹x, d=det(A)⁻¹/² and KD=Γ(65/2)/(π^(7/2)Γ(32)). The executable uses SciPy's regularized beta and the cached C·D=KD. The cache file hash and structure were checked before execution. The cache inverse and 32 Cholesky solves at each of the first four points agreed to relative 9.41e-16; these checks reuse the warmup inputs. The minimum computed Q was 8.74e-8; every term was finite and positive and monotone in epsilon. The separate R003 proof provides an analytic Beta derivation and a high-precision fixed-input check. This run did not repeat a full-G optimizer.

| epsilon | smallest of 4096 Sobol empirical means | 8-block range / mean at that point | median of 8-block relative ranges across 4096 | largest such relative range |
|---|---:|---:|---:|---:|
| 1/1024 | 2.592072869e-6 | 0.874–1.124 | 0.420 | 5.180 |
| 1/256 | 2.586716218e-6 | 0.873–1.126 | 0.396 | 3.929 |
| 1/64 | 1.456854417e-6 | 0.898–1.074 | 0.190 | 1.005 |
| 1/16 | 4.403688900e-8 | 0.863–1.152 | 0.357 | 0.462 |

The eight blocks are the original contiguous blocks of 1024, and relative range is (max block − min block)/the full empirical mean. The selected minimizers, their binary64 coordinates, all 8 block means, eight previously selected exact-float S3 points and their block means are retained in results.json. These are pointwise descriptive summaries, not bounds on true Hε or its minimum. In particular, 4096 finite face points do not cover the face, and even stable blocks cannot certify a shape expectation error. Near S3, the finite shape sample can miss rare contributions; the old S3 values are diagnostics of truncated Hε only.

The fixed finite budget was exhausted exactly once: 4104 points, 4 epsilons and 8192 shapes, 134,479,872 kernel terms. All 257 chunks have terminal records; first 16 points were the warmup and are included, then 15 processes with BLAS threads fixed at 1. Wall time 1.809 s, sum of task CPU 11.896 s (about 6.57 core-equivalents averaged over total wall time, not 15); peak reported RSS for any one process was 124,993,536 bytes. This is not simultaneous aggregate memory or system-wide peak utilization. No resources were filled artificially. Input generation, process setup and ledger writes account for part of wall but are not all included in task CPU. Single-process peak RSS is the OS high-water mark, not a sum of allocations.

For the next stage, pursue the separate theory's deterministic positive-mass cell / interval kernel pilot: explicitly bound cell mass and shape variation as well as spatial coverage. This run suggests the epsilon tradeoff is material and block behavior is heterogeneous; the larger cutoff suppresses H substantially while the smaller cutoffs show large block spread at some points. It does not select a provably optimal cutoff or establish feasibility of a full-face certificate. In particular no deterministic expectation error, lower g, or C7 bound is claimed.
