# R023 finite threshold-selection parameters

The eight registered pairs give their smallest strictly separated value at **a=36, τ=7**, yielding

    C7 < 46650446233303 / 10^12 = 46.650446233303.

The underlying 512-bit interval is approximately 46.6504462333024540772. It is strictly below the complete R022 a=63 interval, approximately 77.4333403222526020323. This is a finite-grid result: it does not claim optimality over all real a or τ.

The bound uses the separately proved original-operator inequality

    C7 ≤ (1+A/τ)(τ+B),
    A=(1+2/a)^7,  B=(1+1/a)^7(1+7 log a).

The proof applies to general L¹ inputs and every positive common radius. The weights in the finite selection remain the original increment-volume weights; τ changes eligibility and the resulting charge factors. The scalar computation does not establish that covering theorem by itself.

| a | τ=6 | τ=7 |
|---:|---:|---:|
| 32 | 46.842738226839 | 46.702391397400 |
| 36 | 46.748940177898 | **46.650446233303** |
| 40 | 46.767339437887 | 46.699178854626 |
| 64 | 47.743137383160 | 47.751907606611 |

Table entries are rational upper ceilings on the 10^-12 grid, not rounded-to-nearest approximations. Each of the eight records in terminal.jsonl also retains the full outward 512-bit dyadic interval. The unique winner's upper endpoint is below every other row's lower endpoint. The independently derived analytic comparison at a=64,τ=6, namely K<1040858/21375<49, is also confirmed against its full interval; it is not the source of the improved numerical winner.

R022's completed resumed result and its manifest were verified by hash. Its initial failed run remains separately preserved; it is not represented as the successful result. The eight rows share four computations of A and B and four log evaluations, and each matches the algebraically equivalent A+B+τ+AB/τ form by interval overlap. The scalar run used 0.00223 s wall, 0.00210 s CPU and peak RSS 30,638,080 bytes in one process; no CPU saturation is claimed.

The checker repeats only the same eight scalar rows at 512 bits, using eight further log calls for verification, verifies exact per-row ceilings, finite-grid ordering, the analytic comparison and the complete old interval. It evaluates no old a=64,τ=1 baseline, extra parameter, optimizer, MC, Green phase or lower-bound event. No MATLAB, commit or push was used. The theorem and its attribution to the original Stein–Strömberg selection mechanism remain in the separately hash-locked proof.
