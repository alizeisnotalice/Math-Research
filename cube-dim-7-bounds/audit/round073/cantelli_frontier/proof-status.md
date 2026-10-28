# R073 preparation — full theorem PENDING

Historical preparation stage: superseded by the full-chain ACCEPT in
review.md for frozen source 50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b.
The stage receipt was not rerun or relabeled as a frozen-source run.

The author has supplied a fixed tuple and parameter certificate, but no
frozen full proof has yet been read. Nothing here activates a C_n bound.
R071's two actual run identities remain in its own directory; their PASS
results are not evidence for this new tuple.

## New sensitive interfaces

Here qmax=11569/12800 differs from 1-sigma=231/256. The cap exponent
eta=(1-qmax)/(2qmax) and its logarithmic band mass must use qmax, whereas
x=(1-sigma)L and the scalar log price use sigma. The cap rate cS and tau
rate cT are distinct; their positive signs and prefix domains require
separate checks. The R071 equal-rate simplification cannot be imported.

At L0=8192, delta/h=12. The source-band endpoints, complete matching
interval and active R>1 therefore have more room than in R071. The count
drift 2h+1/n fits the new eN budget. The product budget is 6Lh=1/128,
not zero, and must be included in the Ulo/Uhi energy ratio.

Let a=chi-.04422/(1-eJ)>0 and p=.183. The fixed Cantelli retention is
computed as (a x-p)^2/[v x+(a x-p)^2], v=.00531. Its base value is not
alone a dimension-uniform proof. On x>p/a,

    d/dx [(a x-p)^2/x] = a^2-p^2/x^2 >0,

so this retention increases over the entire larger-x tail. Each face
must first satisfy the same mean/variance bounds; only tail probabilities
can then be mixed.

The tight cap loss must use n log(qmax/qmin), not the old loose source
mass. If eta*x/n<=1/1025, the elementary exp(u)<=1/(1-u) gives the
finite-n correction 1025/1024. This is valid for u>=0 and u<1; dropping
the correction or using a rounded decimal log is not allowed. It must
be extended with the appropriate prefix monotonicity to all L>=8192.

## Actual new narrow checks

preflight.py used a different n lower envelope L^32/32!, eight positive
atanh terms with an exact geometric log tail and exponential P40. Actual
run PID 1890: seven compressed Fraction groups PASS. These check the
distinct qmax/x roles, trim/count margin, fixed Cantelli retention and
monotone domain, tight target mass and gamma ratio, finite-n cap factor,
separate cap/tau base bounds, and the exact price/slack identity.

The author 53-gate verifier, R071 programs and older passed baselines
were not run. This is not yet a full root-cost or operator audit.

## Pending full audit

Read and lock the frozen actual proof; check all local guards and six
target event probabilities anew; reconstruct whole-matching profiles,
pin handling/J denominators, per-face Cantelli then clipping; verify
radial degrees and exact cap/tau whole-tail budgets; retain ca^2/ca in
the genuine same-history finite-a deficit and general L1/all-radius
cover transfer; pay every positive scalar A/cross/self buffer beneath
the 1/5000 slack. No candidate frontier/model optimum or display
decimal may substitute for any of these steps.
