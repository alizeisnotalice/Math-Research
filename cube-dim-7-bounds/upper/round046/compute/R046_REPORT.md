# R046: inverse-column-weight one-dimensional diagnostics

These are numerical tests of the one-dimensional row formula, not a rigorous interval proof or a C_n bound. Set R=1 by scaling; choose g∈{1/128,1/64,1/32,1/16,1/8} and ρ∈{4,16,256}. For each of the 15 (g,ρ) tasks, h/g=i/16 for i=1,…,16 and r=1+2h k/15 for k=0,…,15, yielding exactly 256 registered (h,r) combinations. Identical r values within a task are cached, so no completed integral was repeated simply to fill CPU time.

For each of 97 distinct r per task, the code independently evaluates the defining piecewise z integral for J(r)−1, then integrates the proposed closed form for the face mean and second moment over w∈[−r,r]. All calculations use mpmath at 80 decimal digits, and log1p for the small E term. Breakpoints of both the captured interval and ℓ are explicitly included. The resulting 4365 one-dimensional definite integral calls cover all 3840 combinations. Every direct/closed comparison, mean bound, variance bound, and r=R mean check passed at the registered numerical tolerances. The largest observed absolute difference between the defining and closed-form J−1 calculations was about 1.20×10⁻⁸¹. This is a numerical discrepancy, not an interval error bound.

Over just these registered cases, the largest observed normalized mean deficit ρ(1−EJ)/(hg) was approximately 1.967303, and the largest normalized second moment ρ²E[(J−1)²]/g³ was approximately 2.211449. Both occurred for g=1/128, ρ=256, h=g and r=65/64. These finite-sample maxima do not establish globally valid or sharp constants; the proposed bounds 6 and 31 were not contradicted in this sample.

A separate SymPy check verifies E(u)=∫₀ᵘv/(2−v)dv and the exact algebraic face-mean cancellation EJ=1 at r=R under the normalized-domain assumptions. It does not prove the general uniform second-moment bound. No additional ρ=n/x asymptotic experiments or covariance grid were run.

The 15 cases were executed on 15 distinct worker PIDs, parent PID 79924; wall time was 2.76 s and summed worker CPU was 36.57 s, equivalent to approximately 13.25 CPU cores averaged over that wall period. Peak worker RSS was about 51 MB. The 0.098 s startup `ps` sample contained 15 newly launched workers but their reported instantaneous %CPU was zero; the jobs finished before the 5 s steady-state sample. Thus the CPU-time ratio is reported, not an unobserved instantaneous full-load claim. Per-task values, checks, PIDs, times, RSS and source hashes are in the case files and `summary.json`; the registration fixes every parameter and stopping rule. No MATLAB session, old parameter grid, or production Git operation was used.

Reproduce from clone root:

```
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round046/compute/r046_row.py --run
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round046/compute/symbolic_check.py
```
