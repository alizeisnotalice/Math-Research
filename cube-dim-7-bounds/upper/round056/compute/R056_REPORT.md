# R056: inverse-source moments and factor-two column diagnostics

This is a finite numerical and exact-point diagnostic of the new inverse-source construction. It does not prove the R051/R056 flexible-budget theorem or a new C_n bound.

## Parameter and scale correction

The initial task description briefly wrote g=R−κρ. That conflicts with g>0, r≥R and κρ>r. The coordinator confirmed the correct original-geometry identity

\[
 g=R-(\kappa-1)\rho=R+\rho-\kappa\rho,
 \qquad \kappa\rho=\rho+q-\zeta\rho.
\]

Only the corrected definition was used in the frozen run.

The 15 d values are 1/10, 1/9, 1/8, 1/7, 1/6, 3/16, 1/5, 7/32, 1/4, 17/64, 9/32, 19/64, 3/10, 5/16 and 1/3. For each d, q=1−3d/4+di/254, i=0,…,127; h=d/1024; R=q+d−h; ρ/R∈{4,8,16,64,256,2²⁰,2⁴⁰,2⁸⁰}; n=ceil(64ρ/q), t=ρ/q, ζ=n⁻³ and g=d−h+ζρ. Four targets r=q+d+θh use θ∈{0,1/3,2/3,1}. All seeds and geometrical domain checks are rational. These are finite diagnostics, not an assertion that these choices are uniform theorem parameters.

## Per-face moments and column checks

For the free coordinate of a fixed source face, with ℓ(z)=min(1,(R+κρ−|z|)/(2R)) and φ=−log ℓ, direct 80-digit piecewise integration was compared to

\[
\mathbb E(\phi/\ell)=\frac R\rho (\log\ell_{\rm pin})^2,\qquad
\mathbb E(\phi^2/\ell)=\frac{2R}{3\rho}(-\log\ell_{\rm pin})^3.
\]

The inverse-weighted **per-face** law is normalized by E(1/ℓ); its mean and variance are stored for each source state. This does not substitute a mixture variance. A separate SymPy check confirms the elementary integral primitives. The proof of all geometric hypotheses and their use in the operator still belongs to the theoretical audit.

The column factor is checked with exact fractions on a fixed finite set of points including the threshold from both sides. For κρ>r and the correct target-coordinate ℓ, the algebraic mechanism is: N=2 implies ℓ=1; when N=1, ℓ≥1/2 implies N/(2ℓ)≤1. The sampled exact checks support the implementation, while this two-branch statement is not a replacement for a proof of all geometry and integration bounds.

All 60 dynamic chunks completed: 15,360 distinct source states, 30,720 high-precision quadratures, 61,440 target combinations and 552,960 rational point checks. Every moment-formula comparison and sampled column check passed; no case was stopped for a candidate counterexample. The largest observed direct-versus-formula absolute difference was about 6.59×10⁻⁸³. The largest recorded inverse-weighted φ variance was approximately 0.000722307 at d=1/3, q=3/4 and ρ=1109/256. These maxima are observations on this finite list, not certified global extrema or rigorous interval error bounds.

## Resources and provenance

Parent PID 88450 launched 15 distinct workers. Total elapsed time was 12.80 s, summed worker CPU 167.96 s, and peak individual worker RSS about 226 MB. A `ps` sample at 5.58 s observed all 15 workers active at 1440.0% aggregate CPU, about 1.75 GiB combined RSS; across the three captured samples peak aggregate CPU was 1446.8% and peak combined RSS about 2.83 GiB. These are samples, not a claim that every core was fully utilized throughout the run. No memory was allocated merely to increase utilization.

`registration.json` fixes the exact finite list and code/theory hashes; each `chunk_*.json` records its entire data, completion, PID, CPU/wall, RSS and exact exception location if one occurred. `summary.json`, `monitor.jsonl` and `symbolic_identities.json` preserve run and symbolic-check metadata. Raw outputs remain in this private clone. No MATLAB, R047/R052 grid replay, Git operation or original user file modification was used.

To reproduce explicitly (this repeats the finite 30,720-integral batch):

```
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round056/compute/r056_inverse.py --run
```
