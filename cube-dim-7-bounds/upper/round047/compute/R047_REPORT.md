# R047: complete cross-radius pinned-face sum diagnostics

This batch tests the complete face sum for two source radii, including same-face and distinct-face cross terms. It is an 80-digit numerical diagnostic, not an interval enclosure or a proof of the R045 weighted-source theorem.

For a target face with pinned coordinate w₀=r and n−1 iid uniform coordinates W∈[−r,r], write J_j=J_{ρ_j} and s_j(w)=N_j(w)/ℓ_{*,j}. Define

\[
\mu=\mathbb E[J_1(W)J_2(W)],\quad
A=\mathbb E[s_1(W)J_2(W)],\quad
B=\mathbb E[s_2(W)J_1(W)],\quad
C=\mathbb E[s_1(W)s_2(W)].
\]

Expanding the double sum over source-face indices, with k=n−1 and pinned quantities indicated by p, gives

\[
4\mathbb E(S_1S_2)=s_{1p}s_{2p}\mu^k
+k(s_{1p}J_{2p}B+s_{2p}J_{1p}A+J_{1p}J_{2p}C)\mu^{k-1}
+k(k-1)J_{1p}J_{2p}AB\mu^{k-2}.
\]

The four groups are pinned/pinned, pinned/free, same free coordinate, and distinct free coordinates. This expansion was independently checked in a four-coordinate discrete microexample by enumerating all eight free-coordinate configurations; the difference was approximately 1.69×10⁻⁸⁰. The E(J−1) primitives were also compared to direct piecewise quadrature in a separate small test.

The registered cases are m=log₂n∈{32,64,128,256,512} and d∈{2⁻⁸,2⁻¹²,2⁻¹⁶}, giving 15 independent tasks. Each uses n=2^m, x=(log n)/2, t=n/x, ζ=n⁻³ and h=2⁻⁴⁰/(xd). For each task, u/d=1/4+i/30 for i=0,…,15 and θ₁,θ₂=−1+2j/15 for j=0,…,15, so all 4096 registered triples were evaluated. The fixed formulas are r=1+u, q_j=r−d+θ_jh, ρ_j=tq_j, R_j=q_j+d−h and g_j=d−h+ζρ_j. These parameters are diagnostic and are not claimed to be the R045 theorem's selected asymptotic regime.

One-dimensional first moments and face indicators are evaluated through exact piecewise primitives; the shared-coordinate cross term E[(J₁−1)(J₂−1)] is computed by 80-digit piecewise quadrature. Exchanging source radii gives the same cross integral, so each task uses 2176 unique cross-radius integrals for its 4096 ordered triples. The run completed 61440 triples, 32640 unique cross integrals and 152064 nonempty quadrature segments. Every displayed normalized value E[S₁S₂]/n² was far below 4096; no numerical counterexample occurred. The stored 32-digit values ranged from approximately 0.9999999890 to 1.0. At large m, values displayed as exactly 1.0 can reflect finite precision and do not establish a strict ≤1 inequality.

To avoid losing a 1/ρ-sized deviation at n=2^512, the computation forms μ−1 directly and uses log1p before multiplying by n; it does not construct a length-n vector. Two already-computed triples, m=32 and 512 with d=1/256 and opposite θ values, were separately recomputed from four direct 250-digit one-dimensional moment integrals; each agrees with its stored 32-digit result within 4×10⁻³³. These checks are numerical and do not replace a uniform analytic or rigorous interval argument.

The 15 workers finished in 61.06 s, with 385.03 s summed worker CPU, 15 distinct PIDs and peak worker RSS about 104 MB. A live `ps` sample at 5.15 s recorded 1461.3% aggregate CPU across all 15 workers and about 542 MiB total RSS; 13 samples were recorded, with aggregate peak RSS about 913 MiB. These are measured samples, not a guarantee of sustained full utilization. Full parameters, values, failure/stop states, source hashes, PID/time/RSS and monitoring are in `registration.json`, `case_*.json`, `summary.json`, `monitor.jsonl`, `micro.json` and `check_moments.json`. No MATLAB, old parameter grid, n-dimensional vector, or production Git operation was used.

Reproduce from clone root:

```
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round047/compute/cross_face.py --run
PYTHONDONTWRITEBYTECODE=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-7-bounds/upper/round047/compute/check_moments.py
```
