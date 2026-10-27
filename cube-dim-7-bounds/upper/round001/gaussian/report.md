# R001 Gaussian-mixture complete-boundary moment result

Completed 2026-10-27. This is a seven-dimensional re-proof and new finite run of the inherited n=6 R123 method, not an original mechanism. Fixed object: Lebesgue centered axis-parallel seven-cube, one common radius, general L1, all positive scales, and the original equal-weight body64 common-time Cauchy generator. No new C7 bound is established.

## Strict result

With q=7/6, Sigma=sum_{v in body64} S_v vv^T and E exp(-lambda S)=exp(-sqrt(lambda)), the exact Gaussian-mixture potential is

G(x)=[1/(2pi^(7/2))] E[det(Sigma)^(-1/2)(x^T Sigma^(-1)x)^(-3)].

For I=integral_[0,1]^6 G(1,y)^(7/6)dy and M=E det(Sigma)^(-1/14), full-boundary Minkowski, before using symmetry of the averaged G, proves

I <= M^(7/6)/(105*2^(25/6)*pi^(13/12)).

The exact radial reduction is M=D E Z, D=2 Gamma(65/2)/Gamma(32), Z=det(sum vv^T/P_v)^(-1/14), P~Dirichlet(1/2)^64. Cauchy--Binet and weighted AM--GM prove 0<Z<=1/64 and the analytic endpoint

M <= M_Gamma=(1/4)[Gamma(65/128)/Gamma(1/2)]^64.

The deterministic checker used 256-bit Arb with exact rational arguments. Its enclosing balls give

- M_Gamma = 0.094562363737025712424627825917389342081619594676634688671371258372615934554 +/- 4.09e-76.
- I_Gamma upper expression = 9.7936533902956432088164651795416774612676259744431036453303013898515961914e-6 +/- 4.60e-80.

The ball is an enclosing evaluation of the exact analytic upper expression; its upper endpoint bounds I. The symbolic coefficient ratio is exactly 1. These statements require no empirical expectation estimate.

For c>0 the whole-space positive-part mass is exactly (6/7)c^(-1/6)*896I. A separately certified global cube floor g_L yields, through the accepted original-generator positive-source transfer, U=7^(7/6)I/g_L^(7/6). The present file supplies the complete I estimate and finite L1 mass; it does not reprove that transfer or supply g_L. Using the coarse analytic M_Gamma would require g_L at least approximately 5.56470466584e-6 for that expression to reach 128; this is a conditional target, not an established floor or operator result.

## Revised fixed point audit

Before any shape run, the parent revised the point component to avoid evaluating seven anchors at which G is infinite. Every a_k=(1,1^k,0^(6-k)), k<6, is the mean of two body vectors; a_6 is a body vector. Exact SymPy rank checks show that the remaining body vectors span R7 after removing those one or two columns. The convolution lower estimate p_t(a_k)>=const*t^(2r-7), r=1 or 2, proves divergence of the potential. The seven anchors were audited, with no Monte Carlo spent on them.

The registered replacement consists of eight rational face points x_j=(1,(j+1)/(j+7),...,(j+6)/(j+7)), j=0,...,7. Their seven distinct positive absolute coordinate values exclude all spans of at most three body lines, since such combinations have at most four absolute coordinate values. All eight were evaluated with exactly the same shapes used for M. No point values are treated as a global g lower bound.

## Actual one-batch numerical evidence

The registration SHA256 is `4c5b636158930bf5b5039f858828f9e44c32015b04ec09acee0a94bbdf2785f9`. The sole shape run used 8192 samples, eight blocks of 1024, NumPy PCG64 with block SeedSequence([202610277,b]), b=0,...,7. R was integrated out exactly. A single process ran with all five BLAS/OpenMP thread limits set to 1. There was no automatic expansion, extra grid, new seed, clipping, dropping or resampling.

Scientific draw/evaluation loop: actual wall 0.280443500 s; actual process CPU 0.272683000 s. Full tool-launched process elapsed 0.566097250 s, including imports and file writes. The receipt records the loop timer and exact environment. NumPy 2.5.4, SciPy 1.18.1, Python 3.15.0, macOS arm64.

All 8192 rows passed Gamma positivity/finiteness, eigenvalue SPD, Cholesky, positive finite output, and Z<=1/64 checks. Zero counts for every failure category; no 80-digit recomputation was triggered. Maximum scaled-matrix condition estimate was 6.24836869457e9, minimum scaled eigenvalue 1.12029240682e-9, and maximum Cholesky relative residual 3.19373520813e-16. The high condition numbers remain a limitation of floating point diagnostic samples; these are not interval-certified sample values.

Observed Z mean 0.002821860941509685, empirical standard deviation 0.0011190970982084244, empirical SEM 0.000012364392921083821. Thus M_hat=0.031801250614156534. The eight M block means were

| block | M mean |
|---:|---:|
| 0 | 0.03140861715408672 |
| 1 | 0.032003060921412425 |
| 2 | 0.031422039915456196 |
| 3 | 0.031497851306177695 |
| 4 | 0.031820868582931074 |
| 5 | 0.03211348429005199 |
| 6 | 0.03183786170505116 |
| 7 | 0.03230622103808502 |

Using the exact radial factor for G, the fixed point diagnostics were:

| j | G estimate | empirical SEM | max single sample contribution |
|---:|---:|---:|---:|
| 0 | 9.36128496156e-6 | 4.63321487975e-7 | 0.00162120678 |
| 1 | 9.29533426921e-6 | 6.20236251574e-7 | 0.00218607452 |
| 2 | 1.04030478566e-5 | 8.77618099696e-7 | 0.00350992687 |
| 3 | 1.25057045955e-5 | 1.25965145417e-6 | 0.00528707680 |
| 4 | 1.55758196233e-5 | 1.79900739488e-6 | 0.00741422016 |
| 5 | 1.96714196099e-5 | 2.53390172876e-6 | 0.00981689537 |
| 6 | 2.48938574910e-5 | 3.50781513364e-6 | 0.01438907061 |
| 7 | 3.13673893481e-5 | 4.76809730965e-6 | 0.02070107572 |

These point observables are visibly heavy-tailed. The reported SEM values are descriptive empirical summaries, not proved confidence intervals. Sampling all points with common shapes also correlates their estimates. The finite-point minimum, 9.29533426921e-6, does not bound G from below over the cube.

The complete block means below are in units of 1e-6:

| block | j=0 | j=1 | j=2 | j=3 | j=4 | j=5 | j=6 | j=7 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 8.525216 | 7.923595 | 8.309272 | 9.581540 | 11.757020 | 14.955001 | 19.351367 | 25.153801 |
| 1 | 9.531002 | 8.798865 | 9.007946 | 10.002785 | 11.774539 | 14.382595 | 17.909422 | 22.443517 |
| 2 | 12.139556 | 14.099639 | 17.735819 | 22.885712 | 29.606553 | 38.059061 | 48.446393 | 60.979831 |
| 3 | 8.036823 | 8.469369 | 10.223052 | 13.195808 | 17.421744 | 23.014645 | 30.139621 | 38.993686 |
| 4 | 9.946746 | 9.384393 | 9.927309 | 11.346699 | 13.496007 | 16.308242 | 19.749787 | 23.801938 |
| 5 | 8.531266 | 8.111860 | 8.941509 | 10.758724 | 13.471602 | 17.061885 | 21.546691 | 26.966740 |
| 6 | 9.515307 | 9.393973 | 10.461643 | 12.449935 | 15.317710 | 19.136927 | 24.045582 | 30.224280 |
| 7 | 8.664363 | 8.180982 | 8.617834 | 9.824434 | 11.761382 | 14.453000 | 17.961998 | 22.375321 |

The block max/min ratios are respectively 1.5105, 1.7794, 2.1345, 2.3885, 2.5182, 2.6462, 2.7051 and 2.7253. One largest sample accounts for respectively 2.114%, 2.871%, 4.119%, 5.161%, 5.811%, 6.092%, 7.056% and 8.056% of the total sample sum at each point. These figures are read directly from the saved run, without any additional draw.

The empirical complete-moment expression I_M(M_hat)=2.74657580200e-6. Inserting only the minimum of these eight sampled point estimates into the conditional cost expression gives 19.7286236802 (0.154129872502 of 128). This diagnoses the cost expression at the specified eight points only. The actual global g may be smaller by an order of magnitude or more; therefore this number must not be called a competitive upper bound. It is not a C7 bound, a proof of competitive feasibility, or an exclusion of any parameter range. Both an expectation certificate stronger than sampling if needed, and a global g certificate, are distinct tasks. No further samples were taken.

## Files, verification and reproduction

- `proof.md`: complete seven-dimensional derivation, including the correct full-boundary factor 896, Gamma/Dirichlet reduction, positive-part mass, exact anchor divergence and the finite-point limitations.
- `registration.json`: frozen pre-run design, point revision, exact rational inputs and source hashes.
- `probe.py`: guarded single-batch program; existing result files cause it to refuse repeat science execution. Reproduction would require a separate empty output copy of the frozen sources and registration, deliberately not run again here.
- `observables.npz`: retained 8192 Z values, 8192x8 point observables, all per-shape eigen/condition/residual diagnostics and method flags. No generated shape is omitted.
- `results.json`, `receipt.json`: actual summaries, per-block results, failure counters, timestamps, runtime, package versions and run hashes.
- `check.py`, `check-results.json`: deterministic symbolic/Arb checks and checks of the saved finite results. Re-run with `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check.py` from this directory; it takes no random draws.
- `cache_rebuild.py`, `cache_rebuild-registration.json`, `cache_rebuild-receipt.json`: separately authorized engineering recovery of the original shapes, described below.
- `sources-hash.json`: final deliverable manifest and read-only method/source provenance.

The original science program retained observables but did not persist P or inverse matrices. This omission was disclosed after the run. The parent explicitly authorized one deterministic replay of the **same** frozen RNG stream solely to restore a reusable cache. Before replay, a separate recovery registration locked the original source, results, receipt and observables hashes. Every one of the 8192 rows then reproduced Z, all eight T values, condition, min/max eigenvalues and Cholesky residual **bit for bit**; any discrepancy would have stopped recovery. The original scientific code, results and receipt were not edited. This is not another independent 8192-shape experiment, and it added zero independent science samples and zero new points.

The cache is outside Git at `/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/cache/gaussian-round001/original-8192-shapes.npz` (8,757,232 bytes), SHA256 `6dfa48bf8954a7d02270af5bda9c1d33d4fcaf0ba4f5e486589186b7f7c09552`. It contains original P, A inverse, det(A)^(-1/2), logdet A, scaled Cholesky, its A scale, and the exact-formula-linked floating radial D and potential C factors. R was analytically eliminated, so A inverse is the appropriate shape inverse; no fictitious sampled R or Sigma was introduced. Future arbitrary-point candidate evaluation can use these saved factorizations directly without reconstructing shapes. The recovery itself used actual wall 17.372359417 s and process CPU 17.093915 s, separately from the scientific loop above. No g optimization or further science was started.

Research deliverables are confined to `cube-dim-7-bounds/upper/round001/gaussian/` in the dedicated clone, with the explicitly authorized engineering cache in the stated workspace scratch path. No commit, push, MATLAB call, recursive agent work, other-dimension write, or source-project edit was performed. Nothing was deleted or cleaned outside these owned paths.
