# R005 Radau1792: proposed frozen finite certificate for Terra

2026-10-27. Sol theory source: `merge4_proof.md`. This owner's `upper/round005/theory/` contains no numerical execution. Terra owns a separate numeric directory. Fixed original n7 body64/common-time/unit-weight Cauchy model, same8 Hadamard blocks. No MC, old optimizer, MATLAB, g/C7 claim, commit/push or automatic scope extension.

Target: improve the complete I budget by four-block merging and a certified two-moment Radau correction. **Reuse** the complete exact R004 q(1/8) interval from its hashed certificate.json; do not recompute it. q(0)=1 by probability normalisation, never a Mellin integral at0.

New positive real pair exponents and fixed grids:
- r=1/16: a=1/32, halfwidth1/64, Cauchy radius3/128,1024 panels.
- r=3/16: a=3/32, halfwidth1/32, Cauchy radius3/64,512 panels.
- r=1/4: a=1/8, halfwidth1/16, Cauchy radius3/32,256 panels.
Every grid covers[0,32], centres(2j+1)halfwidth. Total1792 unique(r,index) panels. All:512-bit acb_series outward arithmetic, degree128/order129, Gamma tail recurrence order16. Each h/rho=2/3. Raw Cauchy error=32 M_star(r)(2/3)^129/(1-2/3); normalise bypi Gamma(r). Tail=m(a)^2Gamma(a+16)^2/[pi Gamma(r)·31·32^31]. M_star(r)=Gamma(a-rho)Gamma(a+rho)m(a-rho)m(a+rho). These constants are per-r, not copied from1/8.

Before scientific run: code/config/source SHA256, existing R004 certificate/theory hash and versions locked. Interface probe changes no panel ledger. Pilot=first4 panels of each of3r,12 total; compare direct constant/series evaluation and Gamma recurrence, validate series order/radii/real coefficients/finiteness. Mathematical identities are independent of these overlap diagnostics.

After pilot pass:15 processes, BLAS/OpenMP1, remaining1780 panels only; reuse12 pilot results. No wall-time cap. Do not rerun completed terminal keys. Record actual CPU/wall/peak memory/occupancy; no false full-load claim.

Assembly: exact coverage of1792keys, correct h/rho/precision/degree perr, complete q intervals, all source hashes unchanged. Use merge4_proof equations(9)-(11): mu_L>0, nu_U>0 finite and nu_U>=mu_L^2; compute downward0<delta_L<1, M_upper=(q_a_upper)^4/2*(1-delta_L)^2, I_upper=M_upper^(7/6)/(105·2^(25/6)pi^(13/12)). Compare to complete R004 budget-expression lower endpoints. Retain exact finer I upper endpoint and a declared outward rational rounded endpoint.

Acceptance: all panels/parameters/hashes valid, strict delta and quantitative budget comparison separated. If failure/interval overlap, report exact issue; no midpoint certificate, no automatic precision/degree/cutoff/exponent extension. Stop after assembly and main audit. No g or C7 upper bound follows solely from this I certificate.
