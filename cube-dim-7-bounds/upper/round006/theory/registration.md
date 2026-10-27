# R006 four-moment Radau4096 proposed finite registration

2026-10-27. Source `radau4_proof.md`. Sol theory only; **do not start numerics until explicit main-owner release**. Terra's later numeric directory owns execution. Original7D body64/equal-weight/common-time model unchanged. No MC, MATLAB, old optimiser, g/C7 claim, commit/push or automatic extension.

Target: three-node Radau upper expectation for g(W)=T16(1+W)^(-1/4), W=cosh(Delta/32)-1 under the same a=1/8 tilt; use moments through4, roots of w²+alpha w+beta and strictly positive weights. Reuse q(1/8) R004 and q(1/16),q(3/16),q(1/4) R005 complete hashed endpoint certificates. q0=1 exactly, no MB at0.

New positive real pair exponents r=1/32,3/32,5/32,7/32. Each T=32,Ntail16,512-bit acb_series,degree128/order129. Set a_r=r/2, rho=3a_r/4. Halfwidths1/128,1/64,1/32,1/32; panel counts2048,1024,512,512; eta ratios2/3,4/9,8/15,8/21. Total4096; exact centres(2j+1)halfwidth. Per-r raw Cauchy error32 Mstar eta^129/(1-eta), Mstar=Gamma(a-rho)Gamma(a+rho)m(a-rho)m(a+rho). Normalized tail=m(a)^2Gamma(a+16)^2/[pi Gamma(r)31·32^31]. All values exact rational except enclosed transcendental factors.

Before run: source/config/code/reused endpoint SHA256 and versions locked. Interface probe no scientific panel claim. Pilot first4 panels perr,16 total, direct/series constant overlap, Gamma recurrence, seriesorder, real coefficients/finiteness and parameter checks. Reuse16 after pass. Explicit main release required for15 single-thread workers to compute remaining4080; no wall cap, finite parameter cap is4096 unique keys. Record actual CPU/wall/RSS/occupancy, not assumed saturation.

Assembly: complete4moments via R_l and exact formulas in proof; certify D>0, alpha<0,beta>0, discriminant>0; stable ordered root intervals positive/disjoint; all3weight intervals positive and mass/moment residual diagnostics. Mathematical exactness/global majorant rests on proof, not naive interval interpolation. Compute correction lower via positive lowerweights and lowernodes, require0<delta6_L<1; buildM/I upper, compare against accepted R005 deterministic budget-expression lower endpoints, retain exact finer endpoint and outward rational I upper. No confusion with actual I/M lower bounds.

If any separation/hash/coverage gate fails, report and stop. No automatic higher precision/degree/cutoff, shifted roots, clipped weights or extra exponent. Stop after certificate/main review. No new g/C7 claim from I alone.
