# A009: positive short-potential lower barrier near every S2 plane

Model: the existing alpha=1 equal-weight body32 Cauchy semigroup. This is
not the alpha=3/2 model. Let T=4/125 and theta=T/4096=1/128000. The task
is a lower barrier for K_T=integral_0^T p_t near S2; B/F own the source
and weak-type transfer and those are not re-proved here.

Only this owned directory may change. No kernel evaluations, MC, spatial
grid, MATLAB, agents, commits, pushes, or changes to completed work.

Use exactly five canonical pair geometries (Hamming k=1,...,5). Use the
self-contained fourth-moment Rademacher bound c=8sqrt(6), not an unverified
sharp Khintchine constant. The orthogonal R4 quotient after deleting an
additional variable has lower coefficient 7sqrt(6).

Use the proved weak truncation tail Ctail*t^-3/L. The stronger suggested
t^-2/L^2 bound is not assumed. Choose kappa=sqrt(5/27) so the Fourier
quadratic lower bound loses at most a/4, and tau<=aL/(4Ctail) so planar
truncation loses at most another a/4.

Eight scalar budgets only: L in {1/4,1/2,5/8,3/4,7/8,1,3/2,2}, with
tau=min(T,aL/(4Ctail)). These fixed physical cutoffs balance tau proportional
to L against the pair density loss proportional to (sqrt(6)+L)^(-4).
Compute an explicit rational rho threshold, preferentially 1/80000 if
strictly validated, otherwise a downward rational threshold from the best
of these eight budgets. Do not add a ninth parameter case.

Every step must use positive coarea/Tonelli or density domination, never
conditioning on a zero-probability fiber. Preserve the pair Jacobian
1/sqrt(det Gram). Prove that the uniform finite-plane barrier includes
plane intersections and S1. Record the remaining compact complement as
unproved, and stop after the finite batch.
