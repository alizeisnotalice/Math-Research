# Astra015: one Cauchy Young family, staged exact certificates

Parent authorization after locked Astra014, 2026-10-26. This directory alone is writable; no previous evidence or other agent directory changes; no commits.

Fixed alpha=1, beta=1/2, same 16 unit-weight unnormalized body directions and shape law as014. Fixed u=(W+2)/12, W=sum_(i<j) r_ij x_i x_j. On the cube W>=-2, hence u>=0. The only Young family is Z=a*u with a>=0. No parameter/point/grid search.

Stage 1, necessary corner only: calculate exactly Eu², E[d u²], E[B u³] for x=(1,1,1,1,1). The true positive term U=E[sqrt(d)u²] is at most sqrt(E[d u²]*E[u²]). Optimize a algebraically, giving corner-family upper budget c1*(E[d u²]*E[u²])^(3/2)/(E[B u³])². Compare it rigorously against 014's required g for the current I expression and against 1/10000. This is only an upper bound on this certificate family at a necessary point, never an upper bound on Green. If it excludes the required target, stop and record that exact fixed-family obstruction. No high determinant powers or full-space production after such an exclusion.

Stage 2 is authorized only if Stage 1 does not exclude the current required target: calculate new-beta corner E[d^k u²], k=1..4. Choose a unique two-node Gauss pair for the weighted measure d*u² dR, rationalize its nodes, and form a cubic Hermite lower interpolant of d^(-1/2); multiply by d for a valid sqrt(d) lower bound. Direct cubic Hermite of sqrt(d) is an upper bound and will not be used. Enclose all polynomial coefficients downward, retain possible negative coefficients, choose a single rational a near U_lower/V, and stop early if even this candidate's necessary corner is inadequate to warrant full-space work.

Stage 3, only if warranted by Stage 2: produce complete new-beta spatial moment coefficients for E[d^k W^n], k1..4,n0..2 and E[B W^n],n0..3. Preserve orbit multiplicities and B-shape correlation. Rebuild the one chosen lower polynomial and verify all 81 face Bernstein coefficients with rational arithmetic. Do not infer a global value from the corner alone. Geometry/source expressions may be reused, never beta=.75 moment values.

All stages use new beta=1/2 cache and source fingerprint. The copied GMP backend changes C_n to product_(j=1)^(n-1)(2j-1), D_n=product_(j=0)^(n-1)(2j+1), with C1=D1=1. Five-bit keys, max total degree31 guards, arbitrary-precision integer accumulators and source-specific cache header remain. Stage1 max degree7; any Stage2/3 max22. Independently verify the single-edge moments E r^(2j)=1/(2j+1) and odd0. Include at least one multiedge rational Stieltjes coefficient comparison so the test is not only one-dimensional.

Each process gets a 120-second checkpoint. The moment cache survives a backend checkpoint; stage accumulators are rebuilt from an empty state after interruption so no partial sums are duplicated. Stage1 is small. For later production, save completed W powers or deterministic stages before proceeding. Record actual runtime, sources and cache provenance. No lifting of fixed-width guard or graph degree limit.

Any positive auxiliary g must be paired with a separate alpha=1 five-dimensional positive-part transfer audit before claiming a C5 result. S1/S2 infinity is permitted on null sets but must not be replaced with four-dimensional continuity. Current C5 upper32 remains unchanged until all proof obligations are met.
