# A051: all 31 unsigned middle domains, short-time upper budget

## Result and exact scope

For the original body32 Cauchy semigroup, let rho=1/64, epsilon=1/1000,
q=6/5, and S_epsilon(x)=integral_0^epsilon p_s(x) ds. For each nonempty
proper A containing coordinate1, put B=A complement and

    D_A={x=(1,y), y in[0,1]^5:
           dist(x,span(1_A,1_B))<1/64,
           mean_B(x) in[1/4,3/4]}.

The complete positive two-parameter decomposition and true two-coordinate
marginal give the certified bound

    ||S_epsilon||_(L^(6/5)(union_A D_A))
       < 4689838525803/1000000000000000000
       = 0.000004689838525803.

This covers all 31 specified middle domains. Their overlaps are allowed,
not declared absent. Both retained-small complementary events are fully
paid. No spatial endpoint, longer time, outside-tube region, global Green
floor, full I integral, or new C6 bound follows. The current C6 upper bound
remains 64. The displayed budget is about 0.0571574 of the conditional
B target 4/48750; that comparison does not establish the other conditions.

## Five geometric types and correct multiplicities

The entries below are approximate displays of complete saved rational
endpoints. The strict union bound above, not these decimals, is the
certificate.

| k=|A| | m=|B| | n_k | joint marginal C_m | bad density M_k | individual short budget U_k |
|---:|---:|---:|---:|---:|---:|
|1|5|1|15.2397366514774|192822.817598020|3.2763373938090e-7|
|2|4|5|15.2001317608193|67325.9938867904|2.8087225339879e-7|
|3|3|10|15.1313163804814|44480.1218451240|2.6400908910704e-7|
|4|2|10|14.9804427041524|42856.6033456880|2.6093618375933e-7|
|5|1|5|14.0401092308809|69683.2300042113|2.6579552389855e-7|

The counts are binom(5,k-1)=(1,5,10,10,5), summing to31. The common function
S_epsilon obeys the positive union estimate

    ||S_epsilon||_union^q <= sum_k n_k U_k^q.

The computed sum is approximately4.03078346107939e-7; taking its5/6 power
gives approximately4.68983852580102e-6. No identification of k with6-k
removes any weight. The first coordinate is fixed in A, which is essential
to the domain geometry.

The k3 row is copied byte-for-byte from complete frozen043/049 intervals,
including their original finite precision uncertainty. It was not
recomputed. Only k1,k2,k4,k5 have new scalar configurations.

## New proof obligations resolved

Theorem.md proves the actual block perturbation |h|<=sqrt(k)rho, the
nearest-line identity including negative B-block residuals, and every
signed-pair plane gap. A different unsigned partition contains a mixed
1/u block; signs cannot improve its projection at a nonnegative point.
A different signed plane on the same partition has an opposite-sign pair
inside a nonsingleton block. These arguments include k1 and m1.

All first and second residual distance bounds and the original-point pair
gap are positive at the fixed rho for all five types. No false uniform
distance is imposed on the shifted third quotient. Its actual Cauchy
parameter density is integrated using the043 averaging identity, so its
isolated zeros remain charged.

The genuine two-column coarea factor has squared Gram4km. It cancels
against the033 projected peak normalization, giving exactly
D_big,k=64 R2,min(k,m)/pi^6. The deleted-parameter complement retains the
(10/9)^3 peak loss, actual counts31/30/29 and Gram bounds, and the complete
26970 ordered-triple coefficient for each small branch. It does not
pretend a frame with removed parameters is tight.

For the new L1 estimate, one fixed coordinate j in B suffices. The true
full two-coordinate marginal, followed by the remaining four-coordinate
positive-orthant factor1/16, gives C_m s^2. This symmetry is applied to the
full density only. The actual nonnegative big and disjoint bad densities
are each dominated by it. The B-coordinate residual norm is
sqrt((m-1)/m), exactly zero when m1. No extra factor m is inserted.

Thus all bad events and all five weights are included. This is an
extension of the049 interface to four further types, not another use of
the old one-coordinate L1 budget from043.

## Actual computation and verification

Registration fixed the scope, sources, parameters, acceptance gates and
20-minute deadline before execution. The producer ran once, PID51121,
exit0, elapsed0.07358712499262765 seconds (tool chunk449644). It completed
seven configurations: two exact structural/source checks, four new Arb768
types and one union computation. Kernel, native, full spatial/fan/frame
enumeration, grid, optimizer and MATLAB call counts are all zero.

The eighth configuration is a rational-only saved-data audit. It imports
no FLINT/Arb and makes no new kernel or parameter evaluation. For each of
the four new types it verifies square-root gap inequalities, uses a fixed
96-term positive atanh series with an explicit geometric remainder to
upper-bound the logarithmic primitive, and checks the short-time budget
by sixth powers. It checks exact source reuse for k3. Finally it proves
that the stated rational union upper bound raised to the sixth power
exceeds the fifth power of the sum of certified weighted powers.

The first audit invocation (chunk6f63f8) stopped before case checks: its
sanity assertion used a strict inequality against the exact same decimal
floor of the saved pi lower endpoint. Changing '<' to '<=' corrected
that implementation assertion. The second invocation (chunka4892d)
passed in0.09829420899040997 seconds. No production source or result was
changed, no parameter was added and no scalar configuration was rerun.
This is one fixed postprocessing configuration with one recorded repair,
not two scientific configurations. All eight allowed configurations are
now exhausted and the batch stops.

The saved-data checker is an arithmetic and provenance audit; it does not
replace the mathematical proof or an independent review of that proof.
All pre-existing source hashes, including the F049 audit manifest, were
checked unchanged. The complete dyadic intervals are in results.json;
checked_results.json records the exact rational certificate.

## Files and reproducibility

- registration.json: fixed inputs, source/code hashes, caps and stopping rule.
- theorem.md: proof of all five types and union interface.
- evaluate.py: one-shot producer; refuses an existing run receipt.
- results.json: complete dyadic endpoints and all counts.
- run_receipt.json: actual producer execution and source-result binding.
- check_saved.py: rational-only saved-data verifier.
- checked_results.json: verified rational upper bounds and weights.
- postprocess_receipt.json: the one repaired assertion and final pass.
- manifest.json: final whitelist and unchanged source hashes.

No subsequent numeric batch is authorized here. The next mathematical
object remains the unpaid middle-domain long-time contribution and the
regions outside these short-time middle domains; this file does not
supply their costs. The parent retains final proof and integration
responsibility.
