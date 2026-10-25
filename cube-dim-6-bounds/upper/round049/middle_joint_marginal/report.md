# A049: true two-coordinate L1 sharpening on the043 middle tube

The new joint-marginal interface is valid and strictly improves the
certified short-time budget on exactly the same representative domain

    D={x=(1,y), y in[0,1]^5: dist(x,span(1_A,1_B))<1/64,
                                      mean_B(x) in[1/4,3/4]},
    A={1,2,3}, B={4,5,6}.

For the actual original body32 density,

    ||integral_0^(1/1000) p_s ds||_L^(6/5)(D)
      <66002272277/250000000000000000
      =2.64009089108e-7.

The saved043 short-budget formula was about4.92508528870e-6. The new
upper-budget formula is about0.0536049781 of it, an improvement factor
about18.6549838. This compares two certified upper budgets on the same
domain and time range; it does not assert that the actual norm decreased
or provide a lower bound on that norm. The old formula was not recomputed.

No larger time interval, spatial endpoint, other plane, whole-face I,
global g or new C6 bound is certified. C6 remains64.

## True marginal, symmetry and domain

Fix the single coordinate j=4. Its two sign blocks contain16 independent
Cauchy(s) parameters each, with sums U,V that are independent Cauchy(16s).
Since (X1,Xj)=(U+V,U-V), the inverse Jacobian is1/2 and

    p2_s(1,z)=2048s^2/
      [pi^2(((1+z)^2+1024s^2)((1-z)^2+1024s^2))].

The remaining four coordinate reflections leave both retained coordinates
fixed and preserve the full body32 density. Their16 positive orthants
therefore give a factor exactly1/16. There is no additional factor for
X1 or Xj and no sum over the three possible choices of j. The theorem
also supplies a direct positive basis/Tonelli proof of the marginal
slice identity at every fixed(1,z), rather than relying on an a.e. slice.

The B residual q=x_B-mean_B(x)1_B has sum zero and norm less thanrho.
The coordinate functional on this zero-sum plane has norm sqrt(2/3), so
the true domain projects inside

    L=1/4-sqrt(2/3)/64,
    U=3/4+sqrt(2/3)/64.

These are approximately0.237242240923 and0.762757759077, strictly within
(0,1). Extending only the other four nonnegative coordinates gives

    integral_D p_s dy <=C_joint s^2,
    C_joint=(128/pi^2) integral_L^U(1-z^2)^-2 dz.

The128 is the checked2048/16 normalization. The primitive

    z/[2(1-z^2)]+log((1+z)/(1-z))/4

was differentiated exactly. Evaluation and an independent rational
positive-series check certify

    C_joint<1513131638048143112947/10^20,

with display value15.13131638048143113. On the entire registered time
range, this new L1 budget is smaller than the loose s/pi budget used in
043; at epsilon their ratio is approximately0.04753643238.

## Actual branch accounting and the new time powers

The actual043 branches P_s and B_s are nonnegative, disjoint in parameter
space and sum to p_s. Hence each is dominated by p_s on D. The1/16
reflection factor is applied only to the full density; it is never
asserted for either branch. Enlarged events used in043 to bound the
complementary density do not replace the actual disjoint complement in
this L1 argument.

Read the saved043 pointwise coefficients D_big and M_bad from their
complete lower/upper dyadic endpoints. They still give

    P_s<=D_big s^-2, B_s<=M_bad.

The new L1 time power is s^2. Interpolation at q=6/5 consequently gives
density-norm powers s^(4/3) and s^(5/3). Integration yields exactly

    (3/7)D_big^(1/6)C_joint^(5/6)epsilon^(7/3),
    (3/8)M_bad^(1/6)C_joint^(5/6)epsilon^(8/3).

No pointwise coefficient or event remainder is removed, recalculated or
replaced. The new estimate accounts for all mass through the same complete
043 branch partition.

| New paid contribution | Short L6/5 upper budget, display |
|---|---:|
| Both specified parameters large | 4.924770754499966e-8 |
| Entire actual complementary event | 2.147613815620414e-7 |
| Sum | 2.640090891070410e-7 |

The complement still represents about81.3462% of this upper budget. The
total is about0.00321761077 of the conditional whole-face threshold
4/48750. These comparisons do not establish that threshold's global g
assumption or pay the remaining spatial and time regions.

## Actual new work and verification

The registered work comprised two exact structural configurations and
one fixed Arb configuration, with no parameter scan:

1. An explicit six-column basis verifies determinants32,2,16 for the
   full basis, retained-coordinate pivot and four-variable Schur Jacobian.
   The sign-block counts follow from2^4 free signs; no frame/fan full
   enumeration is performed. The B-coordinate deviation norm squared
   is checked as2/3, and the density normalization simplifies exactly.
2. Symbolic differentiation verifies the primitive, and exact rational
   arithmetic verifies the interpolation exponents and time prefactors.
3. Arb768 evaluates only the fixed L,U,C_joint and the new branch budgets,
   importing043 D_big,M_bad and its old budget from full saved endpoints.
   Their original512-bit uncertainty is preserved. No decimal display
   is used as an input and no old043 formula is evaluated.

The new driver PID49550 completed all three configurations in0.079950
seconds, exit0 (chunk5fc74f). The old043 whitelist and all registered
source hashes were checked before execution. All raw new Arb endpoints
and exact source keys/intervals are preserved in results.json.

The saved-data checker PID49725 passed in0.023209 seconds, exit0
(chunk6f209e). It makes no Arb or kernel call. A fixed96-term positive
atanh series at the same two endpoints, with its explicit geometric
remainder, independently verifies a rational upper bound for C_joint.
Exact sixth-power inequalities then verify the final short-norm upper
bound using the complete saved043 coefficient intervals. This adds no
new geometric or time choice and does not reevaluate an old budget.

There were no failures or retries, no G/H/K/p/native calls, full
enumerations, spatial grids, optimizers or MATLAB calls. All work is
confined to this new directory; old sources remain unchanged. No
delegation, stage, commit or push occurred.

## Frozen interface

The scope remains the same A3 middle tube, rho1/64 and times0<s<=1/1000.
The remaining long-time, endpoint, other-plane and global integration
objects are outside this batch. There is no claim about their overlap
or a total face budget.

Necessary files: registration.json, theorem.md, evaluate.py, results.json,
run_receipt.json, check_saved.py, checked_results.json, this report and
manifest.json. The driver refuses a second batch. The final manifest
records hashes and wall time; the parent retains final proof review and
integration responsibility.
