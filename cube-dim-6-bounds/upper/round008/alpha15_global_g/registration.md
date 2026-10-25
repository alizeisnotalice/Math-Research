# A008 alpha=3/2: finite analytic whole-cube g batch

Continue the existing equal-weight body32 alpha=3/2 model, not the alpha=1
kernel. Only this owned directory may change. The completed chamber
counterexample and all round007 files remain untouched. No MATLAB, high
dimensional quadrature, MC, kernel grid, agents, commits or pushes.

Read-only inputs: C008 alpha15_coercivity report/results; F008 alpha15_product
report/results; the existing alpha15 general-L1/all-scale interface. Hash them.

Exactly four analytic schemes:

1. Degree-2 global cosine lower polynomial; positive spectral integral uses
   psi<=32 r^(3/2), negative integral uses psi>=c6 r^(3/2).
2. The same degree-2 polynomial, but improve p_t(0)'s positive lower term by
   spherical Jensen with the exact mean mu of psi on S5.
3. Degree-6 global cosine lower polynomial, the same Jensen zero-order term,
   positive higher terms using 32 and negative terms using c6.
4. Degree-10 analogue of scheme 3.

At most eight scalar budget scenarios: four optimized whole-cube g values and
their four associated fixed F-norm cost budgets. Common Gamma evaluations,
root bracketing, and symbolic algebra are certificate checks, not extra
parameter searches. Tail optimization is one-dimensional, with strict Arb384
brackets and monotonicity checks. No unregistered high-order expansion.

The whole-cube radius is Euclidean sqrt(6), never 1. All-time Fourier and
positive-time-tail legality must be proved, including x=0 and possible Green
singularities. If current bounds cannot pay 64, report failure of the allocated
certificate budget, not a lower bound on the actual operator/model cost.
